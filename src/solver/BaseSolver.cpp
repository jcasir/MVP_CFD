#include "solver/BaseSolver.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <filesystem>

BaseSolver::BaseSolver(const ConfigParser& config) : m_cfg(config)
{
    // Read the config file

    output_freq             = m_cfg.getInt("OUTPUT_FREQUENCY");
    time_iter               = 0;

    D                       = m_cfg.getDouble("DIFFUSION_COEFF");
    dt                      = m_cfg.getDouble("TIME_STEP");
    tEnd                    = m_cfg.getDouble("END_TIME");
    t                       = m_cfg.getDouble("INITIAL_TIME");

    spatialScheme           = m_cfg.getEnum("SPATIAL_SCHEME",spatialSchemeMap);
    timeScheme              = m_cfg.getEnum("TIME_SCHEME",timeSchemeMap);
    bcType                  = m_cfg.getEnum("BOUNDARY_CONDITIONS",boundaryConditionMap);

    verbose                 = m_cfg.getBool("VERBOSE");

    output_file             = m_cfg.getString("OUTPUT_FILE");
    output_dir              = m_cfg.getString("OUTPUT_DIR");
    // Accept OUTPUT_DIR with or without a trailing '/'
    if (!output_dir.empty() && output_dir.back() != '/') output_dir += '/';

    // Create directory for the results if it doesn't exists. Makefile already creates the results/ dir,
    // but if the config includes additional subdirectiories they are not present when the program runs.
    // This section ensures that the folder is always present in every case.
    std::filesystem::path output_dir_path = output_dir;
    if (std::filesystem::create_directories(output_dir_path)) {
        std::cout << "Results directory succesfully created: " << output_dir_path << '\n';
    } else {
        std::cout << "Results directory already exists: " << output_dir_path << '\n';
    }


    std::cout << "Base parameters of the solver initialized" << std::endl;
    
    if (verbose) m_cfg.print();
}


bool BaseSolver::getExplicitSchemeFlag() const {
    return (timeScheme == TimeScheme::EULER_EXPLICIT || 
            timeScheme == TimeScheme::RK4);
}

void BaseSolver::checkStability(bool report) const {

    if (report) std::cout << "\nVerifying stability conditions (von Neumann)\n";

    if (!getExplicitSchemeFlag()) {
        if (report) std::cout << "Implicit time scheme: no explicit stability limits.\n";
        return;
    }

    // dx and dy here are the diffusion numbers in the x and y directions, not the grid spacing.
    const StabilityNumbers s = getStabilityNumbers();
    const double d_tot = s.dx + s.dy;
    const double C_tot = s.Cx + s.Cy;

    // NOTE: the scheme-specific limits below are derived for explicit Euler.
    if (timeScheme == TimeScheme::RK4) {
        // RK4 has a much larger stability region that includes a wide stretch of the
        // imaginary axis (|z| <= 2.83): central/QUICK advection is stable even with
        // D = 0, and the Euler-derived conditions below would wrongly reject valid
        // setups (e.g. pure advection with QUICK + RK4). Safe practical bounds:
        if (report) std::cout << "  RK4:  CFL = " << C_tot << " (limit 1.0), "
                              << "d = " << d_tot << " (limit 0.5)\n";
        if (C_tot > 1.0)
            throw StabilityException("RK4 CFL", C_tot, 1.0);
        if (d_tot > 0.5)
            throw StabilityException("RK4 diffusion number", d_tot, 0.5);
        if (report) std::cout << "Stability conditions satisfied.\n\n";
        return;
    }

    if (spatialScheme == SpatialScheme::UPWIND) {
        // Necessary and sufficient (verified numerically): advection and diffusion
        // share the same stability budget.
        const double budget = C_tot + 2.0 * d_tot;
        if (report) std::cout << "  UPWIND:  C + 2d = " << budget << " (limit 1.0)\n";
        if (budget > 1.0)
            throw StabilityException("Upwind budget C+2d", budget, 1.0);
    }
    else { // CENTRAL and QUICK share the low-frequency coupled condition
        // 1) diffusive / high-frequency limit
        if (spatialScheme == SpatialScheme::CENTRAL) {
            if (report) std::cout << "  Diffusion number = " << d_tot << " (limit 0.5)\n";
            if (d_tot > 0.5)
                throw StabilityException("Diffusion number", d_tot, 0.5);
        } else { // QUICK: its own dissipation adds to the physical one at theta = pi
            const double hf = C_tot + 4.0 * d_tot;
            if (report) std::cout << "  QUICK high-freq:  C + 4d = " << hf << " (limit 2.0)\n";
            if (hf > 2.0)
                throw StabilityException("QUICK budget C+4d", hf, 2.0);
        }

        // 2) coupled low-frequency condition: Cx^2/dx + Cy^2/dy <= 2.
        //    Centered (and QUICK) advection with no diffusion is unconditionally
        //    unstable with explicit Euler: physical diffusion is what stabilizes it.
        if (d_tot <= 0.0) {
            if (C_tot > 0.0)
                throw StabilityException(
                    "Coupled condition: CENTRAL/QUICK advection with D = 0 "
                    "is unconditionally unstable with explicit Euler", C_tot, 0.0);
        } else {
            double coupled = 0.0;
            if (s.Cx > 0.0) coupled += (s.dx > 0.0) ? s.Cx*s.Cx/s.dx
                                                    : std::numeric_limits<double>::infinity();
            if (s.Cy > 0.0) coupled += (s.dy > 0.0) ? s.Cy*s.Cy/s.dy
                                                    : std::numeric_limits<double>::infinity();
            if (report) std::cout << "  Coupled:  Cx^2/dx + Cy^2/dy = " << coupled << " (limit 2.0)\n";
            if (coupled > 2.0)
                throw StabilityException("Coupled advection-diffusion", coupled, 2.0);
        }

        // 3) cell Peclet: not a stability limit, a spatial-quality warning (CENTRAL only)
        if (report && spatialScheme == SpatialScheme::CENTRAL) {
            const double pe_x = (s.dx > 0.0) ? s.Cx / s.dx : 0.0;  // = c*dx/D
            const double pe_y = (s.dy > 0.0) ? s.Cy / s.dy : 0.0;
            const double pe   = std::max(pe_x, pe_y);
            if (pe > 2.0)
                std::cerr << "WARNING: cell Peclet = " << pe << " > 2 with CENTRAL: "
                          << "the steady solution may show node-to-node oscillations. "
                          << "Consider UPWIND/QUICK or a finer grid.\n";
        }
    }

    if (report) std::cout << "Stability conditions satisfied.\n\n";
}

double BaseSolver::getStabilityMargin() const {
    // Same conditions as checkStability(), expressed as the worst ratio value/limit.
    const StabilityNumbers s = getStabilityNumbers();
    const double d_tot = s.dx + s.dy;
    const double C_tot = s.Cx + s.Cy;

    if (timeScheme == TimeScheme::RK4)
        return std::max(C_tot, d_tot / 0.5);

    if (spatialScheme == SpatialScheme::UPWIND)
        return C_tot + 2.0 * d_tot;                       // limit 1

    double margin = (spatialScheme == SpatialScheme::CENTRAL)
                  ? d_tot / 0.5
                  : (C_tot + 4.0 * d_tot) / 2.0;          // QUICK high-frequency

    double coupled = 0.0;
    if (s.Cx > 0.0) coupled += (s.dx > 0.0) ? s.Cx*s.Cx/s.dx
                                            : std::numeric_limits<double>::infinity();
    if (s.Cy > 0.0) coupled += (s.dy > 0.0) ? s.Cy*s.Cy/s.dy
                                            : std::numeric_limits<double>::infinity();
    return std::max(margin, coupled / 2.0);
}