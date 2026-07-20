#include "solver/AdvectionDiffusionSolver1D.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

AdvectionDiffusionSolver1D::AdvectionDiffusionSolver1D(const ConfigParser& cfg) : BaseSolver(cfg), outputWriter(cfg)
{
    nx          = m_cfg.getInt("GRID_POINTS");
    L           = m_cfg.getDouble("DOMAIN_LENGHT");
    c           = m_cfg.getDouble("ADVECTION_SPEED");
    mesh_file   = m_cfg.getString("MESH_FILE");

    bcLeft      = m_cfg.getDouble("BC_LEFT");
    bcRight     = m_cfg.getDouble("BC_RIGHT");

    initialCondition = makeIC1D(m_cfg);

    // Grid initialization
    dx = (bcType == BoundaryCondition::PERIODIC) ? L / nx : L / (nx - 1);
    x.resize(nx);
    u.resize(nx, 0.0);

    if (timeScheme == TimeScheme::RK4){
        u_temp.resize(nx, 0.0);

        k1.resize(nx, 0.0);
        k2.resize(nx, 0.0);
        k3.resize(nx, 0.0);
        k4.resize(nx, 0.0);
    }

    std::ofstream file(output_dir + mesh_file + ".csv");
    if (!file){
        throw CannotOpenFile(output_dir + mesh_file,"mesh");
    }
    // Handle the "fencepost problem": write the first element outside the loop
    // to ensure the comma acts only as a separator between elements.    
    file << "x0";
    for (int i = 1; i < nx; ++i){
        file << ",x" << i;
    }
    file << '\n';

    if (verbose) std::cout << "Printing mesh grid:" << '\n' << '\n';
    
    for (int i = 0; i < nx - 1; ++i) {
        x[i] = i * dx;
        file << x[i] << ",";
        if (verbose) std::cout << "coordinate x" << i << ": "<< x[i] << '\n';
    }
    x[nx - 1] = (nx - 1) * dx;
    file << x[nx - 1];
    if (verbose) std::cout << "coordinate x" << nx - 1 << ": "<< x[nx - 1] << '\n';

    file.close();

    std::cout << "Advection-Diffusion Solver 1D Initialized" << std::endl;
}

void AdvectionDiffusionSolver1D::setInitialCondition() {

    initialCondition->setIC(u,x);
    std::cout << "Initial condition set." << std::endl;

    // Checking whether the CFL or the diffusion number are too high (for explicit schemes)
    checkStability();
}

void AdvectionDiffusionSolver1D::setBoundaryConditions()
{
    
    std::cout << "Boundary conditions set: ";
    if (bcType == BoundaryCondition::DIRICHLET) {
        std::cout << "Dirichlet (u[0]=" << bcLeft << ", u[N]=" << bcRight << ")";
    } else if (bcType == BoundaryCondition::NEUMANN) {
        std::cout << "Neumann";
    } else {
        std::cout << "Periodic";
    }
    std::cout << std::endl;
}

void AdvectionDiffusionSolver1D::solve() {
    std::cout << "\nStarting solver:" << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    std::cout << "  t_end = " << tEnd << std::endl;
    {
        const StabilityNumbers s = getStabilityNumbers();
        std::cout << "  CFL = " << s.Cx << std::endl;
        std::cout << "  Diffusion number = " << s.dx << std::endl;
        std::cout << "  Stability margin = " << getStabilityMargin() << " (< 1 = stable)" << std::endl;
    }
    
    int nSteps = 0;
    while (t < tEnd) {
        double dt_actual = std::min(dt, tEnd - t);
        step(dt_actual);
        nSteps++;
        
        if (nSteps % 100 == 0) {
            std::cout << "  Step " << nSteps << ", t = " << t << std::endl;
        }
        if (nSteps % output_freq == 0){
            outputWriter.saveCurrentTimeStep(nSteps,u);
        }
    }
    
    std::cout << "Resolution completed after " << nSteps << " time steps." << std::endl;
}

void AdvectionDiffusionSolver1D::step(double dt) {

    // Setting the start and end points for the loops. The boundary points are
    // handled by applyBoundaryConditions() (for Dirichlet and Neumann BCs), and
    // computeRHS leaves rhs = 0 there anyway, so we only iterate on the internal
    // points. For PERIODIC BCs every point is internal, so the loops cover all points.
    int start = (bcType == BoundaryCondition::PERIODIC) ? 0 : 1;
    int end   = (bcType == BoundaryCondition::PERIODIC) ? nx : nx - 1;

    if (timeScheme == TimeScheme::EULER_EXPLICIT) {
        // Explicit Euler: u^(n+1) = u^n + dt * RHS(u^n)
        std::vector<double> rhs = computeRHS(u);

        for (int i = start; i < end; ++i) {
            u[i] += dt * rhs[i];
        }
        
        applyBoundaryConditions(u);
        
    } else if (timeScheme == TimeScheme::RK4) {
        // 4th order Runge-Kutta
        k1 = computeRHS(u);
        
        for (int i = start; i < end; ++i) {
            u_temp[i] = u[i] + 0.5 * dt * k1[i];
        }
        applyBoundaryConditions(u_temp);
        k2 = computeRHS(u_temp);
        
        for (int i = start; i < end; ++i) {
            u_temp[i] = u[i] + 0.5 * dt * k2[i];
        }
        applyBoundaryConditions(u_temp);
        k3 = computeRHS(u_temp);
        
        for (int i = start; i < end; ++i) {
            u_temp[i] = u[i] + dt * k3[i];
        }
        applyBoundaryConditions(u_temp);
        k4 = computeRHS(u_temp);
        
        for (int i = start; i < end; ++i) {
            u[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
        
        applyBoundaryConditions(u);
    }
    
    t += dt;
    time_iter += 1;
}

std::vector<double> AdvectionDiffusionSolver1D::computeRHS(
    const std::vector<double>& u_current) const
{
    std::vector<double> rhs(nx, 0.0);
    
    // Calculate RHS for boundary points
    int start = (bcType == BoundaryCondition::PERIODIC) ? 0 : 1;
    int end = (bcType == BoundaryCondition::PERIODIC) ? nx : nx - 1;
    
    // Calculate RHS for internal points
    for (int i = start; i < end; ++i) {
        // RHS = -c * ∂u/∂x + D * ∂²u/∂x²
        rhs[i] = -advectionTerm(u_current, i) + diffusionTerm(u_current, i);
    }

    return rhs;
}

// neighbor() returns neighbor values handling all boundary cases via ghost cells.
// See end of file for derivation details.
double AdvectionDiffusionSolver1D::neighbor(const std::vector<double>& field, int i, int di) const
{
    int ni = i + di;
    if (ni >= 0 && ni < nx) return field[ni];
    if (bcType == BoundaryCondition::PERIODIC)
        return field[(ni + nx) % nx];
    if (bcType == BoundaryCondition::DIRICHLET)
        return (ni < 0) ? 2*bcLeft  - field[-ni]
                        : 2*bcRight - field[2*(nx-1) - ni];
    // NEUMANN
    return (ni < 0) ? field[-ni]           + ni * 2 * bcLeft  * dx
                    : field[2*(nx-1) - ni] + (ni - (nx-1)) * 2 * bcRight * dx;
}

double AdvectionDiffusionSolver1D::advectionTerm(const std::vector<double>& field, int i) const
{
    if (spatialScheme == SpatialScheme::CENTRAL) {
        return centralDifference(field, i);
    } else if (spatialScheme == SpatialScheme::UPWIND) {
        return upwindDifference(field, i);
    } else { // QUICK
        return quickDifference(field, i);
    }
}

double AdvectionDiffusionSolver1D::diffusionTerm(const std::vector<double>& field, int i) const
{
    // Lambda function to improve readability of the equation below
    auto Ux = [&](int di) { return neighbor(field, i, di); };
    // Central differencies of the second order for the diffusion.
    return D * (Ux(+1) - 2.0*Ux(0) + Ux(-1)) / (dx * dx);
}

double AdvectionDiffusionSolver1D::centralDifference(const std::vector<double>& field, int i) const
{
    auto Ux = [&](int di) { return neighbor(field, i, di); };
    return c * (Ux(+1) - Ux(-1)) / (2.0 * dx);
}

double AdvectionDiffusionSolver1D::upwindDifference(const std::vector<double>& field, int i) const
{
    auto Ux = [&](int di) { return neighbor(field, i, di); };
    return c * ((c >= 0) ? (Ux(0) - Ux(-1)) / dx : (Ux(+1) - Ux(0)) / dx);
}

double AdvectionDiffusionSolver1D::quickDifference(const std::vector<double>& field, int i) const
{
    // QUICK Scheme (Quadratic Upstream Interpolation for Convective Kinematics)

    auto Ux = [&](int di) { return neighbor(field, i, di); };
    double dux = (c >= 0)
        ? (Ux(-2) - 7*Ux(-1) + 3*Ux(0) + 3*Ux(+1))  / (8.0 * dx)
        : ( -3*Ux(-1) - 3*Ux(0) +7*Ux(+1) - Ux(+2)) / (8.0 * dx);

    return c * dux;
}

void AdvectionDiffusionSolver1D::applyBoundaryConditions(std::vector<double>& u_vec) {
    if (bcType == BoundaryCondition::DIRICHLET) {
        u_vec[0] = bcLeft;
        u_vec[nx - 1] = bcRight;
        
    } else if (bcType == BoundaryCondition::NEUMANN) {
        // du/dx = bcLeft at the left boundary
        u_vec[0] = u_vec[1] - bcLeft * dx;
        // du/dx = bcRight at the right boundary
        u_vec[nx - 1] = u_vec[nx - 2] + bcRight * dx;
    } 
}

BaseSolver::StabilityNumbers AdvectionDiffusionSolver1D::getStabilityNumbers() const {
    StabilityNumbers s;
    s.Cx = std::abs(c) * dt / dx;
    s.dx = D * dt / (dx * dx);
    return s;
}

/*
    neighbor(di) returns neighbor values handling all boundary cases via ghost cells.
    If inside the domain it returns the value directly (e.g. neighbor with di = +1
    returns u[i+1]). This keeps the scheme equations clean and readable while handling
    boundary cases separately according to the boundary condition type.
    Inside each spatial scheme, a thin local lambda Ux(di) aliases this function so the
    stencil formulas stay compact (Ux(+1) instead of neighbor(field, i, +1)).

    PERIODIC BCs Take the value on the other side
    (ni + nx) % nx es. ni = -1 => (-1 + nx)/nx = nx - 1
                       ni = -1 --------------->  nx - 1
                   es. ni = nx => (nx + nx)/nx = 0
                       ni = nx --------------->  0

    Ghost points for DIRICHLET BCs
    Bc_Diriclet = (u(i+1) + u(i-1)) / 2

    i = 0 => bcLeft         u[-1] = 2*bcLeft - u[1] 
    ni = -1 → 2*bcLeft - u[-(-1)] = 2*bcLeft - u[1] 
    ni = -2 → 2*bcLeft - u[-(-2)] = 2*bcLeft - u[2] 

    i = nx - 1 => bcRight               u[nx]  = 2*bcRight - u[nx - 2]
    ni = nx    → 2*bcRight - u[2*nx-2-nx]      = 2*bcRight - u[nx - 2]
    ni = nx+1  → 2*bcRight - u[2*nx-2-(nx+1)]  = 2*bcRight - u[nx - 3]

    Ghost points for NEUMANN BCs
    Bc_Neumann = (u(i+1) - u(i-1)) / (2 * dx)
    Bc_Neumann = (u(i+2) - u(i-2)) / (4 * dx)

    i = 0 => bcLeft                  u[-1] = u[1] - 2 * bcLeft * dx
    ni = -1 → u[-(-1)] + (-1)*2*bcLeft*dx  = u[1] - 2 * bcLeft * dx 
    ni = -2 → u[-(-2)] + (-2)*2*bcLeft*dx  = u[2] - 4 * bcLeft * dx 

    i = nx - 1 => bcRight                      u[nx] = u[nx-2] + 2 * bcRight * dx 
    ni = nx   → u[2*nx-2-nx] + (1)*2*bcRight*dx      = u[nx-2] + 2 * bcRight * dx 
    ni = nx+1 → u[2*nx-2-(nx+1)] + (2)*2*bcRight*dx  = u[nx-3] + 4 * bcRight * dx 

    Performance note: neither the member-function calls nor the adapter lambda Ux add
    any overhead. Everything lives in the same translation unit, so under -O2/-O3 the compiler
    inlines both layers completely — the generated machine code is the same as writing the
    ghost-cell logic directly inside the stencil formulas.
*/
