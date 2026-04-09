#include "AdvectionDiffusionSolver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

AdvectionDiffusionSolver::AdvectionDiffusionSolver(const std::string& filename)
    // int nx, double L, double c, double D,
    // SpatialScheme spatialScheme, TimeScheme timeScheme)
    // : nx(nx), L(L), c(c), D(D), 
    //   spatialScheme(spatialScheme), timeScheme(timeScheme),
    //   t(0.0), bcType(BoundaryCondition::DIRICHLET), bcLeft(0.0), bcRight(0.0)
{
    // Lettura config
    ConfigParser cfg(filename);

    const int nx                                = cfg.getInt("GRID_POINTS");
    const int L                                 = cfg.getInt("DOMAIN_LENGHT");
    const double c                              = cfg.getDouble("ADVECTION_SPEED");
    const double D                              = cfg.getDouble("DIFFUSION_SPEED");
    const double bcLeft                         = cfg.getDouble("BC_LEFT",0.0);
    const double bcRight                        = cfg.getDouble("BC_RIGHT",0.0);
    const double dt                             = cfg.getDouble("TIME_STEP");
    const double tEnd                           = cfg.getDouble("END_TIME");
    const double t                              = cfg.getDouble("INITIAL_TIME",0.0);
    const SpatialScheme spatialScheme           = cfg.getEnum("SPATIAL_SCHEME",spatialSchemeMap);
    const TimeScheme timeScheme                 = cfg.getEnum("TIME_SCHEME",timeSchemeMap);
    const BoundaryCondition bcType              = cfg.getEnum("BOUNDARY_CONDITIONS",boundaryConditionMap);
    const InitialCondition initialCondition     = cfg.getEnum("INITIAL_CONDITIONS",initialConditionMap);

    const bool verbose                          = cfg.getBool  ("VERBOSE");

    // Grid initialization
    dx = L / (nx - 1);
    x.resize(nx);
    u.resize(nx, 0.0);
    
    for (int i = 0; i < nx; ++i) {
        x[i] = i * dx;
    }

    std::cout << "Advection-Diffusion Solver 1D Initialised" << std::endl;
    
    if (verbose) cfg.print();

    // Setting initial conditions
    solver.setInitialCondition(retriveInitialConditionFunction(initialCondition));

    // Setting boundary conditions
    solver.setBoundaryConditions();
}

void AdvectionDiffusionSolver::setInitialCondition(std::function<double(double)> ic) {
    for (int i = 0; i < nx; ++i) {
        u[i] = ic(x[i]);
    }
    std::cout << "Condizione iniziale impostata." << std::endl;
}

void AdvectionDiffusionSolver::setBoundaryConditions()
{
    
    std::cout << "Condizioni al contorno impostate: ";
    if (bcType == BoundaryCondition::DIRICHLET) {
        std::cout << "Dirichlet (u[0]=" << bcLeft << ", u[N]=" << bcRight << ")";
    } else if (bcType == BoundaryCondition::NEUMANN) {
        std::cout << "Neumann";
    } else {
        std::cout << "Periodiche";
    }
    std::cout << std::endl;
}

void AdvectionDiffusionSolver::solve() {
    std::cout << "\nInizio risoluzione:" << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    std::cout << "  t_end = " << tEnd << std::endl;
    std::cout << "  CFL = " << getCFL() << std::endl;
    std::cout << "  Numero diffusione = " << getDiffusionNumber() << std::endl;
    
    int nSteps = 0;
    while (t < tEnd) {
        double dt_actual = std::min(dt, tEnd - t);
        step(dt_actual);
        nSteps++;
        
        if (nSteps % 100 == 0) {
            std::cout << "  Step " << nSteps << ", t = " << t << std::endl;
        }
    }
    
    std::cout << "Risoluzione completata dopo " << nSteps << " passi temporali." << std::endl;
}

void AdvectionDiffusionSolver::step(double dt) {
    if (timeScheme == TimeScheme::EULER_EXPLICIT) {
        // Eulero esplicito: u^(n+1) = u^n + dt * RHS(u^n)
        std::vector<double> rhs = computeRHS(u);
        
        for (int i = 0; i < nx; ++i) {
            u[i] += dt * rhs[i];
        }
        
        applyBoundaryConditions(u);
        
    } else if (timeScheme == TimeScheme::RK4) {
        // Runge-Kutta 4° order
        std::vector<double> k1 = computeRHS(u);
        
        std::vector<double> u_temp(nx);
        for (int i = 0; i < nx; ++i) {
            u_temp[i] = u[i] + 0.5 * dt * k1[i];
        }
        applyBoundaryConditions(u_temp);
        std::vector<double> k2 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            u_temp[i] = u[i] + 0.5 * dt * k2[i];
        }
        applyBoundaryConditions(u_temp);
        std::vector<double> k3 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            u_temp[i] = u[i] + dt * k3[i];
        }
        applyBoundaryConditions(u_temp);
        std::vector<double> k4 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            u[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }
        
        applyBoundaryConditions(u);
    }
    
    t += dt;
}

std::vector<double> AdvectionDiffusionSolver::computeRHS(
    const std::vector<double>& u_current) const
{
    std::vector<double> rhs(nx, 0.0);
    
    // Calcola RHS per i punti interni
    int start = (bcType == BoundaryCondition::PERIODIC) ? 0 : 1;
    int end = (bcType == BoundaryCondition::PERIODIC) ? nx : nx - 1;
    
    for (int i = start; i < end; ++i) {
        // RHS = -c * ∂u/∂x + D * ∂²u/∂x²
        rhs[i] = -advectionTerm(u_current, i) + diffusionTerm(u_current, i);
    }
    
    return rhs;
}

double AdvectionDiffusionSolver::advectionTerm(
    const std::vector<double>& u_current, int i) const
{
    if (spatialScheme == SpatialScheme::CENTRAL) {
        return c * centralDifference(u_current, i);
    } else if (spatialScheme == SpatialScheme::UPWIND) {
        return c * upwindDifference(u_current, i);
    } else { // QUICK
        return c * quickDifference(u_current, i);
    }
}

double AdvectionDiffusionSolver::diffusionTerm(
    const std::vector<double>& u_current, int i) const
{
    // Differenze centrali del secondo ordine per la diffusione
    int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : i - 1;
    int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : i + 1;
    
    return D * (u_current[ip1] - 2.0*u_current[i] + u_current[im1]) / (dx * dx);
}

double AdvectionDiffusionSolver::centralDifference(
    const std::vector<double>& u_current, int i) const
{
    int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : i - 1;
    int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : i + 1;
    
    return (u_current[ip1] - u_current[im1]) / (2.0 * dx);
}

double AdvectionDiffusionSolver::upwindDifference(
    const std::vector<double>& u_current, int i) const
{
    if (c > 0) {
        // Upwind in direzione negativa
        int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : i - 1;
        return (u_current[i] - u_current[im1]) / dx;
    } else {
        // Upwind in direzione positiva
        int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : i + 1;
        return (u_current[ip1] - u_current[i]) / dx;
    }
}

double AdvectionDiffusionSolver::quickDifference(
    const std::vector<double>& u_current, int i) const
{
    // Schema QUICK (QUadratic Upstream Interpolation for Convective Kinematics)
    // Richiede 3 punti upstream
    
    if (c > 0) {
        int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : std::max(0, i - 1);
        int im2 = (im1 == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : std::max(0, im1 - 1);
        int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : std::min(nx - 1, i + 1);
        
        // Interpolazione quadratica
        return (-u_current[im2] + 8.0*u_current[ip1] - 8.0*u_current[im1] + u_current[i]) / (12.0 * dx);
    } else {
        int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : std::min(nx - 1, i + 1);
        int ip2 = (ip1 == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : std::min(nx - 1, ip1 + 1);
        int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : std::max(0, i - 1);
        
        return (-u_current[i] + 8.0*u_current[ip1] - 8.0*u_current[im1] + u_current[ip2]) / (12.0 * dx);
    }
}

void AdvectionDiffusionSolver::applyBoundaryConditions(std::vector<double>& u_vec) {
    if (bcType == BoundaryCondition::DIRICHLET) {
        u_vec[0] = bcLeft;
        u_vec[nx - 1] = bcRight;
        
    } else if (bcType == BoundaryCondition::NEUMANN) {
        // du/dx = bcLeft al contorno sinistro
        u_vec[0] = u_vec[1] - bcLeft * dx;
        // du/dx = bcRight al contorno destro
        u_vec[nx - 1] = u_vec[nx - 2] + bcRight * dx;
        
    } else if (bcType == BoundaryCondition::PERIODIC) {
        // Le condizioni periodiche sono già gestite negli indici
        // Ma forziamo l'uguaglianza per sicurezza
        double avg = 0.5 * (u_vec[0] + u_vec[nx - 1]);
        u_vec[0] = avg;
        u_vec[nx - 1] = avg;
    }
}

double AdvectionDiffusionSolver::getCFL() const {
    return std::abs(c) * dx;  // CFL number (moltiplicare per dt per ottenere condizione)
}

double AdvectionDiffusionSolver::getDiffusionNumber() const {
    return D / (dx * dx);  // Diffusion number (moltiplicare per dt)
}

void AdvectionDiffusionSolver::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
        return;
    }
    
    file << "# x u(x,t=" << t << ")" << std::endl;
    file << std::scientific << std::setprecision(10);
    
    for (int i = 0; i < nx; ++i) {
        file << x[i] << " " << u[i] << std::endl;
    }
    
    file.close();
    std::cout << "Soluzione salvata in " << filename << std::endl;
}

void AdvectionDiffusionSolver::printStats() const {
    double uMin = *std::min_element(u.begin(), u.end());
    double uMax = *std::max_element(u.begin(), u.end());
    double uMean = 0.0;
    for (double val : u) {
        uMean += val;
    }
    uMean /= nx;
    
    std::cout << "\nStatistiche soluzione (t = " << t << "):" << std::endl;
    std::cout << "  u_min  = " << uMin << std::endl;
    std::cout << "  u_max  = " << uMax << std::endl;
    std::cout << "  u_mean = " << uMean << std::endl;
}
