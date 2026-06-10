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
    std::cout << "  CFL = " << getCFL() << std::endl;
    std::cout << "  Diffusion number = " << getDiffusionNumber() << std::endl;
    
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
    if (timeScheme == TimeScheme::EULER_EXPLICIT) {
        // Explicit Euler: u^(n+1) = u^n + dt * RHS(u^n)
        std::vector<double> rhs = computeRHS(u);
        
        for (int i = 0; i < nx; ++i) {
            u[i] += dt * rhs[i];
        }
        
        applyBoundaryConditions(u);
        
    } else if (timeScheme == TimeScheme::RK4) {
        // 4th order Runge-Kutta
        k1 = computeRHS(u);
        
        for (int i = 0; i < nx; ++i) {
            u_temp[i] = u[i] + 0.5 * dt * k1[i];
        }
        applyBoundaryConditions(u_temp);
        k2 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            u_temp[i] = u[i] + 0.5 * dt * k2[i];
        }
        applyBoundaryConditions(u_temp);
        k3 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            u_temp[i] = u[i] + dt * k3[i];
        }
        applyBoundaryConditions(u_temp);
        k4 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
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
        // Ux(di) returns neighbor values handling all boundary cases via ghost cells.
        // See end of file for derivation details.
        auto Ux = [&](int di) -> double {
            int ni = i + di;
            if (ni >= 0 && ni < nx) return u_current[ni];
            if (bcType == BoundaryCondition::PERIODIC)
                return u_current[(ni + nx) % nx];
            if (bcType == BoundaryCondition::DIRICHLET)
                return (ni < 0) ? 2*bcLeft  - u_current[-ni]
                                : 2*bcRight - u_current[2*(nx-1) - ni];
            // NEUMANN
            return (ni < 0) ? u_current[-ni]           + ni * 2 * bcLeft  * dx
                            : u_current[2*(nx-1) - ni] + (ni - (nx-1)) * 2 * bcRight * dx;
        };
        const double ui = u_current[i];
        // RHS = -c * ∂u/∂x + D * ∂²u/∂x²
        rhs[i] = -advectionTerm(Ux, ui) + diffusionTerm(Ux, ui);
    }
    
    return rhs;
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

double AdvectionDiffusionSolver1D::getCFL() const {
    return (std::abs(c) * dt) / dx;  // CFL number 
}

double AdvectionDiffusionSolver1D::getDiffusionNumber() const {
    return (D * dt) / (dx * dx);  // Diffusion number 
}

/*
    Ux(di) return neighbor values handling all boundary cases via ghost cells.
    If inside the domain they return the value directly (e.g. Ux(+1) returns u[idx(i+1,j)]).
    This keeps the scheme equations clean and readable while handling boundary cases
    separately according to the boundary condition type.

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

    The lambda function Ux is defined inside the inner loop so that the current
    index i is captured by reference [&] rather than passed as an explicit argument.
    This keeps the scheme equations clean (e.g. Ux(-1), Ux(+1) instead of Ux(u, i, -1)).
    Redefining the lambdas at each iteration has no performance cost: they are stack-allocated
    objects with no heap allocation, and the compiler inlines them completely under -O2/-O3.
*/
