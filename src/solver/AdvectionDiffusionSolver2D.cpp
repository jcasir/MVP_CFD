#include "solver/AdvectionDiffusionSolver2D.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

AdvectionDiffusionSolver2D::AdvectionDiffusionSolver2D(const ConfigParser& config)
 : BaseSolver(config), pvdWriter(output_dir + output_file + ".pvd")
{

    nx          = m_cfg.getInt("GRID_POINTS_X");
    ny          = m_cfg.getInt("GRID_POINTS_Y");
    Lx          = m_cfg.getDouble("DOMAIN_LENGHT_X");
    Ly          = m_cfg.getDouble("DOMAIN_LENGHT_Y");

    c_x         = m_cfg.getDouble("ADVECTION_SPEED_X");
    c_y         = m_cfg.getDouble("ADVECTION_SPEED_Y");

    bcTop       = m_cfg.getDouble("BC_TOP");
    bcBottom    = m_cfg.getDouble("BC_BOTTOM");


    // Grid initialization
    dx = (bcType == BoundaryCondition::PERIODIC) ? Lx / nx : Lx / (nx - 1);
    dy = (bcType == BoundaryCondition::PERIODIC) ? Ly / ny : Ly / (ny - 1);

    x.resize(nx);
    y.resize(ny);

    u.resize(nx * ny, 0.0);

    if (timeScheme == TimeScheme::RK4){
        u_temp.resize(nx * ny, 0.0);

        k1.resize(nx * ny, 0.0);
        k2.resize(nx * ny, 0.0);
        k3.resize(nx * ny, 0.0);
        k4.resize(nx * ny, 0.0);
    }

    initialCondition = makeIC2D(m_cfg);


    if (verbose) std::cout << "Printing mesh grid:" << '\n' << '\n';
    
    for (int i = 0; i < nx; ++i) {
        x[i] = i * dx;
        if (verbose) std::cout << "coordinate x" << i << ": "<< x[i] << '\n';
    }

    for (int i = 0; i < ny; ++i) {
        y[i] = i * dy;
        if (verbose) std::cout << "coordinate y" << i << ": "<< y[i] << '\n';
    }

    std::cout << "Advection-Diffusion Solver 2D Initialised" << std::endl;

}

AdvectionDiffusionSolver2D::~AdvectionDiffusionSolver2D() {
    try {
        pvdWriter.write();
    } catch (const std::exception& e) {
        std::cerr << "ERROR writing PVD: " << e.what() << "\n";
    }
}

int AdvectionDiffusionSolver2D::idx(int i,int j) const{
    return (i * ny + j);
}

void AdvectionDiffusionSolver2D::setInitialCondition() {

    initialCondition->setIC(u,x,y);
    std::cout << "Condizione iniziale impostata." << std::endl;

    // Checking whether the CFL or the diffusion number are too high (for explicit schemes)
    checkStability();

    //print initial condition
    // VTUWriter wants the number of cells so it must be given (nx - 1) because nx is the number of points
    // 1 for nz is the default to set the dimension to 2D.
    std::string outputFile = makeVTUFilename(output_file,0);
    VTUWriter vtuWriter(output_dir + outputFile,(nx - 1),(ny - 1),1,dx,dy,0.0);
    vtuWriter.addScalar("u",u);
    vtuWriter.write();
    pvdWriter.addStep(0.0,outputFile);
}

void AdvectionDiffusionSolver2D::setBoundaryConditions()
{
    
    std::cout << "Boundary conditions set: ";
    if (bcType == BoundaryCondition::DIRICHLET) {
        std::cout << "Dirichlet (left=" << bcLeft << ", right=" << bcRight 
          << ", bottom=" << bcBottom << ", top=" << bcTop << ")";
    } else if (bcType == BoundaryCondition::NEUMANN) {
        std::cout << "Neumann";
    } else {
        std::cout << "Periodic";
    }
    std::cout << std::endl;
}

void AdvectionDiffusionSolver2D::solve() {
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
            std::string outputFile = makeVTUFilename(output_file,nSteps);
            VTUWriter outputWriter(output_dir + outputFile,(nx - 1),(ny - 1),1,dx,dy,0.0);
            outputWriter.addScalar("u",u);
            outputWriter.write();
            pvdWriter.addStep(t,outputFile);
        }
    }
    
    std::cout << "Risoluzione completata dopo " << nSteps << " passi temporali." << std::endl;
}

void AdvectionDiffusionSolver2D::step(double dt) {
    if (timeScheme == TimeScheme::EULER_EXPLICIT) {
        // Esplicit Euler: u^(n+1) = u^n + dt * RHS(u^n)
        std::vector<double> rhs = computeRHS(u);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u[idx(i,j)] += dt * rhs[idx(i,j)];
            }
        }
        
        applyBoundaryConditions(u);
        
    } else if (timeScheme == TimeScheme::RK4) {
        // Runge-Kutta 4° order
        k1 = computeRHS(u);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u_temp[idx(i,j)] = u[idx(i,j)] + 0.5 * dt * k1[idx(i,j)];
            }
        }
        applyBoundaryConditions(u_temp);
        k2 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u_temp[idx(i,j)] = u[idx(i,j)] + 0.5 * dt * k2[idx(i,j)];
            }
        }
        applyBoundaryConditions(u_temp);
        k3 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u_temp[idx(i,j)] = u[idx(i,j)] + dt * k3[idx(i,j)];
            }
        }
        applyBoundaryConditions(u_temp);
        k4 = computeRHS(u_temp);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u[idx(i,j)] += (dt / 6.0) * (k1[idx(i,j)] + 2.0*k2[idx(i,j)] + 2.0*k3[idx(i,j)] + k4[idx(i,j)]);
            }
        }
        
        applyBoundaryConditions(u);
    }
    
    t += dt;
}

std::vector<double> AdvectionDiffusionSolver2D::computeRHS(
    const std::vector<double>& u_current) const
{
    std::vector<double> rhs(nx * ny, 0.0);
    
    // Computes RHS for the internal points
    int start = (bcType == BoundaryCondition::PERIODIC) ? 0 : 1;
    int end_x = (bcType == BoundaryCondition::PERIODIC) ? nx : nx - 1;
    int end_y = (bcType == BoundaryCondition::PERIODIC) ? ny : ny - 1;
    
    for (int i = start; i < end_x; ++i) {
        for (int j = start; j < end_y ; ++j){
            // Ux(di) and Uy(dj) return neighbor values handling all boundary cases via ghost cells.
            // See end of file for derivation details.
            auto Ux = [&](int di) -> double {
                int ni = i + di;
                if (ni >= 0 && ni < nx) return u_current[idx(ni, j)];
                if (bcType == BoundaryCondition::PERIODIC)
                    return u_current[idx((ni + nx) % nx, j)];  
                if (bcType == BoundaryCondition::DIRICHLET)
                    return (ni < 0) ? 2*bcLeft  - u_current[idx(-ni, j)]
                                    : 2*bcRight - u_current[idx(2*(nx-1) - ni, j)];
                // NEUMANN
                return (ni < 0) ? u_current[idx(-ni,           j)] + ni * 2 * bcLeft  * dx
                                : u_current[idx(2*(nx-1) - ni, j)] + (ni - (nx-1)) * 2 * bcRight * dx;
            };
            auto Uy = [&](int dj) -> double {
                int nj = j + dj;
                if (nj >= 0 && nj < ny) return u_current[idx(i, nj)];
                if (bcType == BoundaryCondition::PERIODIC)
                    return u_current[idx(i, (nj + ny) % ny)];
                if (bcType == BoundaryCondition::DIRICHLET)
                    return (nj < 0) ? 2*bcBottom - u_current[idx(i, -nj)]
                                    : 2*bcTop    - u_current[idx(i, 2*(ny-1) - nj)];
                // NEUMANN
                return (nj < 0) ? u_current[idx(i, 0)]    + nj * 2 * bcBottom * dy
                                : u_current[idx(i, ny-1)] + (nj - (ny-1)) * 2 * bcTop * dy;
            };
            const double uij = u_current[idx(i,j)];
            // RHS = -(c_x*∂u/∂x + c_y*∂u/∂y) + D*(∂²u/∂x² + ∂²u/∂y²)
            rhs[idx(i,j)] = -advectionTerm(Ux, Uy, uij) + diffusionTerm(Ux, Uy, uij);
        }
    }
    
    return rhs;
}

void AdvectionDiffusionSolver2D::applyBoundaryConditions(std::vector<double>& u_vec) {
    if (bcType == BoundaryCondition::DIRICHLET) {
        for (int i = 0; i < nx; ++i) {
            u_vec[idx(i, 0)] = bcBottom;
            u_vec[idx(i, ny - 1)] = bcTop;
        }
        for (int j = 0; j < ny; ++j){
            u_vec[idx(0, j)] = bcLeft;
            u_vec[idx(nx - 1,j)] = bcRight;
        }
    } 

    else if (bcType == BoundaryCondition::NEUMANN) {

        for (int i = 0; i < nx; ++i) {
            u_vec[idx(i, 0)] = u_vec[idx(i, 1)] - bcBottom * dy;
            u_vec[idx(i,ny - 1)] = u_vec[idx(i,ny - 2)] + bcTop * dy;
        }
        for (int j = 0; j < ny; ++j){
            u_vec[idx(0, j)] = u_vec[idx(1, j)] - bcLeft * dx;
            u_vec[idx(nx - 1,j)] = u_vec[idx(nx - 2,j)] + bcRight * dx;
        }
    }
}

double AdvectionDiffusionSolver2D::getCFL() const {
    double termX = std::abs(c_x) / dx;
    double termY = std::abs(c_y) / dy;
    return dt * (termX + termY);  // CFL number
}

double AdvectionDiffusionSolver2D::getDiffusionNumber() const {
    return D * dt * (1.0/(dx*dx) + 1.0/(dy*dy)); // Diffusion number 
}

/*
    Ux(di) and Uy(dj) return neighbor values handling all boundary cases via ghost cells.
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

    The lambda functions Ux and Uy are defined inside the inner loop so that the current
    indexes i and j are captured by reference [&] rather than passed as explicit arguments.
    This keeps the scheme equations clean (e.g. Ux(-1), Ux(+1) instead of Ux(u, i, j, -1)).
    Redefining the lambdas at each iteration has no performance cost: they are stack-allocated
    objects with no heap allocation, and the compiler inlines them completely under -O2/-O3.
*/
