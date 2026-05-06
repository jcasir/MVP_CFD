#include "solver/AdvectionDiffusionSolver2D.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

AdvectionDiffusionSolver2D::AdvectionDiffusionSolver2D(const ConfigParser& config) : BaseSolver(config)
{

    ny = m_cfg.getInt("GRID_POINTS_X");
    nx = m_cfg.getInt("GRID_POINTS_Y");
    Lx = m_cfg.getDouble("DOMAIN_LENGHT_X");
    Ly = m_cfg.getDouble("DOMAIN_LENGHT_Y");


    // Grid initialization
    dx = Lx / (nx - 1);
    dy = Ly / (ny - 1);

    x.resize(nx);
    y.resize(ny);

    u.resize(nx * ny, 0.0);

    initialCondition = makeIC2D(m_cfg);

    std::ofstream file(mesh_file);
    
    // Handle the "fencepost problem": write the first element outside the loop
    // to ensure the comma acts only as a separator between elements.    
    file << "x0";
    for (int i = 1; i < nx; ++i){
        file << ",x" << i;
    }
    for (int i = 0; i < ny; ++i){
        file << ",y" << i;
    }
    file << '\n';

    if (verbose) std::cout << "Printing mesh grid:" << '\n' << '\n';
    
    for (int i = 0; i < nx; ++i) {
        x[i] = i * dx;
        file << x[i] << ",";
        if (verbose) std::cout << "coordinate x" << i << ": "<< x[i] << '\n';
    }
    x[nx - 1] = (nx - 1) * dx;
    file << x[nx - 1];
    if (verbose) std::cout << "coordinate x" << nx - 1 << ": "<< x[nx - 1] << '\n';

    for (int i = 0; i < ny; ++i) {
        y[i] = i * dy;
        file << y[i] << ",";
        if (verbose) std::cout << "coordinate y" << i << ": "<< y[i] << '\n';
    }
    y[ny - 1] = (ny - 1) * dy;
    file << y[ny - 1];
    if (verbose) std::cout << "coordinate y" << ny - 1 << ": "<< y[ny - 1] << '\n';


    std::cout << "Advection-Diffusion Solver 2D Initialised" << std::endl;
}

void AdvectionDiffusionSolver2D::setInitialCondition() {

    initialCondition->setIC(u,x,y);
    std::cout << "Condizione iniziale impostata." << std::endl;
}

void AdvectionDiffusionSolver2D::setBoundaryConditions()
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
    }
    
    std::cout << "Risoluzione completata dopo " << nSteps << " passi temporali." << std::endl;
}

void AdvectionDiffusionSolver2D::step(double dt) {
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

std::vector<double> AdvectionDiffusionSolver2D::computeRHS(
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

double AdvectionDiffusionSolver2D::advectionTerm(
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

double AdvectionDiffusionSolver2D::diffusionTerm(
    const std::vector<double>& u_current, int i) const
{
    // Differenze centrali del secondo ordine per la diffusione
    int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : i - 1;
    int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : i + 1;
    
    return D * (u_current[ip1] - 2.0*u_current[i] + u_current[im1]) / (dx * dx);
}

double AdvectionDiffusionSolver2D::centralDifference(
    const std::vector<double>& u_current, int i) const
{
    int im1 = (i == 0 && bcType == BoundaryCondition::PERIODIC) ? nx - 1 : i - 1;
    int ip1 = (i == nx - 1 && bcType == BoundaryCondition::PERIODIC) ? 0 : i + 1;
    
    return (u_current[ip1] - u_current[im1]) / (2.0 * dx);
}

double AdvectionDiffusionSolver2D::upwindDifference(
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

double AdvectionDiffusionSolver2D::quickDifference(
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

void AdvectionDiffusionSolver2D::applyBoundaryConditions(std::vector<double>& u_vec) {
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

double AdvectionDiffusionSolver2D::getCFL() const {
    return std::abs(c) * dx;  // CFL number (moltiplicare per dt per ottenere condizione)
}

double AdvectionDiffusionSolver2D::getDiffusionNumber() const {
    return D / (dx * dx);  // Diffusion number (moltiplicare per dt)
}

void AdvectionDiffusionSolver2D::saveCurrentTimeStep() {
    std::ofstream file(outputfile, std::ios::app);
    if (!file){
        std::cerr << "Error: impossible to open file " << outputfile << std::endl;
        return;
    }
    file << t << ",";
    for (int i = 0; i < nx; i++){
        file << x[i] << "," << u[i] << "," << i;
    }
    file << '\n';
    file.close();
    std::cout << "Current Time step saved into " << outputfile << std::endl;
}

void AdvectionDiffusionSolver2D::createOutputFile(){
    std::ofstream file(outputfile);
    if (!file){
        std::cerr << "Error: impossible to create file " << outputfile << std::endl;
        return;
    }
    file << "t,";
    for (int i = 0; i < nx; i++){
        file << "x" << i << ",u" << i;
    }
    file << '\n';
    file.close();
}

// void AdvectionDiffusionSolver1D::saveToFile(const std::string& filename) const {
//     std::ofstream file(filename);
//     if (!file.is_open()) {
//         std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
//         return;
//     }
    
//     file << "# x u(x,t=" << t << ")" << std::endl;
//     file << std::scientific << std::setprecision(10);
    
//     for (int i = 0; i < nx; ++i) {
//         file << x[i] << " " << u[i] << std::endl;
//     }
    
//     file.close();
//     std::cout << "Soluzione salvata in " << filename << std::endl;
// }

// void AdvectionDiffusionSolver1D::printStats() const {
//     double uMin = *std::min_element(u.begin(), u.end());
//     double uMax = *std::max_element(u.begin(), u.end());
//     double uMean = 0.0;
//     for (double val : u) {
//         uMean += val;
//     }
//     uMean /= nx;
    
//     std::cout << "\nStatistiche soluzione (t = " << t << "):" << std::endl;
//     std::cout << "  u_min  = " << uMin << std::endl;
//     std::cout << "  u_max  = " << uMax << std::endl;
//     std::cout << "  u_mean = " << uMean << std::endl;
// }
