#include "solver/IncNS2D.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

IncNS2D::IncNS2D(const ConfigParser& config)
 : BaseSolver(config), pvdWriter(output_dir + output_file + ".pvd")
{

    nx          = m_cfg.getInt("GRID_POINTS_X");
    ny          = m_cfg.getInt("GRID_POINTS_Y");
    Lx          = m_cfg.getDouble("DOMAIN_LENGHT_X");
    Ly          = m_cfg.getDouble("DOMAIN_LENGHT_Y");


    // Grid initialization
    // Non-overlapping grid for PERIODIC (domain [0,L)), overlapping for DIRICHLET/NEUMANN (nodes at exact boundaries)
    dx = (bcType == BoundaryCondition::PERIODIC) ? Lx / nx : Lx / (nx - 1);
    dy = (bcType == BoundaryCondition::PERIODIC) ? Ly / ny : Ly / (ny - 1);

    x.resize(nx);
    y.resize(ny);

    u.resize(nx * ny, 0.0);
    v.resize(nx * ny, 0.0);
    p.resize(nx * ny, 0.0);

    if (timeScheme == TimeScheme::RK4){
	    u_temp.resize(nx * ny, 0.0);
	    v_temp.resize(nx * ny, 0.0);

	    k1_x.resize(nx * ny, 0.0);
	    k2_x.resize(nx * ny, 0.0);
	    k3_x.resize(nx * ny, 0.0);
	    k4_x.resize(nx * ny, 0.0);

	    k1_y.resize(nx * ny, 0.0);
	    k2_y.resize(nx * ny, 0.0);
	    k3_y.resize(nx * ny, 0.0);
	    k4_y.resize(nx * ny, 0.0);
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

    std::cout << "Incompressible Navier Stokes Solver 2D Initialised" << std::endl;

}

IncNS2D::~IncNS2D() {
    try {
        pvdWriter.write();
    } catch (const std::exception& e) {
        std::cerr << "ERROR writing PVD: " << e.what() << "\n";
    }
}

int IncNS2D::idx(int i,int j) const{
    return (i * ny + j);
}

void IncNS2D::setInitialCondition() {

    // TODO: pressure should have a dedicated IC (e.g. solved from ∇²p=0 given the initial velocity field).
    // For Burgers this was acceptable since there is no pressure variable, but for incompressible NS
    // a consistent pressure initialization matters — though in practice most simulations start from rest
    // (u=v=p=0) so this rarely causes issues.
    // Any fix should live inside InitialCondition, not here.
    initialCondition->setIC(u,x,y);
    initialCondition->setIC(v,x,y);
    initialCondition->setIC(p,x,y);
    std::cout << "Initial condition set." << std::endl;

    // Checking whether the CFL or the diffusion number are too high (for explicit schemes)
    checkStability();

    //print initial condition
    // VTUWriter wants the number of cells so it must be given (nx - 1) because nx is the number of points
    // 1 for nz is the default to set the dimension to 2D.
    std::string outputFile = makeVTUFilename(output_file,0);
    VTUWriter vtuWriter(output_dir + outputFile,(nx - 1),(ny - 1),1,dx,dy,0.0);
    vtuWriter.addVector("u",u,v);
    vtuWriter.addScalar("Pressure",p);
    vtuWriter.write();
    pvdWriter.addStep(0.0,outputFile);
}

void IncNS2D::setBoundaryConditions()
{
    // Lambda function to set the Boundary conditions for each direction of the velocity.
    // Repetition is reduced since the process of printing and setting the boundary conditions is similar for
    // both direction of the velocity.
    auto makeBCs = [](const std::array<double, 4>& arrVal, 
                    const std::array<char, 4>& arrType,
                    const std::string& velName)
        -> BoundaryConditionValues
    {
        constexpr std::string_view sides[] = {"Bottom", "Top", "Left", "Right"};
        // For loop to print the boundary condition values and type for each side
        for (int i = 0; i < 4; ++i)
            std::cout << varName << " " << sides[i] << ": "
                      << arrType[i] << " = " << arrVal[i] << "\n";

        // this line here is the one that really sets the boundary condition inside the designated struct
        return { {arrVal[0],arrType[0]}, {arrVal[1],arrType[1]},
                 {arrVal[2],arrType[2]}, {arrVal[3],arrType[3]} };
    };

    std::cout << "Boundary conditions set: ";

    // For periodic boundary conditions, all individual boundary settings are skipped as they are not needed.
    if (bcType == BoundaryCondition::PERIODIC) {
        std::cout << "Periodic\n";
        return;
    }

    // First argument: array of boundary condition values on each side.
    // Second argument: array of boundary condition types on each side.
    // Third argument: name of the variable.
    // The return value from getBCs() is passed directly as an argument to improve compactness
    // and reduce memory allocation, even though this is not necessary since the temporary arrays
    // inside this function would be destroyed at scope exit.
    u_bcs = makeBCs(m_cfg.getBCs<double,4>("BOUNDARY_CONDITIONS_VALUES_U"),
                    m_cfg.getBCs<char,4>("BOUNDARY_CONDITIONS_TYPES_U"),
                    "U");
    v_bcs = makeBCs(m_cfg.getBCs<double,4>("BOUNDARY_CONDITIONS_VALUES_V"),
                    m_cfg.getBCs<char,4>("BOUNDARY_CONDITIONS_TYPES_V"),
                    "V");
    p_bcs = makeBCs(m_cfg.getBCs<double,4>("BOUNDARY_CONDITIONS_VALUES_P"),
                    m_cfg.getBCs<char,4>("BOUNDARY_CONDITIONS_TYPES_P"),
                    "P");
    std::cout << std::endl;
}

void IncNS2D::solve() {
    std::cout << "\nStarting solver:" << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    std::cout << "  t_end = " << tEnd << std::endl;
    std::cout << "  Diffusion number = " << getDiffusionNumber() << std::endl;
    
    int nSteps = 0;
    while (t < tEnd) {
        const double cfl = getCFL();
        if (cfl > 1.0) {
            std::cerr << "CFL exceeded limit at step " << nSteps << "\n";
            throw StabilityException("CFL", cfl, 1.0);
        }
        double dt_actual = std::min(dt, tEnd - t);
        step(dt_actual);
        nSteps++;
        
        if (nSteps % 100 == 0) {
            std::cout << "  Step " << nSteps << ", t = " << t << std::endl;
        }
        if (nSteps % output_freq == 0){
            std::string outputFile = makeVTUFilename(output_file,nSteps);
            VTUWriter outputWriter(output_dir + outputFile,(nx - 1),(ny - 1),1,dx,dy,0.0);
            outputWriter.addVector("u",u,v);
            outputWriter.write();
            pvdWriter.addStep(t,outputFile);
        }
    }
    
    std::cout << "Resolution completed after " << nSteps << " time steps." << std::endl;
}

void IncNS2D::step(double dt_actual) {
    if (timeScheme == TimeScheme::EULER_EXPLICIT) {
        // Esplicit Euler: u^(n+1) = u^n + dt * RHS(u^n)
        std::vector<double> rhs_x = computeRHS(u,u,v,u_bcs);
        std::vector<double> rhs_y = computeRHS(v,u,v,v_bcs);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u[idx(i,j)] += dt_actual * rhs_x[idx(i,j)];
                v[idx(i,j)] += dt_actual * rhs_y[idx(i,j)];
            }
        }
        
        applyBoundaryConditions(u,u_bcs);
        applyBoundaryConditions(v,v_bcs);
        
    } else if (timeScheme == TimeScheme::RK4) {
        // Runge-Kutta 4° order
        k1_x = computeRHS(u,u,v,u_bcs);
        k1_y = computeRHS(v,u,v,v_bcs);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u_temp[idx(i,j)] = u[idx(i,j)] + 0.5 * dt_actual * k1_x[idx(i,j)];
                v_temp[idx(i,j)] = v[idx(i,j)] + 0.5 * dt_actual * k1_y[idx(i,j)];
            }
        }
        applyBoundaryConditions(u_temp,u_bcs);
        applyBoundaryConditions(v_temp,v_bcs);

        k2_x = computeRHS(u_temp,u_temp,v_temp,u_bcs);
        k2_y = computeRHS(v_temp,u_temp,v_temp,v_bcs);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u_temp[idx(i,j)] = u[idx(i,j)] + 0.5 * dt_actual * k2_x[idx(i,j)];
                v_temp[idx(i,j)] = v[idx(i,j)] + 0.5 * dt_actual * k2_y[idx(i,j)];
            }
        }
        applyBoundaryConditions(u_temp,u_bcs);
        applyBoundaryConditions(v_temp,v_bcs);

        k3_x = computeRHS(u_temp,u_temp,v_temp,u_bcs);
        k3_y = computeRHS(v_temp,u_temp,v_temp,v_bcs);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u_temp[idx(i,j)] = u[idx(i,j)] + dt_actual * k3_x[idx(i,j)];
                v_temp[idx(i,j)] = v[idx(i,j)] + dt_actual * k3_y[idx(i,j)];
            }
        }
        applyBoundaryConditions(u_temp,u_bcs);
        applyBoundaryConditions(v_temp,v_bcs);

        k4_x = computeRHS(u_temp,u_temp,v_temp,u_bcs);
        k4_y = computeRHS(v_temp,u_temp,v_temp,v_bcs);
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j){
                u[idx(i,j)] += (dt_actual / 6.0) * (k1_x[idx(i,j)] + 2.0*k2_x[idx(i,j)] + 2.0*k3_x[idx(i,j)] + k4_x[idx(i,j)]);
                v[idx(i,j)] += (dt_actual / 6.0) * (k1_y[idx(i,j)] + 2.0*k2_y[idx(i,j)] + 2.0*k3_y[idx(i,j)] + k4_y[idx(i,j)]);
            }
        }
        
        applyBoundaryConditions(u,u_bcs);
        applyBoundaryConditions(v,v_bcs);
    }
    
    t += dt_actual;
}

std::vector<double> IncNS2D::computeRHS(
    const std::vector<double>& field,
    const std::vector<double>& u_curr,
    const std::vector<double>& v_curr,
    const BoundaryConditionValues& bcs
    ) const
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
                if (ni >= 0 && ni < nx) return field[idx(ni, j)];
                if (bcType == BoundaryCondition::PERIODIC)
                    return field[idx((ni + nx) % nx, j)];  
                if (bcType == BoundaryCondition::DIRICHLET)
                    return (ni < 0) ? 2*bcs.left  - field[idx(-ni, j)]
                                    : 2*bcs.right - field[idx(2*(nx-1) - ni, j)];
                // NEUMANN
                return (ni < 0) ? field[idx(-ni,           j)] + ni * 2 * bcs.left  * dx
                                : field[idx(2*(nx-1) - ni, j)] + (ni - (nx-1)) * 2 * bcs.right * dx;
            };
            auto Uy = [&](int dj) -> double {
                int nj = j + dj;
                if (nj >= 0 && nj < ny) return field[idx(i, nj)];
                if (bcType == BoundaryCondition::PERIODIC)
                    return field[idx(i, (nj + ny) % ny)];
                if (bcType == BoundaryCondition::DIRICHLET)
                    return (nj < 0) ? 2*bcs.bottom - field[idx(i, -nj)]
                                    : 2*bcs.top    - field[idx(i, 2*(ny-1) - nj)];
                // NEUMANN
                return (nj < 0) ? field[idx(i, -nj)]           + nj * 2 * bcs.bottom * dy
                                : field[idx(i, 2*(ny-1) - nj)] + (nj - (ny-1)) * 2 * bcs.top * dy;
            };
            const double uij = field[idx(i,j)];
            const double vel_x = u_curr[idx(i,j)];
            const double vel_y = v_curr[idx(i,j)];
            // RHS = -(c_x*∂u/∂x + c_y*∂u/∂y) + D*(∂²u/∂x² + ∂²u/∂y²)
            rhs[idx(i,j)] = -advectionTerm(Ux, Uy, uij, vel_x, vel_y) + diffusionTerm(Ux, Uy, uij);
        }
    }
    
    return rhs;
}

void IncNS2D::applyBoundaryConditions(std::vector<double>& u_vec,const BoundaryConditionValues& bcs) {
    // Setting bottom boundary conditions
    if (bcs.bottom.type == 'D'){
        for (int i = 0; i < nx; ++i) {
            u_vec[idx(i, 0)] = bcs.bottom.value;  
        }
    }
    else{
        for (int i = 0; i < nx; ++i) {
            u_vec[idx(i, 0)] = u_vec[idx(i, 1)] - bcs.bottom.value * dy; 
        }   
    }

    // Setting top boundary conditions
    if (bcs.top.type == 'D'){
        for (int i = 0; i < nx; ++i) {
            u_vec[idx(i, ny - 1)] = bcs.top.value; 
        }
    }
    else{
        for (int i = 0; i < nx; ++i) {
            u_vec[idx(i,ny - 1)] = u_vec[idx(i,ny - 2)] + bcs.top.value * dy; 
        }   
    }

    // Setting left boundary conditions
    if (bcs.left.type == 'D'){
        for (int j = 0; j < ny; ++j) {
            u_vec[idx(0, j)] = bcs.left.value;  
        }
    }
    else{
        for (int j = 0; j < ny; ++j) {
            u_vec[idx(0, j)] = u_vec[idx(1, j)] - bcs.left.value * dx;
        }   
    }

    // Setting right boundary conditions
    if (bcs.right.type == 'D'){
        for (int j = 0; j < ny; ++j) {
            u_vec[idx(nx - 1,j)] = bcs.right.value; 
        }
    }
    else{
        for (int j = 0; j < ny; ++j) {
            u_vec[idx(nx - 1,j)] = u_vec[idx(nx - 2,j)] + bcs.right.value * dx; 
        }   
    }
}

double IncNS2D::getCFL() const {
	double maxU = *std::max_element(u.begin(), u.end(), 
	    [](double a, double b){ return std::abs(a) < std::abs(b); });
	double maxV = *std::max_element(v.begin(), v.end(),
	    [](double a, double b){ return std::abs(a) < std::abs(b); });
	return dt * (std::abs(maxU)/dx + std::abs(maxV)/dy);  // CFL number
}

double IncNS2D::getDiffusionNumber() const {
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