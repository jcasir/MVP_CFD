#pragma once

#include "MathSettings.hpp"
#include "ConfigParser.hpp"
#include "InitialConditions.hpp"
#include "BaseSolver.hpp"
#include "outputWriter.hpp"
#include "ErrorHandler.hpp"
#include <stdexcept>
#include <vector>
#include <memory>
#include <array>
#include <string>
#include <functional>

/**
 * @brief 2D solver for the Incompressible Navier-Stokes equations.
 * Continuity equation (Incompressibility constraint):
 * div(u) = 0                                                                       [1]
 * Momentum equations:
 * ddt(u) + div(u,u) = -1/rho * grad(p) + D * lap(u)                                [2]
 * 
 * Discretized in time (note: the pressure term carries Δt like every other RHS term, and the
 * pressure is taken at n+1 because there is no equation to advance it explicitly — it will be
 * determined by the incompressibility constraint):
 *
 * u[n+1] = u[n] + Δt * (−div(u[n],u[n]) + D * lap(u[n])) - Δt/rho * grad(p[n+1])   [3]
 *
 * From this equation we can compute an intermediate velocity. This intermediate velocity is computed
 * ignoring the pressure component. So we are basically advancing in time with everything, except
 * the pressure. We call this velocity u_star
 *
 * u_star = u[n] + Δt * (−div(u[n],u[n]) + D * lap(u[n]))                           [4]
 *
 * So the correct velocity on the next step can be written as:
 *
 * u[n+1] = u_star - Δt/rho * grad(p[n+1])                                          [5]
 *
 * This equation has two unknowns: u[n+1] and p[n+1]. We first compute the pressure by imposing the
 * incompressibility constraint on this equation. In particular we take the divergence of this equation
 * and we impose div(u[n+1]) = 0. This leads to the Pressure-Poisson Equation (PPE):
 *
 *  lap(p[n+1]) = (rho/Δt) * div(u_star)                                            [6]
 * 
 * Since the PPE was derived by imposing that the divergence of the velocity at the next step is zero,
 * solving it and applying the corrector enforces the incompressibility constraint: the corrected field 
 * is divergence-free up to the tolerance reached by the iterations. Discretizing the equation yields a 
 * linear system, which we solve iteratively (Gauss-Seidel).
 *
 * After p[n+1] is determined we go back to the correct velocity equation and finally compute u[n+1]:
 *
 * u[n+1] = u_star - Δt/rho * grad(p[n+1])                                          [7]
 *
 * The Helmholtz–Hodge theorem says that any vector field (with suitable BCs) splits uniquely
 * into a divergence-free part plus a gradient part:
 *
 * u_star = u_div-free + grad(phi)                                                  [8]
 *
 * So the predictor [4] produces a generic u_star. The PPE computes precisely that scalar phi
 * (here phi = (Δt/rho) * p) whose gradient is the "compressible" part of u_star. The corrector
 * subtracts it, leaving the divergence-free part. Subtracting the gradient component is an
 * orthogonal projection onto the space of divergence-free fields, hence the name projection method.
 * 
 */

class IncNS2D : public BaseSolver{
public:
	IncNS2D(const ConfigParser& config);
	~IncNS2D() override; 

    // Configuration methods
    void setInitialCondition() override;
    void setBoundaryConditions() override;
    
    // Resolution methods
    void solve() override;
    void step(double dt) override;
    
    // Retrieve the stability numbers for the checkStability() method in BaseSolver.
    StabilityNumbers getStabilityNumbers() const override;

private:
    //Additional parameters for the domain setup
    int ny;                         // Grid points
    int nx;                         
    double Lx;                      // Domain lenght
    double Ly;                      
    double dx;                      // Grid spacing
    double dy;                     
    std::vector<double> x;          // Grid coordinates
    std::vector<double> y;          

    struct SideBC { double value; char type; };
    struct BoundaryConditionValues { SideBC bottom, top, left, right; };
    BoundaryConditionValues u_bcs, v_bcs, p_bcs;  

    // Physical parameters   
    double rho;                     // Density is constant since the solver is for the incompressible Navier Stokes equations

    // PPE settings and tolerance
    int ppe_max_iter;
    double ppe_toll;
    std::string ppe_toll_type;      // "D": ppe_toll is a divergence tolerance, "E": effective residual tolerance
    double omega_sor;               // SOR relaxation factor; 1.0 = plain Gauss-Seidel

   	// Solution
	std::vector<double> u;     // x-component of velocity field. 
	std::vector<double> v;     // y-component of velocity field.
    std::vector<double> p;     // pressure component

    // Predictor velocity.
    std::vector<double> u_star; 
    std::vector<double> v_star;

    //PPE Source term.
    std::vector<double> b;

    // Intermediate state vectors for RK4
    std::vector<double> u_temp; 
    std::vector<double> v_temp;

    // Intermediate RK4 stages 
    std::vector<double> k1_x, k2_x, k3_x, k4_x;
	std::vector<double> k1_y, k2_y, k3_y, k4_y;

    //Output class handler
    PVDWriter pvdWriter;

    // Flag to handle all-Neumann pressure BCs
    bool Neumann_pressure_flag;
    
    // Private methods for calculations

    // Indexing for 1D array
    int idx(int i,int j) const;
    // Handling of the neighbor values, ghost cells included
    double neighborX(const std::vector<double>& f, const BoundaryConditionValues& bcs,
                 int i, int j, int di) const;
    double neighborY(const std::vector<double>& f, const BoundaryConditionValues& bcs,
                 int i, int j, int dj) const;
    double computeDiv(const std::vector<double>& u, const std::vector<double>& v, double dx, double dy, int i, int j) const;
    std::vector<double> computeRHS(const std::vector<double>& field,
                                const std::vector<double>& u_curr,
                                const std::vector<double>& v_curr,
                                const BoundaryConditionValues& bcs) const;
    void solvePressurePoisson(double dt_actual);
    double ppeResidual() const;
    void applyBoundaryConditions(std::vector<double>& field, const BoundaryConditionValues& bcs);
    
    // Spatial schemes
    double advectionTerm(const std::vector<double>& field, const BoundaryConditionValues& bcs,
                 double cx, double cy, int i, int j) const;
    double diffusionTerm(const std::vector<double>& field, const BoundaryConditionValues& bcs,
                 int i, int j) const;
    double centralDifference(const std::vector<double>& field, const BoundaryConditionValues& bcs,
                 double cx, double cy, int i, int j) const;
    double upwindDifference(const std::vector<double>& field, const BoundaryConditionValues& bcs,
                 double cx, double cy, int i, int j) const;
    double quickDifference(const std::vector<double>& field, const BoundaryConditionValues& bcs,
                 double cx, double cy, int i, int j) const;

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition2D> initialCondition;
};