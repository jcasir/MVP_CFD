#pragma once

#include "MathSettings.hpp"
#include "ConfigParser.hpp"
#include "InitialConditions.hpp"
#include "BaseSolver.hpp"
#include "outputWriter.hpp"
#include <stdexcept>
#include <vector>
#include <memory>
#include <array>
#include <string>
#include <functional>

/**
 * @brief 2D solver for the Incompressible Navier-Stokes equations.
 * * Continuity equation (Incompressibility constraint):
 * ∂u/∂x + ∂v/∂y = 0
 * * Momentum equations (x and y components):
 * ∂u/∂t + (u·∇)u = -1/ρ·∇p + ν·∇²u
 * * Pressure-Poisson Equation (PPE) for mass conservation:
 * ∇²p = -ρ∇·[(u·∇)u]
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
    
    // Getters
    double getCFL() const override;
    double getDiffusionNumber() const override;

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

   	// Solution
	std::vector<double> u;     // x-component of velocity field. 
	std::vector<double> v;     // y-component of velocity field.
    std::vector<double> p;     // pressure component


    // Intermediate state vectors for RK4
    std::vector<double> u_temp; 
    std::vector<double> v_temp;

    // Intermediate RK4 stages 
    std::vector<double> k1_x, k2_x, k3_x, k4_x;
	std::vector<double> k1_y, k2_y, k3_y, k4_y;

    //Output class handler
    PVDWriter pvdWriter;
    
    // Private methods for calculations
    std::vector<double> computeRHS(const std::vector<double>& field,
                                const std::vector<double>& u_curr,
                                const std::vector<double>& v_curr,
                                const BoundaryConditionValues& bcs) const;
    void applyBoundaryConditions(std::vector<double>& u_vec, const BoundaryConditionValues& bcs);
    
    // Spatial schemes
    template<typename FuncX, typename FuncY>
    double advectionTerm(FuncX Ux, FuncY Uy, double uij, double vel_x, double vel_y) const
    {
        if (spatialScheme == SpatialScheme::CENTRAL) {
            return centralDifference(Ux, Uy, vel_x, vel_y);
        } else if (spatialScheme == SpatialScheme::UPWIND) {
            return upwindDifference(Ux, Uy, uij, vel_x, vel_y);
        } else { // QUICK
            return quickDifference(Ux, Uy, uij, vel_x, vel_y);
        }
    }

    template<typename FuncX, typename FuncY>
    double diffusionTerm(FuncX Ux, FuncY Uy, double uij) const
    {
        // Central differencies of the second order for the diffusion.
        return D * ((Ux(+1) - 2.0*uij + Ux(-1)) / (dx * dx)
                  + (Uy(+1) - 2.0*uij + Uy(-1)) / (dy * dy));

    }

    template<typename FuncX, typename FuncY>
    double centralDifference(FuncX Ux, FuncY Uy, double vel_x, double vel_y) const
    {
        return vel_x * (Ux(+1) - Ux(-1)) / (2.0 * dx) 
             + vel_y * (Uy(+1) - Uy(-1)) / (2.0 * dy);
    }

    template<typename FuncX, typename FuncY>
    double upwindDifference(FuncX Ux, FuncY Uy, double uij, double vel_x, double vel_y) const
    {
        double dux = (vel_x >= 0) ? (uij - Ux(-1)) / dx : (Ux(+1) - uij) / dx;
        double duy = (vel_y >= 0) ? (uij - Uy(-1)) / dy : (Uy(+1) - uij) / dy;
        return vel_x * dux + vel_y * duy;
    }

    template<typename FuncX, typename FuncY>
    double quickDifference(FuncX Ux, FuncY Uy, double uij, double vel_x, double vel_y) const
    {
        // QUICK Scheme (Quadratic Upstream Interpolation for Convective Kinematics)

        double dux = (vel_x >= 0)
            ? (Ux(-2) - 7*Ux(-1) + 3*uij + 3*Ux(+1))   / (8.0 * dx)
            : ( -3*Ux(-1) - 3*uij +7*Ux(+1) - Ux(+2)) / (8.0 * dx);

        double duy = (vel_y >= 0)
            ? (Uy(-2) - 7*Uy(-1) + 3*uij + 3*Uy(+1))   / (8.0 * dy)
            : ( -3*Uy(-1) - 3*uij +7*Uy(+1) - Uy(+2))  / (8.0 * dy);

        return vel_x * dux + vel_y * duy;
    }

    //Indexing for 1D array
    int idx(int i,int j) const;

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition2D> initialCondition;
};