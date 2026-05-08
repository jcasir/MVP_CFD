#pragma once

#include "MathSettings.hpp"
#include "ConfigParser.hpp"
#include "InitialConditions.hpp"
#include "BaseSolver.hpp"
#include "vtuWriter.hpp"
#include <stdexcept>
#include <vector>
#include <string>
#include <functional>

/**
 * 2D solver for the advection diffusion equation:
 * // ∂u/∂t + c_x*∂u/∂x + c_y*∂u/∂y = D*(∂²u/∂x² + ∂²u/∂y²)
 */
class AdvectionDiffusionSolver2D : public BaseSolver {
public:
    AdvectionDiffusionSolver2D(const ConfigParser& config);
    
    // Configuration methods
    void setInitialCondition() override;
    void setBoundaryConditions() override;
    
    // Resolution methods
    void solve() override;
    void step(double dt) override;

    //Output
    void finalOutput() override;
    
    // Getters
    double getCFL() const override;
    double getDiffusionNumber() const override;
    // const std::vector<double>& getSolution() const override { return u; }
    // const std::vector<double>& getGrid() const override { return x; }
    // double getCurrentTime() const override { return t; }
    
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
    double bcTop;                   // Boundary conditions exclusive to 2D problem
    double bcBottom;                     

    // Physical parameters
    double c_x;                       // Advection speed in x
    double c_y;                       // Advection speed in y

    //Output class handler
    PVDWriter pvdWriter;
    
    // Private methods for calculations
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);
    
    // Spatial schemes
    template<typename FuncX, typename FuncY>
    double advectionTerm(FuncX Ux, FuncY Uy, double uij) const
    {
        if (spatialScheme == SpatialScheme::CENTRAL) {
            return centralDifference(Ux, Uy);
        } else if (spatialScheme == SpatialScheme::UPWIND) {
            return upwindDifference(Ux, Uy, uij);
        } else { // QUICK
            return quickDifference(Ux, Uy, uij);
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
    double centralDifference(FuncX Ux, FuncY Uy) const
    {
        return c_x * (Ux(+1) - Ux(-1)) / (2.0 * dx) 
             + c_y * (Uy(+1) - Uy(-1)) / (2.0 * dy);
    }

    template<typename FuncX, typename FuncY>
    double upwindDifference(FuncX Ux, FuncY Uy, double uij) const
    {
        double dux = (c_x >= 0) ? (uij - Ux(-1)) / dx : (Ux(+1) - uij) / dx;
        double duy = (c_y >= 0) ? (uij - Uy(-1)) / dy : (Uy(+1) - uij) / dy;
        return c_x * dux + c_y * duy;
    }

    template<typename FuncX, typename FuncY>
    double quickDifference(FuncX Ux, FuncY Uy, double uij) const
    {
        // QUICK Scheme (Quadratic Upstream Interpolation for Convective Kinematics)

        double dux = (c_x >= 0)
            ? (-Ux(-2) + 8*Ux(-1) - 8*Ux(+1) + uij)   / (12.0 * dx)
            : ( -uij   + 8*Ux(+1) - 8*Ux(-1) + Ux(+2)) / (12.0 * dx);

        double duy = (c_y >= 0)
            ? (-Uy(-2) + 8*Uy(-1) - 8*Uy(+1) + uij)   / (12.0 * dy)
            : ( -uij   + 8*Uy(+1) - 8*Uy(-1) + Uy(+2)) / (12.0 * dy);

        return c_x * dux + c_y * duy;
    }

    //Indexing for 1D array
    int idx(int i,int j) const;

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition2D> initialCondition;
};