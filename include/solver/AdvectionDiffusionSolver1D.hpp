#pragma once

#include "MathSettings.hpp"
#include "ConfigParser.hpp"
#include "InitialConditions.hpp"
#include "BaseSolver.hpp"
#include "ErrorHandler.hpp"
#include "outputWriter.hpp"
#include <stdexcept>
#include <vector>
#include <memory>
#include <string>
#include <functional>

/**
 * 1D solver for the advection-diffusion equation:
 * ∂u/∂t + c*∂u/∂x = D*∂²u/∂x²
 */
class AdvectionDiffusionSolver1D : public BaseSolver {
public:
    
    // Constructor
    AdvectionDiffusionSolver1D(const ConfigParser& cfg);
    
    // Configuration methods
    void setInitialCondition() override;
    void setBoundaryConditions() override;
    
    // Solver methods
    void solve() override;
    void step(double dt) override;
    
    // Output
    void saveCurrentTimeStep();
    void createOutputFile();
    
    // Getters
    double getCFL() const override;
    double getDiffusionNumber() const override;
    
private:

    // Domain parameters
    double dx;                      // Grid spacing
    std::vector<double> x;          // Grid coordinates
    int nx;                         // Number of grid points
    double L;                       // Domain lenght

    // Solution
    std::vector<double> u;          // Current solution

    // Intermediate state vectors for RK4
    std::vector<double> u_temp;
    
    // Intermediate RK4 stages 
    std::vector<double> k1, k2, k3, k4;

    // Physical parameters
    double c;                       // Advection speed

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition1D> initialCondition;

    // Output handler
    OutputWriter1D outputWriter;
    
    // Private methods for calculations
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);
    
    // Spatial schemes
    template<typename FuncX>
    double advectionTerm(FuncX Ux, double ui) const
    {
        if (spatialScheme == SpatialScheme::CENTRAL) {
            return centralDifference(Ux);
        } else if (spatialScheme == SpatialScheme::UPWIND) {
            return upwindDifference(Ux, ui);
        } else { // QUICK
            return quickDifference(Ux, ui);
        }
    }

    template<typename FuncX>
    double diffusionTerm(FuncX Ux, double ui) const
    {
        // Central differencies of the second order for the diffusion.
        return D * (Ux(+1) - 2.0*ui + Ux(-1)) / (dx * dx);
    }

    template<typename FuncX>
    double centralDifference(FuncX Ux) const
    {
        return c * (Ux(+1) - Ux(-1)) / (2.0 * dx);
    }

    template<typename FuncX>
    double upwindDifference(FuncX Ux, double ui) const
    {
        return c * ((c >= 0) ? (ui - Ux(-1)) / dx : (Ux(+1) - ui) / dx);
    }

    template<typename FuncX>
    double quickDifference(FuncX Ux, double ui) const
    {
        // QUICK Scheme (Quadratic Upstream Interpolation for Convective Kinematics)

        double dux = (c >= 0)
            ? (Ux(-2) - 7*Ux(-1) + 3*ui + 3*Ux(+1))  / (8.0 * dx)
            : ( -3*Ux(-1) - 3*ui +7*Ux(+1) - Ux(+2)) / (8.0 * dx);

        return c * dux;
    }

};