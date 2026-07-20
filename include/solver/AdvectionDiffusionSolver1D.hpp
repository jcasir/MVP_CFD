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
    
    // Retrieve the stability numbers for the checkStability() method in BaseSolver.
    StabilityNumbers getStabilityNumbers() const override;
    
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

    // Handling of the neighbor values, ghost cells included
    double neighbor(const std::vector<double>& f, int i, int di) const;
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);

    // Spatial schemes
    double advectionTerm(const std::vector<double>& field, int i) const;
    double diffusionTerm(const std::vector<double>& field, int i) const;
    double centralDifference(const std::vector<double>& field, int i) const;
    double upwindDifference(const std::vector<double>& field, int i) const;
    double quickDifference(const std::vector<double>& field, int i) const;

};