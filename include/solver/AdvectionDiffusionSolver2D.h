#pragma once

#include "MathSettings.h"
#include "ConfigParser.h"
#include "InitialConditions.h"
#include "BaseSolver.h"
#include <stdexcept>
#include <vector>
#include <string>
#include <functional>

/**
 * 2D solver for the linear convection equation:
 * ∂u/∂t + c*∂u/∂x + c*∂u/∂y = 0
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
    
    // Output
    void saveCurrentTimeStep();
    void createOutputFile() override;
    // void saveToFile(const std::string& filename) const override;
    // void printStats() const override;
    
    // Getters
    double getCFL() const override;
    double getDiffusionNumber() const override;
    // const std::vector<double>& getSolution() const override { return u; }
    // const std::vector<double>& getGrid() const override { return x; }
    // double getCurrentTime() const override { return t; }
    
private:

    //Additional parameters for the domain setup
    int ny;
    int nx;                         // Grid points
    double Lx;                       
    double Ly;                      // Domain lenght
    double dx;                      
    double dy;                      // Grid spacing
    std::vector<double> x;   
    std::vector<double> y;          // Grid coordinates
    
    // Private methods for calculations
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    double advectionTerm(const std::vector<double>& u_current, int i) const;
    double diffusionTerm(const std::vector<double>& u_current, int i) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);
    
    // Spatial schemes
    double centralDifference(const std::vector<double>& u_current, int i) const;
    double upwindDifference(const std::vector<double>& u_current, int i) const;
    double quickDifference(const std::vector<double>& u_current, int i) const;

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition2D> initialCondition;
};