#pragma once

#include "MathSettings.hpp"
#include "ConfigParser.hpp"
#include "InitialConditions.hpp"
#include "BaseSolver.hpp"
#include "outputWriter.hpp"
#include <stdexcept>
#include <vector>
#include <memory>
#include <string>
#include <functional>

/**
 * 2D solver for the advection diffusion equation:
 * // ∂u/∂t + c_x*∂u/∂x + c_y*∂u/∂y = D*(∂²u/∂x² + ∂²u/∂y²)
 */
class AdvectionDiffusionSolver2D : public BaseSolver {
public:
    AdvectionDiffusionSolver2D(const ConfigParser& config);
    ~AdvectionDiffusionSolver2D() override; 
    
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
    double bcTop;                   // Boundary conditions exclusive to 2D problem
    double bcBottom;                     

    // Physical parameters
    double c_x;                       // Advection speed in x
    double c_y;                       // Advection speed in y

    // Solution
    std::vector<double> u;          // Current solution

    // Intermediate state vectors for RK4
    std::vector<double> u_temp;

    // Intermediate RK4 stages 
    std::vector<double> k1, k2, k3, k4;

    //Output class handler
    PVDWriter pvdWriter;
    
    // Private methods for calculations

    // Indexing for 1D array
    int idx(int i,int j) const;
    // Handling of the neighbor values, ghost cells included
    double neighborX(const std::vector<double>& f, int i, int j, int di) const;
    double neighborY(const std::vector<double>& f, int i, int j, int dj) const;
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);

    // Spatial schemes
    double advectionTerm(const std::vector<double>& field, int i, int j) const;
    double diffusionTerm(const std::vector<double>& field, int i, int j) const;
    double centralDifference(const std::vector<double>& field, int i, int j) const;
    double upwindDifference(const std::vector<double>& field, int i, int j) const;
    double quickDifference(const std::vector<double>& field, int i, int j) const;

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition2D> initialCondition;
};