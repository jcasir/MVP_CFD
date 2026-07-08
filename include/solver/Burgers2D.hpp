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
 * 2D solver for the Burgers' equations:
 * ∂u/∂t + (u·∇)u = D·∇²u
 */

class Burgers2D : public BaseSolver{
public:
	Burgers2D(const ConfigParser& config);
	~Burgers2D() override; 

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

    struct BoundaryConditionValues {
        double bottom, top, left, right;
    };
    BoundaryConditionValues u_bcs, v_bcs;

   	// Solution
	std::vector<double> u;     // x-component of velocity field. 
	std::vector<double> v;     // y-component of velocity field.

    // Intermediate state vectors for RK4
    std::vector<double> u_temp; 
    std::vector<double> v_temp;

    // Intermediate RK4 stages 
    std::vector<double> k1_x, k2_x, k3_x, k4_x;
	std::vector<double> k1_y, k2_y, k3_y, k4_y;

    //Output class handler
    PVDWriter pvdWriter;
    
    // Private methods for calculations

    // Indexing for 1D array
    int idx(int i,int j) const;
    // Handling of the neighbor values, ghost cells included
    double neighborX(const std::vector<double>& f, const BoundaryConditionValues& bcs,
                 int i, int j, int di) const;
    double neighborY(const std::vector<double>& f, const BoundaryConditionValues& bcs,
                 int i, int j, int dj) const;
	std::vector<double> computeRHS(const std::vector<double>& field,
                                const std::vector<double>& u_curr,
                                const std::vector<double>& v_curr,
                                const BoundaryConditionValues& bcs) const;
    void applyBoundaryConditions(std::vector<double>& u_vec, const BoundaryConditionValues& bcs);

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