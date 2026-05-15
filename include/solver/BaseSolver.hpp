#pragma once

#include "MathSettings.hpp"
#include "ConfigParser.hpp"
#include "ErrorHandler.hpp"
#include <stdexcept>
#include <vector>
#include <string>
#include <functional>



class BaseSolver {
public:
    
    // Constructor
    BaseSolver(const ConfigParser& config);

    virtual ~BaseSolver() = default;
    
    // Configuration methods
    virtual void setInitialCondition() = 0;
    virtual void setBoundaryConditions() = 0;
    
    // Solver methods
    virtual void solve() = 0;
    virtual void step(double dt) = 0;

    virtual double getCFL() const = 0;
    virtual double getDiffusionNumber() const = 0;

protected:
    //config
    const ConfigParser& m_cfg;

    // Domain parameters
    std::string config_file;        // Config file
    double dt;                      // Time step
    double tEnd;                    // End time
    int time_iter;
    double t;                       // Current time
    
    // Physical parameters
    double D;                       // Diffusion coefficient
    
    // Numerical schemes
    SpatialScheme spatialScheme;
    TimeScheme timeScheme;
    
    // Boundary Conditions
    BoundaryCondition bcType;
    double bcLeft, bcRight;

    //output flag for debug
    bool verbose;

    //output file
    std::string output_file;
    std::string output_dir;
    std::string mesh_file;
    int output_freq;

    //Methods to check stability
    void checkStability() const;
    bool getExplicitSchemeFlag() const;
};