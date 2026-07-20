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

    // Largest ratio condition/limit among the von Neumann conditions of the current
    // scheme: < 1 means stable, >= 1 unstable. Cheap enough to print at every step.
    double getStabilityMargin() const;

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

    // Methods to check stability. The logic is implemented in the base solver,
    // each specific solver implements his own stability numbers which will then be used by
    // checkStability() with the logic for the condition only in this function.
    // With report = false it stays silent unless a condition is violated (throws):
    // meant for the per-step calls inside the solve() loops, where printing the
    // full block at every step would flood the log.
    void checkStability(bool report = true) const;
    bool getExplicitSchemeFlag() const;

    // Per-direction stability numbers (1D solvers leave the y entries at zero).
    // Filled by each solver, consumed by checkStability().
    struct StabilityNumbers {
        double Cx = 0.0, Cy = 0.0;   // per-direction Courant numbers  c_i*dt/dx_i
        double dx = 0.0, dy = 0.0;   // per-direction diffusion numbers D*dt/dx_i^2
    };
    virtual StabilityNumbers getStabilityNumbers() const = 0;
};