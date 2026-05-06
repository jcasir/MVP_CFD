#pragma once

#include "MathSettings.h"
#include "ConfigParser.h"
#include "InitialConditions.h"
#include "BaseSolver.h"
#include "ErrorHandler.h"
#include <stdexcept>
#include <vector>
#include <string>
#include <functional>

/**
 * Risolutore 1D per l'equazione di avvezione-diffusione:
 * ∂u/∂t + c*∂u/∂x = D*∂²u/∂x²
 */
class AdvectionDiffusionSolver1D : public BaseSolver {
public:
    
    // Costruttore
    AdvectionDiffusionSolver1D(const ConfigParser& cfg);
    
    // Metodi di configurazione
    void setInitialCondition() override;
    void setBoundaryConditions() override;
    
    // Metodi di risoluzione
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
    
    // Private methods for calculations
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    double advectionTerm(const std::vector<double>& u_current, int i) const;
    double diffusionTerm(const std::vector<double>& u_current, int i) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);
    
    // Schemi spaziali
    double centralDifference(const std::vector<double>& u_current, int i) const;
    double upwindDifference(const std::vector<double>& u_current, int i) const;
    double quickDifference(const std::vector<double>& u_current, int i) const;

    // smart pointer for initial condition class
    std::unique_ptr<IInitialCondition1D> initialCondition;
};

