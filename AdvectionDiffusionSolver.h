#pragma once

#include "math_settings.h"
#include "ConfigParser.h"
#include "InitialConditions.h"
#include <vector>
#include <string>
#include <functional>

/**
 * Risolutore 1D per l'equazione di avvezione-diffusione:
 * ∂u/∂t + c*∂u/∂x = D*∂²u/∂x²
 */
class AdvectionDiffusionSolver {
public:
    
    // Costruttore
    AdvectionDiffusionSolver(
        int nx,                     // Numero di punti griglia
        double L,                   // Lunghezza dominio
        double c,                   // Velocità avvezione
        double D,                   // Coefficiente diffusione
        SpatialScheme spatialScheme = SpatialScheme::CENTRAL,
        TimeScheme timeScheme = TimeScheme::EULER_EXPLICIT
    );
    
    // Metodi di configurazione
    void setInitialCondition(std::function<double(double)> ic);
    void setBoundaryConditions(
        BoundaryCondition bcType,
        double leftValue = 0.0,
        double rightValue = 0.0
    );
    
    // Metodi di risoluzione
    void solve(double dt, double tEnd);
    void step(double dt);
    
    // Output
    void saveToFile(const std::string& filename) const;
    void printStats() const;
    
    // Getters
    const std::vector<double>& getSolution() const { return u; }
    const std::vector<double>& getGrid() const { return x; }
    double getCurrentTime() const { return t; }
    double getCFL() const;
    double getDiffusionNumber() const;
    
private:
    // Parametri del dominio
    int nx;                         // Numero di punti
    double L;                       // Lunghezza dominio
    double dx;                      // Spaziatura griglia
    std::vector<double> x;          // Coordinate griglia
    
    // Parametri fisici
    double c;                       // Velocità avvezione
    double D;                       // Coefficiente diffusione
    
    // Soluzione
    std::vector<double> u;          // Soluzione corrente
    double t;                       // Tempo corrente
    
    // Schemi numerici
    SpatialScheme spatialScheme;
    TimeScheme timeScheme;
    
    // Condizioni al contorno
    BoundaryCondition bcType;
    double bcLeft, bcRight;
    
    // Metodi privati per i calcoli
    std::vector<double> computeRHS(const std::vector<double>& u_current) const;
    double advectionTerm(const std::vector<double>& u_current, int i) const;
    double diffusionTerm(const std::vector<double>& u_current, int i) const;
    void applyBoundaryConditions(std::vector<double>& u_vec);
    
    // Schemi spaziali
    double centralDifference(const std::vector<double>& u_current, int i) const;
    double upwindDifference(const std::vector<double>& u_current, int i) const;
    double quickDifference(const std::vector<double>& u_current, int i) const;
};

