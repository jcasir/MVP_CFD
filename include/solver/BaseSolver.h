#pragma once

#include "math_settings.h"
#include "ConfigParser.h"
#include <vector>
#include <string>
#include <functional>



class BaseSolver {
public:
    
    // Costruttore
    BaseSolver(const ConfigParser& cconfigfg);

    virtual ~BaseSolver() = default;
    
    // Metodi di configurazione
    virtual void setInitialCondition() = 0;
    virtual void setBoundaryConditions() = 0;
    
    // Metodi di risoluzione
    virtual void solve() = 0;
    virtual void step(double dt) = 0;
    
    // Output
    virtual void createOutputFile() = 0;
    // virtual void saveToFile(const std::string& filename) const = 0;
    // virtual void printStats() const = 0;
    
    // Getters
    virtual const std::vector<double>& getSolution() const { return u; }
    virtual const std::vector<double>& getGrid() const { return x; }
    virtual double getCurrentTime() const { return t; }
    virtual double getCFL() const = 0;
    virtual double getDiffusionNumber() const = 0;

protected:
    //config
    const ConfigParser& m_cfg;

    // Parametri del dominio
    int nx;                         // Numero di punti griglia
    double L;                       // Lunghezza dominio
    double dx;                      // Spaziatura griglia
    std::vector<double> x;          // Coordinate griglia
    std::string config_file;        // Config file
    double dt;                      // Time step
    double tEnd;                    // End time
    int time_iter;
    
    // Parametri fisici
    double c;                       // Velocità avvezione
    double D;                       // Coefficiente diffusione
    
    // Soluzione
    std::vector<double> u;          // Soluzione corrente
    double t;                       // Tempo corrente
    
    // Numerical schemes
    SpatialScheme spatialScheme;
    TimeScheme timeScheme;
    
    // Boundary Conditions
    BoundaryCondition bcType;
    double bcLeft, bcRight;

    //Initial conditions
    InitialCondition initialCondition;

    //output flag for debug
    bool verbose;

    //output file
    std::string outputfile;
    std::string mesh_file;
    int output_freq;

    //Methods to check stability
    void checkStability() const;
    bool getExplicitSchemeFlag() const;
};

