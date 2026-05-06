#include "solver/BaseSolver.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

BaseSolver::BaseSolver(const ConfigParser& config) : m_cfg(config)
    // int nx, double L, double c, double D,
    // SpatialScheme spatialScheme, TimeScheme timeScheme)
    // : nx(nx), L(L), c(c), D(D), 
    //   spatialScheme(spatialScheme), timeScheme(timeScheme),
    //   t(0.0), bcType(BoundaryCondition::DIRICHLET), bcLeft(0.0), bcRight(0.0)
{
    // Lettura config

    nx                      = m_cfg.getInt("GRID_POINTS");
    output_freq             = m_cfg.getInt("OUTPUT_FREQUENCY");
    time_iter               = 0;

    L                       = m_cfg.getDouble("DOMAIN_LENGHT");
    c                       = m_cfg.getDouble("ADVECTION_SPEED");
    D                       = m_cfg.getDouble("DIFFUSION_COEFF");
    bcLeft                  = m_cfg.getDouble("BC_LEFT");
    bcRight                 = m_cfg.getDouble("BC_RIGHT");
    dt                      = m_cfg.getDouble("TIME_STEP");
    tEnd                    = m_cfg.getDouble("END_TIME");
    t                       = m_cfg.getDouble("INITIAL_TIME");

    spatialScheme           = m_cfg.getEnum("SPATIAL_SCHEME",spatialSchemeMap);
    timeScheme              = m_cfg.getEnum("TIME_SCHEME",timeSchemeMap);
    bcType                  = m_cfg.getEnum("BOUNDARY_CONDITIONS",boundaryConditionMap);

    verbose          = m_cfg.getBool("VERBOSE");

    output_file       = m_cfg.getString("OUTPUT_FILE");
    mesh_file        = m_cfg.getString("MESH_FILE");


    std::cout << "Base parameters of the solver initialized" << std::endl;
    
    if (verbose) m_cfg.print();
}


bool BaseSolver::getExplicitSchemeFlag() const {
    return (timeScheme == TimeScheme::EULER_EXPLICIT || 
            timeScheme == TimeScheme::RK4);
}

void BaseSolver::checkStability() const {

    std::cout << "Verifying time stability" << std::endl;

    double cfl  = getCFL();
    double diff = getDiffusionNumber();
    bool flagExplicit = getExplicitSchemeFlag();

    if (flagExplicit && cfl > 1){
        throw InitialStabilityException("CFL",cfl,1.0); //last argument is the limit for this specific stability check
    }
    else {
        std::cout << ( (flagExplicit) ? "Explicit" : "Implicit" ) << " Time scheme" << '\n';
        std::cout << " CFL = " << cfl << " (Limit for Explicit Schemes = 1.0)" << '\n' << '\n';
    }

    std::cout << "Verifying diffusion stability" << std::endl;

    if (flagExplicit && diff > 0.5){
        throw InitialStabilityException("Diffusion number",diff,0.5); //last argument is the limit for this specific stability check
    }
    else {
        std::cout << ( (flagExplicit) ? "Explicit" : "Implicit" ) << " Time scheme" << '\n';
        std::cout << " Diffusion number = " << diff << " (Limit for Explicit Schemes = 0.5)" << '\n' << '\n';
    }

}