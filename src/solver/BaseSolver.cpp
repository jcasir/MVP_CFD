#include "solver/BaseSolver.h"
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

    nx               = m_cfg.getInt("GRID_POINTS");
    L                = m_cfg.getInt("DOMAIN_LENGHT");
    output_freq      = m_cfg.getInt("OUTPUT_FREQUENCY");
    time_iter        = 0;

    c                = m_cfg.getDouble("ADVECTION_SPEED");
    D                = m_cfg.getDouble("DIFFUSION_COEFF");
    bcLeft           = m_cfg.getDouble("BC_LEFT");
    bcRight          = m_cfg.getDouble("BC_RIGHT");
    dt               = m_cfg.getDouble("TIME_STEP");
    tEnd             = m_cfg.getDouble("END_TIME");
    t                = m_cfg.getDouble("INITIAL_TIME");

    spatialScheme    = m_cfg.getEnum("SPATIAL_SCHEME",spatialSchemeMap);
    timeScheme       = m_cfg.getEnum("TIME_SCHEME",timeSchemeMap);
    bcType           = m_cfg.getEnum("BOUNDARY_CONDITIONS",boundaryConditionMap);
    initialCondition = m_cfg.getEnum("INITIAL_CONDITIONS",initialConditionMap);

    verbose          = m_cfg.getBool("VERBOSE");

    outputfile       = m_cfg.getString("OUTPUT_FILE");
    mesh_file        = m_cfg.getString("MESH_FILE");


    std::cout << "Base parameters of the solver initialized" << std::endl;
    
    if (verbose) m_cfg.print();
}