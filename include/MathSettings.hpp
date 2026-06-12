#pragma once

#include <string>
#include <map>

enum class SpatialScheme {
    CENTRAL,        // Central difference scheme (2nd order)
    UPWIND,         // Upwind scheme (1st order)
    QUICK           // QUICK scheme (3rd order)
};

enum class TimeScheme {
    EULER_EXPLICIT,     // Explicit Euler
    RK4                 // 4th order Runge-Kutta
};

enum class BoundaryCondition {
    DIRICHLET,      // Fixed value at boundary
    NEUMANN,        // Fixed derivative at boundary
    PERIODIC,       // Periodic conditions
    MIXED           // For Incompressible NS 2D when different types of BCs are needed for different sides of the domain
};

//enum maps initialization
inline const std::map<std::string, SpatialScheme> spatialSchemeMap = {
    {"CENTRAL", SpatialScheme::CENTRAL},
    {"UPWIND",  SpatialScheme::UPWIND},
    {"QUICK",   SpatialScheme::QUICK}
};
inline const std::map<std::string, TimeScheme> timeSchemeMap = {
    {"EULER_EXPLICIT", TimeScheme::EULER_EXPLICIT},
    {"RK4",   TimeScheme::RK4}
};
inline const std::map<std::string, BoundaryCondition> boundaryConditionMap = {
    {"DIRICHLET", BoundaryCondition::DIRICHLET},
    {"NEUMANN",  BoundaryCondition::NEUMANN},
    {"PERIODIC",   BoundaryCondition::PERIODIC},
    {"MIXED",   BoundaryCondition::MIXED}
};