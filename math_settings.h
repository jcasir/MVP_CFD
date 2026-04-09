#pragma once

#include <string>
#include <map>

enum class SpatialScheme {
    CENTRAL,        // Schema centrato (2° ordine)
    UPWIND,         // Schema upwind (1° ordine)
    QUICK           // Schema QUICK (3° ordine)
};

enum class TimeScheme {
    EULER_EXPLICIT,     // Eulero esplicito
    RK4                 // Runge-Kutta 4° ordine
};

enum class BoundaryCondition {
    DIRICHLET,      // Valore fisso al contorno
    NEUMANN,        // Derivata fissa al contorno
    PERIODIC        // Condizioni periodiche
};

enum class InitialCondition {
    GAUSSIAN,      // Bell-shaped smooth curve
    SQUARE_WAVE,        // Alternating high-low signal
    SINUSOIDAL        // Smooth periodic oscillation
};

//enum maps initialisation
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
    {"PERIODIC",   BoundaryCondition::PERIODIC}
};
inline const std::map<std::string, InitialCondition> initialConditionMap = {
    {"GAUSSIAN", InitialCondition::GAUSSIAN},
    {"SQUARE_WAVE",  InitialCondition::SQUARE_WAVE},
    {"SINUSOIDAL",   InitialCondition::SINUSOIDAL}
};