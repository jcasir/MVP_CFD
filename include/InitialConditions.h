// condizioni_iniziali.h
#pragma once
#include <cmath>
#include <numbers> 
#include <ErrorHandler.h>
#include <iostream>


inline double gaussian(double x) {
    constexpr double x0    = 0.5;
    constexpr double sigma = 0.05;
    return std::exp(-std::pow(x - x0, 2) / (2 * sigma * sigma));
}

inline double squareWave(double x) {
    return (x >= 0.3 && x <= 0.7) ? 1.0 : 0.0;
}

inline double sinusoidal(double x) {
    return std::sin(2.0 * std::numbers::pi * x);
}

// inline double hat_func(double x, double y, std::pair<double> rangex, std::pair<double> rangey) {
//     // Defining the limits of the indicies
//     int row_start = static_cast<int>(0.5 / dy);
//     int row_end   = static_cast<int>(1.0 / dy + 1);
//     int col_start = static_cast<int>(0.5 / dx);
//     int col_end   = static_cast<int>(1.0 / dx + 1);

//     // Ciclo per assegnare il valore 2 alla sotto-matrice
//     for (int i = row_start; i < row_end; ++i) {
//         for (int j = col_start; j < col_end; ++j) {
//             u[i][j] = 2.0;
//         }
//     }
//     return std::sin(2.0 * std::numbers::pi * x);
// }

template <typename T>
inline std::function<double(double)> retriveInitialConditionFunction (const T& val){
    switch (val) { 
        case T::GAUSSIAN:
            return gaussian;
        case T::SQUARE_WAVE:
            return squareWave;
        case T::SINUSOIDAL:
            return sinusoidal;
        default:
            throw InvalidOption("Error: Unknown initial condition");
    }
}