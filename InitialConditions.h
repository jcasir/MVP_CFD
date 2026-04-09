// condizioni_iniziali.h
#pragma once
#include <cmath>
#include <numbers> 
#include <stdexcept>

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
            std::cerr << "Errore: Condizione iniziale sconosciuta!" << std::endl;
            throw std::invalid_argument("Condizione iniziale sconosciuta");
    }
}