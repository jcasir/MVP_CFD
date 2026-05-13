// condizioni_iniziali.h
#pragma once
#include "ErrorHandler.hpp"
#include "ConfigParser.hpp"
#include <iostream>
#include <cmath>
#include <numbers> 
#include <vector>




class IInitialCondition1D {
public:
    IInitialCondition1D(const ConfigParser& config) {
        amplitude = config.getDouble("AMPLITUDE_IC");
    }
    virtual void setIC(std::vector<double>& u, const std::vector<double>& x) = 0;
    virtual ~IInitialCondition1D() = default;
protected:
    double amplitude;
};

class IInitialCondition2D {
public:
    IInitialCondition2D(const ConfigParser& config) {
        amplitude = config.getDouble("AMPLITUDE_IC");
    }
    virtual void setIC(std::vector<double>& u, const std::vector<double>& x, const std::vector<double>& y) = 0;
    virtual ~IInitialCondition2D() = default;
protected:
    double amplitude;
};





// 1D Initial conditions
class GaussianIC1D : public IInitialCondition1D { 
public:
    GaussianIC1D(const ConfigParser& config) : IInitialCondition1D(config){
        x0 = config.getDouble("X0_IC");
        sigma = config.getDouble("SIGMA_IC");
    }
    void setIC(std::vector<double>& u, const std::vector<double>& x) override {
        for (size_t i = 0; i < x.size(); ++i) {
            u[i] = amplitude * std::exp(-std::pow(x[i] - x0, 2) / (2 * sigma * sigma));
        }
    }
private:
    double x0;
    double sigma;
};
class SquareWaveIC1D : public IInitialCondition1D{ 
public:
    SquareWaveIC1D(const ConfigParser& config) : IInitialCondition1D(config){
        range_start = config.getDouble("RANGE_START_IC");
        range_end = config.getDouble("RANGE_END_IC");
    }
    void setIC(std::vector<double>& u, const std::vector<double>& x) override {
        for (size_t i = 0; i < x.size(); ++i) {
            u[i] = (x[i] >= range_start && x[i] <= range_end) ? amplitude : 0.0;
        }
    }
private:
    double range_start;
    double range_end;
};
class SinusoidalIC1D : public IInitialCondition1D { 
public:
    SinusoidalIC1D(const ConfigParser& config) : IInitialCondition1D(config){}

    void setIC(std::vector<double>& u, const std::vector<double>& x) override {
        for (size_t i = 0; i < x.size(); ++i) {
            u[i] = amplitude * std::sin(2.0 * std::numbers::pi * x[i]);
        }
    }
};





//2D Initial Conditions
class SquareWaveIC2D : public IInitialCondition2D { 
public:
    SquareWaveIC2D(const ConfigParser& config) : IInitialCondition2D(config){
        range_start_x = config.getDouble("RANGE_START_X_IC");
        range_end_x = config.getDouble("RANGE_END_X_IC");
        range_start_y = config.getDouble("RANGE_START_Y_IC");
        range_end_y = config.getDouble("RANGE_END_Y_IC");
        baseline = config.getDouble("BASELINE_VALUE_IC");
    }
    void setIC(std::vector<double>& u, const std::vector<double>& x, const std::vector<double>& y) override {
        for (size_t i = 0; i < x.size(); ++i) {
            for (size_t j = 0; j < y.size(); ++j){
                u[i * y.size() + j] = (x[i] >= range_start_x && x[i] <= range_end_x
                                    && y[j] >= range_start_y && y[j] <= range_end_y ) ? amplitude : baseline;
            }
        }
    }
private:
    double range_start_x;
    double range_end_x;
    double range_start_y;
    double range_end_y;
    double baseline;
};


// factory functions to instanciate the correct class.
inline std::unique_ptr<IInitialCondition1D> makeIC1D(const ConfigParser& config) {
    std::string ic = config.getString("INITIAL_CONDITIONS");
    if      (ic == "GAUSSIAN")    return std::make_unique<GaussianIC1D>(config);
    else if (ic == "SQUARE_WAVE") return std::make_unique<SquareWaveIC1D>(config);
    else if (ic == "SINUSOIDAL")  return std::make_unique<SinusoidalIC1D>(config);
    else throw InvalidInitialCondition(ic,"1D");
}
inline std::unique_ptr<IInitialCondition2D> makeIC2D(const ConfigParser& config){
    std::string ic = config.getString("INITIAL_CONDITIONS");
    if (ic == "SQUARE_WAVE") return std::make_unique<SquareWaveIC2D>(config);
    else throw InvalidInitialCondition(ic,"2D");
}




