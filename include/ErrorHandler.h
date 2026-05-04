#pragma once

#include <stdexcept>
#include <string>
#include <iostream>


class GeneralError : public std::runtime_error {
public:
    explicit GeneralError(const std::string& msg) 
        : std::runtime_error(msg) {}
};

class InitialStabilityException : public GeneralError {
public:
    InitialStabilityException(const std::string& msg, double v, double limit) 
        : GeneralError("Error: " + msg + " = " + std::to_string(v) + '\n'
            + "Limit for Explicit Schemes = " + std::to_string(limit)) {}
};
