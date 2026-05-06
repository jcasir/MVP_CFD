#pragma once

#include <stdexcept>
#include <string>
#include <iostream>


class GeneralRuntimeError : public std::runtime_error {
public:
    explicit GeneralRuntimeError(const std::string& msg) 
        : std::runtime_error(msg) {}
};

class InitialStabilityException : public GeneralRuntimeError {
public:
    InitialStabilityException(const std::string& msg, double v, double limit) 
        : GeneralRuntimeError("Error: " + msg + " = " + std::to_string(v) + '\n'
            + "Limit for Explicit Schemes = " + std::to_string(limit)) {}
};

class InvalidOption : public GeneralRuntimeError {
public:
    InvalidOption(const std::string& msg) 
        : GeneralRuntimeError(msg) {}
};

class KeyNotFound : public GeneralRuntimeError {
public:
    KeyNotFound(const std::string& key) 
        : GeneralRuntimeError("Key: " + key + " not found in the config file.") {}
};

class InvalidValueForKey : public GeneralRuntimeError {
public:
    InvalidValueForKey(const std::string& val, const std::string& key) 
        : GeneralRuntimeError("Invalid value " + val + " for key: " + key) {}
};

class CannotOpenFile : public GeneralRuntimeError {
public:
    CannotOpenFile(const std::string& filename, const std::string& origin = "") 
        : GeneralRuntimeError("Cannot open " + origin + " file: " + filename) {}
};

class InvalidInitialCondition : public GeneralRuntimeError {
public:
    InvalidInitialCondition(const std::string& ic, const std::string& dim) 
        : GeneralRuntimeError("Initial Condition: " + ic + " not valid for " + dim + " initial condition") {}
};


