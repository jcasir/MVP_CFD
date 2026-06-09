#pragma once

#include <ErrorHandler.hpp>
#include <string>
#include <sstream>
#include <array>
#include <unordered_map>
#include <map>
#include <stdexcept>

class ConfigParser {
public:
    explicit ConfigParser(const std::string& filename);

    // Getters — throw KeyNotFound if the key is missing
    std::string getString(const std::string& key) const;
    int         getInt   (const std::string& key) const;
    double      getDouble(const std::string& key) const;
    bool        getBool  (const std::string& key, bool def = false) const;

    //Getter for the BCs for the Burgers and IncNS solvers
    template <size_t N>
    std::array<double,N> getBCs(const std::string& key) const{
        auto it = data_.find(key);
        if (it == data_.end()) throw KeyNotFound(key);

        std::stringstream is(it->second);
        char t; // Variable to throw away '{'; ',' and '}';
        std::array<double,N> results;

        // Start reading the string
        try{
            is >> t; // '{';
            is >> results[0]; // Fencepost problem
            for (size_t i = 1; i < N; ++i){
                is >> t; // ','
                is >> results[i];
            }
            is >> t; // '}';

        }
        // Catch-all block for any unexpected system exceptions during execution
        catch (...) {throw std::runtime_error("Invalid double for key: " + key); }

        // Check if parsing failed (e.g., missing elements, bad formatting, or wrong data types)
        if (is.fail()){
            throw std::runtime_error("Parsing of the Boundary Conditions values failed");
        }
        
        return results;
    }

    // Getter for the enumerations
    template <typename T>
    T getEnum(const std::string& key, const std::map<std::string, T>& mapping) const {
        const std::string& val = getString(key);
        auto it = mapping.find(val);
        if (it == mapping.end())
            throw InvalidValueForKey(val,key);
        return it->second;
    }

    bool hasKey(const std::string& key) const;
    void print() const;  // useful for debug

private:
    std::unordered_map<std::string, std::string> data_;

    void parse(const std::string& filename);
    static std::string trim(const std::string& s);
    double parseDouble(const std::string& key, const std::string& value) const;
};