#pragma once

#include <ErrorHandler.h>
#include <string>
#include <unordered_map>
#include <map>
#include <stdexcept>

class ConfigParser {
public:
    explicit ConfigParser(const std::string& filename);

    // Getter con valore di default se la chiave non esiste
    std::string getString(const std::string& key) const;
    int         getInt   (const std::string& key) const;
    double      getDouble(const std::string& key) const;
    bool        getBool  (const std::string& key, bool def = false) const;

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
    void print() const;  // utile per debug

private:
    std::unordered_map<std::string, std::string> data_;

    void parse(const std::string& filename);
    static std::string trim(const std::string& s);
    double parseDouble(const std::string& key, const std::string& value) const;
};