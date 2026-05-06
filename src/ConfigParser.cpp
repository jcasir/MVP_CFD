#include "ConfigParser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

ConfigParser::ConfigParser(const std::string& filename) {
    parse(filename);
}

void ConfigParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw CannotOpenFile(filename);

    std::string line;
    while (std::getline(file, line)) {
        // Rimuovi commenti (tutto dopo '#')
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
            line = line.substr(0, comment_pos);

        line = trim(line);
        if (line.empty()) continue;

        // Formato atteso: KEY = VALUE
        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key   = trim(line.substr(0, eq_pos));
        std::string value = trim(line.substr(eq_pos + 1));

        if (!key.empty())
            data_[key] = value;
    }
}

std::string ConfigParser::getString(const std::string& key) const {
    auto it = data_.find(key);
    if (it == data_.end()) throw KeyNotFound(key);
    return it->second;
}

int ConfigParser::getInt(const std::string& key) const {
    auto it = data_.find(key);
    if (it == data_.end()) throw KeyNotFound(key);
    try { return std::stoi(it->second); }
    catch (...) { throw std::runtime_error("Invalid int for key: " + key); }
}

double ConfigParser::parseDouble(const std::string& key, const std::string& value) const {
    try { return std::stod(value); }
    catch (...) { throw std::runtime_error("Invalid double for key: " + key); }
}

double ConfigParser::getDouble(const std::string& key) const {
    auto it = data_.find(key);
    if (it == data_.end()) throw KeyNotFound(key);
    return parseDouble(key, it->second);
}

bool ConfigParser::getBool(const std::string& key, bool def) const {
    auto it = data_.find(key);
    if (it == data_.end()) return def;
    std::string v = it->second;
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);
    return (v == "true" || v == "1" || v == "yes");
}

bool ConfigParser::hasKey(const std::string& key) const {
    return data_.count(key) > 0;
}

void ConfigParser::print() const {
    for (const auto& [k, v] : data_)
        std::cout << k << " = " << v << "\n";
    std::cout << '\n';
}

std::string ConfigParser::trim(const std::string& s) {
    const char* ws = " \t\r\n";
    size_t start = s.find_first_not_of(ws);
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(ws);
    return s.substr(start, end - start + 1);
}