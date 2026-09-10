#ifndef SYMBOLOGY_NUMERIC_PARSE_HPP
#define SYMBOLOGY_NUMERIC_PARSE_HPP
#include <string>
#include <stdexcept>

inline bool decimal_integer(const std::string& value) {
    size_t i = !value.empty() && (value[0] == '-' || value[0] == '+') ? 1 : 0;
    if (i == value.size()) return false;
    for (; i < value.size(); ++i)
        if (value[i] < '0' || value[i] > '9') return false;
    return true;
}

inline long long strict_integer(const std::string& value) {
    if (!decimal_integer(value)) throw std::runtime_error("expected an integer: " + value);
    return std::stoll(value);
}

inline void validate_rational(const std::string& value) {
    const auto slash = value.find('/');
    if (!decimal_integer(value.substr(0, slash)))
        throw std::runtime_error("expected an exact rational: " + value);
    if (slash != std::string::npos) {
        const auto den = value.substr(slash + 1);
        if (!decimal_integer(den) || den.find_first_not_of("+-0") == std::string::npos)
            throw std::runtime_error("invalid or zero rational denominator: " + value);
    }
}
#endif
