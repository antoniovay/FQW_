#include "Lib.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


namespace {
    template <typename T>
    T readValue(const std::string& prompt, T defaultValue) {
        std::cout << prompt << " [по умолчанию " << defaultValue << "]: ";
        std::string line;
        if (!std::getline(std::cin, line))
            return defaultValue;
        if (line.empty())
            return defaultValue;
        std::istringstream input(line);
        T value{};
        if (!(input >> value))
            throw std::invalid_argument("Не удалось прочитать значение: " + prompt);
        return value;
    }
}



int main() {
    
}
