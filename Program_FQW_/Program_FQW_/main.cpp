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
    try {
        std::cout << "Моделирование переноса фотонов в сферической атмосфере Земли (Монте-Карло)\n";
        
        const std::uint64_t photonCount = readValue<std::uint64_t>("Число фотонов", 10000);
        const int mediumType = readValue<int>("Среда: 1 — однородная, 2 — слоистая", 1);
        const int groundType = readValue<int>("Поверхность Земли: 1 — поглощает, 2 — диффузно отражает", 1);
        const double atmosphereHeightKm = readValue<double>("Высота атмосферы, км", 100.0);
        const double latitudeDeg = readValue<double>("Широта стартовой точки, градусы", 0.0);
        const double longitudeDeg = readValue<double>("Долгота стартовой точки, градусы", 0.0);
        const double scattering = readValue<double>("Базовый коэффициент рассеяния mu_s, 1/м", 1.0e-5);
        const double absorption = readValue<double>("Базовый коэффициент поглощения mu_a, 1/м", 2.0e-6);

        if (photonCount == 0) throw std::invalid_argument("Число фотонов должно быть > 0");
        if (atmosphereHeightKm <= 0.0 || std::abs(latitudeDeg) > 90.0 ||
            !std::isfinite(longitudeDeg) || (mediumType != 1 && mediumType != 2) ||
            (groundType != 1 && groundType != 2))
            throw std::invalid_argument("Проверьте тип среды, режим поверхности, высоту и координаты");

        constexpr double earthRadius = 6371000.0;
        const double topHeight = atmosphereHeightKm*1000.0;
        
        std::shared_ptr<IAtmosphere> atmosphere;
        
        if (mediumType == 2) {
            std::vector<LayeredAtmosphere::Layer> layers{
                {0.0, std::min(10000.0, topHeight), scattering*3.0, absorption*2.0, 1.0},
                {std::min(10000.0, topHeight), std::min(30000.0, topHeight), scattering, absorption, 2.0},
                {std::min(30000.0, topHeight), topHeight, scattering*0.1, absorption*0.1, 3.0}
            };
            layers.erase(std::remove_if(layers.begin(), layers.end(), [](const auto& layer) {
                return layer.upperHeight <= layer.lowerHeight;
            }), layers.end());
            atmosphere = std::make_shared<LayeredAtmosphere>(0.0, topHeight, std::move(layers));
        }
        else {
            atmosphere = std::make_shared<SimpleAtmosphere>(0.0, topHeight, scattering, absorption);
        }
 
        const double lat = latitudeDeg*FQW_PI/180.0;
        const double lon = longitudeDeg*FQW_PI/180.0;
        
        const Vec3 normal{std::cos(lat)*std::cos(lon), std::cos(lat)*std::sin(lon), std::sin(lat)};
        const Vec3 startPosition = normal*earthRadius;
        const GroundMode groundMode = (groundType == 2) ? GroundMode::DiffuseReflect : GroundMode::Absorb;
        
        Simulation simulation(atmosphere, earthRadius, groundMode);
        
        const auto stats = simulation.runMany(startPosition, normal, photonCount);

        std::cout << "\n\n------------------- Результаты -------------------\n" << std::fixed << std::setprecision(6);
        std::cout << "Запущено фотонов:             " << stats.launched << '\n';
        std::cout << "Вышло из атмосферы:           " << stats.escaped << " (" << 100.0*stats.escapedFraction() << " %)\n";
        std::cout << "Поглощено атмосферой:         " << stats.atmosphereAbsorbed << " (" << 100.0*stats.atmosphereAbsorbedFraction() << " %)\n";
        std::cout << "Поглощено поверхностью Земли: " << stats.groundAbsorbed << " (" << 100.0*stats.groundAbsorbedFraction() << " %)\n";
        std::cout << "Не завершено за лимит событий: " << stats.active << '\n';
        std::cout << "Средняя длина пути, м:        " << stats.meanPathLength << '\n';
        std::cout << "Баланс: " << (stats.escaped + stats.atmosphereAbsorbed + stats.groundAbsorbed + stats.active) << " / " << stats.launched << '\n';
        std::cout << "--------------------------------------------------";
        return 0;
    }
    
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << '\n';
        return 1;
    }
}
