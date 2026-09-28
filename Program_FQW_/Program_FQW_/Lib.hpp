#ifndef Lib_hpp
#define Lib_hpp

#include <cstdio>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <random>



// 3D вектор

struct Vec3 {
    double x{}, y{}, z{};

    Vec3() = default;
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    double length() const { return std::sqrt(x*x + y*y + z*z); }

    Vec3 normalized() const {
        double len = length();
        if (len < 1e-12) return {0.0, 0.0, 1.0};
        return {x/len, y/len, z/len};
    }

    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s)  const { return {x*s, y*s, z*s}; }
};


// Интерфейс

class IAtmosphere {
public:
    virtual ~IAtmosphere() = default;

    virtual double getGroundHeight() const = 0;
    virtual double getTopHeight() const = 0;
    virtual double getAttenuation(double altitude) const = 0; // Ослабление
    virtual double getAbsorption(double altitude) const = 0; // Поглощение

    virtual bool contains(double altitude) const; // Проверка на содержание
};


// Простая однородная атмосфера

class SimpleAtmosphere : public IAtmosphere {
protected:
    double groundHeight_;
    double topHeight_;
    double attenuation_;
    double absorption_;

public:
    SimpleAtmosphere(double groundHeight,
                     double topHeight,
                     double attenuation,
                     double absorption)
        : groundHeight_(groundHeight),
          topHeight_(topHeight),
          attenuation_(attenuation),
          absorption_(absorption)
    {
        if (topHeight_ <= groundHeight_)
            throw std::invalid_argument("Верхняя граница должна быть выше нижней");
    }

    double getGroundHeight() const override;
    double getTopHeight() const override;
    double getAttenuation(double altitude) const override;
    double getAbsorption(double altitude) const override;
};


// Слоистая атмосфера

class LayeredAtmosphere : public IAtmosphere {
public:
    struct Layer {
        double lowerHeight;
        double upperHeight;
        double attenuation;
        double absorption;
        double number;
    };

protected:
    double groundHeight_;
    double topHeight_;
    std::vector<Layer> layers_;

    const Layer* findLayer(double altitude) const;

public:
    LayeredAtmosphere(double groundHeight,
                      double topHeight,
                      std::vector<Layer> layers)
        : groundHeight_(groundHeight),
          topHeight_(topHeight),
          layers_(std::move(layers))
    {
        if (topHeight_ <= groundHeight_)
            throw std::invalid_argument("Верхняя граница должна быть выше нижней");

        // Сортируем слои по высоте и валидируем
        std::sort(layers_.begin(), layers_.end(),
                  [](const Layer& a, const Layer& b) {
                      return a.lowerHeight < b.lowerHeight;
                  });

        for (const auto& l : layers_) {
            if (l.upperHeight <= l.lowerHeight)
                throw std::invalid_argument("Слой " + std::to_string(l.number) +
                    ": верхняя граница ниже нижней");
            if (l.lowerHeight < groundHeight_ || l.upperHeight > topHeight_)
                throw std::invalid_argument("Слой " + std::to_string(l.number) +
                    " выходит за границы атмосферы");
        }
    }

    double getGroundHeight() const override;
    double getTopHeight()    const override;
    double getAttenuation(double altitude) const override;
    double getAbsorption(double altitude) const override;

    const std::vector<Layer>& getLayers() const;
    
    double getLayerNumber(double altitude) const;
};


// Фотон

class Photon {
public:
    enum class State { Alive, Absorbed, Escaped };

    Photon(const Vec3& position, const Vec3& direction, double stepLength);

    const Vec3& getPosition()  const { return position_; }
    const Vec3& getDirection() const { return direction_; }
    State getState() const { return state_; }
    double getPathLength() const { return pathLength_; }
    int getStepCount() const { return stepCount_; }

    // Сделать один шаг. False, если фотон больше не Alive.
    bool step(const IAtmosphere& atm);

private:
    Vec3   position_;
    Vec3   direction_;
    double stepLength_;
    State  state_;
    double pathLength_;
    int    stepCount_;

    static std::mt19937& rng();
    static double uniform(double a, double b);
    static Vec3 randomDirection();      // изотропно по сфере
    static Vec3 randomAboveHorizon();   // изотропно по верхней полусфере
};


// Симуляция

class Simulation {
public:
    struct Result {
        Photon::State state;
        double pathLength;
        int steps;
        Vec3 finalPosition;
        Vec3 finalDirection;
    };

    Simulation(std::shared_ptr<IAtmosphere> atmosphere,
               double stepLength);

    // startPosition (z == getGroundHeight()).
    Result run(const Vec3& startPosition,
               const Vec3& startDirection,
               int maxSteps = 100000) const;

private:
    std::shared_ptr<IAtmosphere> atmosphere_;
    double stepLength_;
};



#endif /* Lib_hpp */
