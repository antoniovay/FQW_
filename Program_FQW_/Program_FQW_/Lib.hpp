#ifndef Lib_hpp
#define Lib_hpp

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

constexpr double FQW_PI = 3.141592653589793238462643383279502884;

struct Vec3 {
    double x{}, y{}, z{};
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    double length() const { return std::sqrt(x*x + y*y + z*z); }
    double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y*o.z-z*o.y, z*o.x-x*o.z, x*o.y-y*o.x};
    }
    Vec3 normalized() const {
        const double len = length();
        if (len < 1e-15) throw std::invalid_argument("Нулевой вектор направления");
        return {x/len, y/len, z/len};
    }
    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};

class IAtmosphere {
public:
    virtual ~IAtmosphere() = default;
    virtual double getGroundHeight() const = 0;
    virtual double getTopHeight() const = 0;
    virtual double getAttenuation(double altitude) const = 0;
    virtual double getAbsorption(double altitude) const = 0;
    virtual bool contains(double altitude) const;
    virtual std::vector<double> getBoundaryHeights() const;
};

class SimpleAtmosphere : public IAtmosphere {
protected:
    double groundHeight_, topHeight_, attenuation_, absorption_;
public:
    SimpleAtmosphere(double groundHeight, double topHeight,
                     double attenuation, double absorption);
    double getGroundHeight() const override;
    double getTopHeight() const override;
    double getAttenuation(double altitude) const override;
    double getAbsorption(double altitude) const override;
};

class LayeredAtmosphere : public IAtmosphere {
public:
    struct Layer {
        double lowerHeight;
        double upperHeight;
        double attenuation;
        double absorption;
        double number = 0.0;
    };
protected:
    double groundHeight_, topHeight_;
    std::vector<Layer> layers_;
    const Layer* findLayer(double altitude) const;
public:
    LayeredAtmosphere(double groundHeight, double topHeight, std::vector<Layer> layers);
    double getGroundHeight() const override;
    double getTopHeight() const override;
    double getAttenuation(double altitude) const override;
    double getAbsorption(double altitude) const override;
    const std::vector<Layer>& getLayers() const;
    double getLayerNumber(double altitude) const;
    std::vector<double> getBoundaryHeights() const override;
};

enum class GroundMode { Absorb, DiffuseReflect };

class Photon {
public:
    enum class State { Alive, Absorbed, Escaped, GroundAbsorbed };
    Photon(const Vec3& position, const Vec3& direction);
    const Vec3& getPosition() const { return position_; }
    const Vec3& getDirection() const { return direction_; }
    State getState() const { return state_; }
    double getPathLength() const { return pathLength_; }
    int getStepCount() const { return stepCount_; }
    bool step(const IAtmosphere& atmosphere, double earthRadius, GroundMode groundMode);
    static Vec3 randomIsotropicDirection();
    static Vec3 randomLambertianDirection(const Vec3& outwardNormal);
    static double uniform(double a = 0.0, double b = 1.0);
private:
    Vec3 position_, direction_;
    State state_;
    double pathLength_;
    int stepCount_;
    static double nextSphereIntersection(const Vec3& p, const Vec3& d, double radius);
    static Vec3 sampleScatterDirection();
};

class Simulation {
public:
    struct Result {
        Photon::State state;
        double pathLength;
        int steps;
        Vec3 finalPosition;
        Vec3 finalDirection;
    };
    struct Statistics {
        std::uint64_t launched = 0;
        std::uint64_t escaped = 0;
        std::uint64_t atmosphereAbsorbed = 0;
        std::uint64_t groundAbsorbed = 0;
        std::uint64_t active = 0;
        double meanPathLength = 0.0;
        double escapedFraction() const;
        double atmosphereAbsorbedFraction() const;
        double groundAbsorbedFraction() const;
    };

    explicit Simulation(std::shared_ptr<IAtmosphere> atmosphere,
                        double earthRadius = 6371000.0,
                        GroundMode groundMode = GroundMode::Absorb);
    Result run(const Vec3& startPosition, const Vec3& startDirection,
               int maxEvents = 100000) const;
    Statistics runMany(const Vec3& startPosition, const Vec3& startDirection,
                       std::uint64_t photonCount, int maxEvents = 100000) const;
private:
    std::shared_ptr<IAtmosphere> atmosphere_;
    double earthRadius_;
    GroundMode groundMode_;
};

#endif
