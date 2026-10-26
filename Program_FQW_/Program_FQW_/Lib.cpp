#include "Lib.hpp"
#include <limits>


namespace {
    constexpr double EPS_DISTANCE = 1e-7;
    constexpr double COEFFICIENT_EPS = 1e-15;
    
    void validateCoefficient(double value, const char* name) {
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument(std::string(name) + " должен быть конечным и >= 0");
    }
}


// ----- IAtmosphere -----

bool IAtmosphere::contains(double altitude) const {
    return altitude >= getGroundHeight() && altitude <= getTopHeight();
}

std::vector<double> IAtmosphere::getBoundaryHeights() const {
    return {getGroundHeight(), getTopHeight()};
}


// ----- SimpleAtmosphere -----

SimpleAtmosphere::SimpleAtmosphere(double groundHeight, double topHeight, double attenuation, double absorption)
    : groundHeight_(groundHeight), topHeight_(topHeight), attenuation_(attenuation), absorption_(absorption) {
    if (!std::isfinite(groundHeight_) || !std::isfinite(topHeight_) || topHeight_ <= groundHeight_)
        throw std::invalid_argument("Верхняя граница атмосферы должна быть выше нижней");
    validateCoefficient(attenuation_, "Коэффициент рассеяния");
    validateCoefficient(absorption_, "Коэффициент поглощения");
}

double SimpleAtmosphere::getGroundHeight() const { return groundHeight_; }

double SimpleAtmosphere::getTopHeight() const { return topHeight_; }

double SimpleAtmosphere::getAttenuation(double altitude) const {
    return contains(altitude) ? attenuation_ : 0.0;
}

double SimpleAtmosphere::getAbsorption(double altitude) const {
    return contains(altitude) ? absorption_ : 0.0;
}


// ----- LayeredAtmosphere -----

LayeredAtmosphere::LayeredAtmosphere(double groundHeight, double topHeight, std::vector<Layer> layers)
    : groundHeight_(groundHeight), topHeight_(topHeight), layers_(std::move(layers)) {
    if (!std::isfinite(groundHeight_) || !std::isfinite(topHeight_) || topHeight_ <= groundHeight_)
        throw std::invalid_argument("Верхняя граница атмосферы должна быть выше нижней");
    std::sort(layers_.begin(), layers_.end(), [](const Layer& a, const Layer& b) {
        return a.lowerHeight < b.lowerHeight;
    });
    double previousUpper = groundHeight_;
    for (const auto& layer : layers_) {
        if (!std::isfinite(layer.lowerHeight) || !std::isfinite(layer.upperHeight) || layer.upperHeight <= layer.lowerHeight)
            throw std::invalid_argument("Некорректные границы слоя");
        if (layer.lowerHeight < groundHeight_ || layer.upperHeight > topHeight_)
            throw std::invalid_argument("Слой выходит за границы атмосферы");
        if (layer.lowerHeight < previousUpper)
            throw std::invalid_argument("Слои атмосферы не должны перекрываться");
        validateCoefficient(layer.attenuation, "Коэффициент рассеяния слоя");
        validateCoefficient(layer.absorption, "Коэффициент поглощения слоя");
        previousUpper = layer.upperHeight;
    }
}

const LayeredAtmosphere::Layer* LayeredAtmosphere::findLayer(double altitude) const {
    for (const auto& layer : layers_) {
        if (altitude >= layer.lowerHeight && altitude < layer.upperHeight) return &layer;
    }
    if (!layers_.empty() && altitude == layers_.back().upperHeight)
        return &layers_.back();
    return nullptr;
}

double LayeredAtmosphere::getGroundHeight() const { return groundHeight_; }

double LayeredAtmosphere::getTopHeight() const { return topHeight_; }

double LayeredAtmosphere::getAttenuation(double altitude) const {
    if (!contains(altitude))
        return 0.0;
    const auto* layer = findLayer(altitude);
    return layer ? layer->attenuation : 0.0;
}

double LayeredAtmosphere::getAbsorption(double altitude) const {
    if (!contains(altitude))
        return 0.0;
    const auto* layer = findLayer(altitude);
    return layer ? layer->absorption : 0.0;
}

const std::vector<LayeredAtmosphere::Layer>& LayeredAtmosphere::getLayers() const {
    return layers_;
}

double LayeredAtmosphere::getLayerNumber(double altitude) const {
    const auto* layer = findLayer(altitude);
    return layer ? layer->number : -1.0;
}

std::vector<double> LayeredAtmosphere::getBoundaryHeights() const {
    std::vector<double> boundaries{groundHeight_, topHeight_};
    for (const auto& layer : layers_) {
        boundaries.push_back(layer.lowerHeight);
        boundaries.push_back(layer.upperHeight);
    }
    std::sort(boundaries.begin(), boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
    return boundaries;
}


// ----- Photon -----

Photon::Photon(const Vec3& position, const Vec3& direction)
    : position_(position), direction_(direction.normalized()), state_(State::Alive), pathLength_(0.0), stepCount_(0)
{}

double Photon::uniform(double a, double b) {
    static thread_local std::mt19937_64 generator(std::random_device{}());
    return std::uniform_real_distribution<double>(a, b)(generator);
}

Vec3 Photon::randomIsotropicDirection() {
    const double z = uniform(-1.0, 1.0);
    const double phi = 2.0 * FQW_PI * uniform();
    const double r = std::sqrt(std::max(0.0, 1.0-z*z));
    return {r*std::cos(phi), r*std::sin(phi), z};
}

Vec3 Photon::randomLambertianDirection(const Vec3& outwardNormal) {
    const Vec3 n = outwardNormal.normalized();
    const double cosTheta = std::sqrt(uniform());
    const double sinTheta = std::sqrt(std::max(0.0, 1.0-cosTheta*cosTheta));
    const double phi = 2.0 * FQW_PI * uniform();
    const Vec3 reference = (std::abs(n.z) < 0.9) ? Vec3{0,0,1} : Vec3{1,0,0};
    const Vec3 tangent = reference.cross(n).normalized();
    const Vec3 bitangent = n.cross(tangent).normalized();
    return (n*cosTheta + tangent*(sinTheta*std::cos(phi)) + bitangent*(sinTheta*std::sin(phi))).normalized();
}

double Photon::nextSphereIntersection(const Vec3& p, const Vec3& d, double radius) {
    const double b = p.dot(d);
    const double c = p.dot(p) - radius*radius;
    double discriminant = b*b-c;
    if (discriminant < 0.0) return std::numeric_limits<double>::infinity();
    discriminant = std::sqrt(std::max(0.0, discriminant));
    const double t1 = -b-discriminant;
    const double t2 = -b+discriminant;
    double result = std::numeric_limits<double>::infinity();
    if (t1 > EPS_DISTANCE) result = t1;
    if (t2 > EPS_DISTANCE) result = std::min(result, t2);
    return result;
}

Vec3 Photon::sampleScatterDirection() {
    return randomIsotropicDirection();
}

bool Photon::step(const IAtmosphere& atmosphere, double earthRadius, GroundMode groundMode) {
    if (state_ != State::Alive)
        return false;
    if (!(earthRadius > 0.0) || !std::isfinite(earthRadius))
        throw std::invalid_argument("Радиус Земли должен быть положительным");
    ++stepCount_;
    const double currentRadius = position_.length();
    const double currentAltitude = currentRadius-earthRadius;
    const double radialDirection = position_.dot(direction_);

    if (currentAltitude <= atmosphere.getGroundHeight()+EPS_DISTANCE && radialDirection < 0.0) {
        if (groundMode == GroundMode::Absorb)
            state_ = State::GroundAbsorbed;
        else direction_ = randomLambertianDirection(position_.normalized());
        return state_ == State::Alive;
    }

    const double u = std::max(uniform(), std::numeric_limits<double>::min());
    double remainingOpticalDepth = -std::log(u);
    const auto boundaries = atmosphere.getBoundaryHeights();

    for (int guard = 0; guard < 10000 && state_ == State::Alive; ++guard) {
        const double altitude = position_.length()-earthRadius;
        if (altitude >= atmosphere.getTopHeight()-EPS_DISTANCE && position_.dot(direction_) >= 0.0) {
            state_ = State::Escaped;
            return false;
        }
        if (altitude < atmosphere.getGroundHeight()-EPS_DISTANCE || altitude > atmosphere.getTopHeight()+EPS_DISTANCE) {
            state_ = altitude < atmosphere.getGroundHeight() ? State::GroundAbsorbed : State::Escaped;
            return false;
        }

        double boundaryDistance = std::numeric_limits<double>::infinity();
        for (double h : boundaries) {
            const double distance = nextSphereIntersection(position_, direction_, earthRadius+h);
            if (distance < boundaryDistance) {
                boundaryDistance = distance;
            }
        }
        const double muS = atmosphere.getAttenuation(altitude);
        const double muA = atmosphere.getAbsorption(altitude);
        const double muT = muS + muA;
        const double collisionDistance = (muT > COEFFICIENT_EPS)
            ? remainingOpticalDepth/muT : std::numeric_limits<double>::infinity();

        if (collisionDistance < boundaryDistance) {
            position_ += direction_*collisionDistance;
            pathLength_ += collisionDistance;
            if (uniform() < muA/muT) {
                state_ = State::Absorbed;
                return false;
            }
            direction_ = sampleScatterDirection();
            return true;
        }

        if (!std::isfinite(boundaryDistance)) {
            state_ = State::Escaped;
            return false;
        }
        position_ += direction_*boundaryDistance;
        pathLength_ += boundaryDistance;
        if (muT > COEFFICIENT_EPS)
            remainingOpticalDepth = std::max(0.0, remainingOpticalDepth-muT*boundaryDistance);

        const double hitAltitude = position_.length()-earthRadius;
        if (std::abs(hitAltitude-atmosphere.getTopHeight()) < 1e-4 && position_.dot(direction_) > 0.0) {
            state_ = State::Escaped;
            return false;
        }
        if (std::abs(hitAltitude-atmosphere.getGroundHeight()) < 1e-4 && position_.dot(direction_) < 0.0) {
            if (groundMode == GroundMode::Absorb) {
                state_ = State::GroundAbsorbed;
                return false;
            }
            direction_ = randomLambertianDirection(position_.normalized());
            return true;
        }
        const double nudge = 1e-6;
        position_ += direction_*nudge;
        pathLength_ += nudge;
    }
    return state_ == State::Alive;
}


// ----- Simulation -----

double Simulation::Statistics::escapedFraction() const {
    return launched ? static_cast<double>(escaped)/launched : 0.0;
}

double Simulation::Statistics::atmosphereAbsorbedFraction() const {
    return launched ? static_cast<double>(atmosphereAbsorbed)/launched : 0.0;
}

double Simulation::Statistics::groundAbsorbedFraction() const {
    return launched ? static_cast<double>(groundAbsorbed)/launched : 0.0;
}

Simulation::Simulation(std::shared_ptr<IAtmosphere> atmosphere, double earthRadius, GroundMode groundMode)
    : atmosphere_(std::move(atmosphere)), earthRadius_(earthRadius), groundMode_(groundMode) {
    if (!atmosphere_)
        throw std::invalid_argument("Атмосфера не может быть nullptr");
    if (!(earthRadius_ > 0.0) || !std::isfinite(earthRadius_))
        throw std::invalid_argument("Радиус Земли должен быть положительным");
}

Simulation::Result Simulation::run(const Vec3& startPosition, const Vec3& startDirection, int maxEvents) const {
    if (maxEvents <= 0)
        throw std::invalid_argument("maxEvents должен быть > 0");
    const double altitude = startPosition.length()-earthRadius_;
    if (std::abs(altitude-atmosphere_->getGroundHeight()) > 1e-3)
        throw std::invalid_argument("Стартовая точка должна лежать на поверхности Земли (границе атмосферы)");
    Photon photon(startPosition, startDirection);
    while (photon.getState() == Photon::State::Alive && photon.getStepCount() < maxEvents)
        photon.step(*atmosphere_, earthRadius_, groundMode_);
    return {photon.getState(), photon.getPathLength(), photon.getStepCount(),
            photon.getPosition(), photon.getDirection()};
}

Simulation::Statistics Simulation::runMany(const Vec3& startPosition, const Vec3& startDirection, std::uint64_t photonCount, int maxEvents) const {
    Statistics stats;
    stats.launched = photonCount;
    double totalPath = 0.0;
    for (std::uint64_t i=0; i<photonCount; ++i) {
        const Result result = run(startPosition, startDirection, maxEvents);
        totalPath += result.pathLength;
        switch (result.state) {
            case Photon::State::Escaped: ++stats.escaped; break;
            case Photon::State::Absorbed: ++stats.atmosphereAbsorbed; break;
            case Photon::State::GroundAbsorbed: ++stats.groundAbsorbed; break;
            case Photon::State::Alive: ++stats.active; break;
        }
    }
    stats.meanPathLength = photonCount ? totalPath/static_cast<double>(photonCount) : 0.0;
    return stats;
}
