#include "Lib.hpp"

// ----- IAtmosphere -----

bool IAtmosphere::contains(double altitude) const {
        return altitude >= getGroundHeight() && altitude <= getTopHeight();
    }


// ----- SimpleAtmosphere -----

double SimpleAtmosphere::getGroundHeight() const {
    return groundHeight_;
}

double SimpleAtmosphere::getTopHeight()    const {
    return topHeight_;
}

double SimpleAtmosphere::getAttenuation(double altitude) const {
    if (!contains(altitude)) return 0.0;
    return attenuation_;
}

double SimpleAtmosphere::getAbsorption(double altitude) const {
    if (!contains(altitude)) return 0.0;
    return absorption_;
}


// ----- LayeredAtmosphere -----

const LayeredAtmosphere::Layer* LayeredAtmosphere::findLayer(double altitude) const {
    for (const auto& layer : layers_) {
        if (altitude >= layer.lowerHeight && altitude < layer.upperHeight)
            return &layer;
    }
    // Если высота ровно на верхней границе последнего слоя
    if (!layers_.empty() &&
        altitude == layers_.back().upperHeight)
        return &layers_.back();
    return nullptr;
}

double LayeredAtmosphere::getGroundHeight() const {
    return groundHeight_;
}

double LayeredAtmosphere::getTopHeight() const {
    return topHeight_;
}

double LayeredAtmosphere::getAttenuation(double altitude) const {
    if (!contains(altitude)) return 0.0;
    const Layer* l = findLayer(altitude);
    return l ? l->attenuation : 0.0;
}

double LayeredAtmosphere::getAbsorption(double altitude) const {
    if (!contains(altitude)) return 0.0;
    const Layer* l = findLayer(altitude);
    return l ? l->absorption : 0.0;
}

const std::vector<LayeredAtmosphere::Layer>& LayeredAtmosphere::getLayers() const { return layers_; }

double LayeredAtmosphere::getLayerNumber(double altitude) const {
    const Layer* l = findLayer(altitude);
    return l ? l->number : -1;
}


// ----- Photon -----

std::mt19937& Photon::rng() {
    static std::mt19937 gen(std::random_device{}());
    return gen;
}

double Photon::uniform(double a, double b) {
    std::uniform_real_distribution<double> d(a, b);
    return d(rng());
}

Vec3 Photon::randomDirection() {
    // Равномерно по ед сфере
    const double u   = uniform(-1.0, 1.0);
    const double phi = uniform(0.0, 2.0 * M_PI);
    const double s   = std::sqrt(std::max(0.0, 1.0 - u*u));
    return { s * std::cos(phi), s * std::sin(phi), u };
}

Vec3 Photon::randomAboveHorizon() {
    // Равномерно по z >= 0
    Vec3 d;
    do { d = randomDirection(); } while (d.z <= 0.0);
    return d;
}

Photon::Photon(const Vec3& position,
               const Vec3& direction,
               double stepLength)
    : position_(position),
      direction_(direction.normalized()),
      stepLength_(stepLength),
      state_(State::Alive),
      pathLength_(0.0),
      stepCount_(0)
{
    if (stepLength_ <= 0.0)
        throw std::invalid_argument("stepLength должен быть > 0");
}

bool Photon::step(const IAtmosphere& atm) {
    if (state_ != State::Alive) return false;

    const double groundZ = atm.getGroundHeight();
    const double topZ    = atm.getTopHeight();

    const Vec3 oldPos = position_;
    Vec3 newPos = oldPos + direction_ * stepLength_;
    ++stepCount_;

    // ---- 1. Вылет за верхнюю границу ----
    if (newPos.z >= topZ) {
        if (direction_.z > 1e-12) {
            double t = (topZ - oldPos.z) / direction_.z;
            t = std::clamp(t, 0.0, stepLength_);
            position_ = oldPos + direction_ * t;
            pathLength_ += t;
        } else {
            position_ = newPos;
            pathLength_ += stepLength_;
        }
        state_ = State::Escaped;
        return false;
    }

    // ---- 2. Удар о землю ----
    if (newPos.z <= groundZ) {
        if (direction_.z < -1e-12) {
            double t = (groundZ - oldPos.z) / direction_.z;
            t = std::clamp(t, 0.0, stepLength_);
            position_ = oldPos + direction_ * t;
            pathLength_ += t;
        } else {
            position_ = newPos;
            position_.z = groundZ;
            pathLength_ += stepLength_;
        }
        direction_ = randomAboveHorizon();
        return true;
    }

    // ---- 3. Внутри атмосферы: розыгрыш ----
    position_ = newPos;
    pathLength_ += stepLength_;

    const double att = atm.getAttenuation(position_.z);
    const double abs = atm.getAbsorption(position_.z);
    const double mu  = att + abs;

    if (mu > 0.0) {
        const double pInteract = 1.0 - std::exp(-mu * stepLength_);
        if (uniform(0.0, 1.0) < pInteract) {
            const double pAbsorb = abs / mu;
            if (uniform(0.0, 1.0) < pAbsorb) {
                state_ = State::Absorbed;
                return false;
            }
            // Рассеяние — изотропно
            direction_ = randomDirection();
        }
    }
    return true;
}


// Simulation

Simulation::Simulation(std::shared_ptr<IAtmosphere> atmosphere,
                       double stepLength)
    : atmosphere_(std::move(atmosphere)),
      stepLength_(stepLength)
{
    if (!atmosphere_)
        throw std::invalid_argument("atmosphere не может быть nullptr");
    if (stepLength_ <= 0.0)
        throw std::invalid_argument("stepLength должен быть > 0");
}

Simulation::Result Simulation::run(const Vec3& startPosition,
                                   const Vec3& startDirection,
                                   int maxSteps) const
{
    Photon photon(startPosition, startDirection, stepLength_);

    while (photon.getState() == Photon::State::Alive &&
           photon.getStepCount() < maxSteps)
    {
        photon.step(*atmosphere_);
    }

    return {
        photon.getState(),
        photon.getPathLength(),
        photon.getStepCount(),
        photon.getPosition(),
        photon.getDirection()
    };
}
