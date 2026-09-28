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

// Полезно: узнать имя слоя на высоте
std::string LayeredAtmosphere::getLayerName(double altitude) const {
    const Layer* l = findLayer(altitude);
    return l ? l->name : "вне атмосферы";
}
