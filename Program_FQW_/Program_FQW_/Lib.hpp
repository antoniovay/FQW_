
#ifndef Lib_hpp
#define Lib_hpp

#include <stdio.h>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
#include <algorithm>


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
        std::string name;
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
                throw std::invalid_argument("Слой '" + l.name +
                    "': верхняя граница ниже нижней");
            if (l.lowerHeight < groundHeight_ || l.upperHeight > topHeight_)
                throw std::invalid_argument("Слой '" + l.name +
                    "' выходит за границы атмосферы");
        }
    }

    double getGroundHeight() const override;
    double getTopHeight()    const override;
    double getAttenuation(double altitude) const override;
    double getAbsorption(double altitude) const override;

    const std::vector<Layer>& getLayers() const;
    
    std::string getLayerName(double altitude) const;
};


#endif /* Lib_hpp */
