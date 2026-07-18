#include "algorithm/ComplexNumerics.h"

#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
std::complex<double> pekerisRoot(std::complex<double> value)
{
    if (value.real() >= 0.0)
    {
        return std::sqrt(value);
    }
    return std::complex<double>(0.0, 1.0) * std::sqrt(-value);
}

std::complex<double> complexSoundSpeed(double soundSpeed,
                                       double attenuation,
                                       double frequency,
                                       char attenuationUnit)
{
    AttenuationContext context;
    context.unit = attenuationUnit;
    context.frequency = frequency;
    context.referenceFrequency = frequency;
    return complexSoundSpeed(0.0, soundSpeed, attenuation, context);
}

double thorpNepersPerMetre(double frequencyHz)
{
    if (!(frequencyHz > 0.0))
    {
        throw std::invalid_argument("frequency must be positive");
    }
    const double f2 = std::pow(frequencyHz / 1000.0, 2.0);
    return (3.3e-3 + 0.11 * f2 / (1.0 + f2) +
            44.0 * f2 / (4100.0 + f2) + 3.0e-4 * f2) /
           8685.8896;
}

double francoisGarrisonNepersPerMetre(
    double frequencyHz, const VolumeAbsorptionParameters &parameters)
{
    if (!(frequencyHz > 0.0) || parameters.salinityPsu < 0.0 ||
        !(parameters.ph > 0.0) ||
        !(parameters.temperatureCelsius > -273.0) ||
        parameters.meanDepthMetres < 0.0)
    {
        throw std::invalid_argument("invalid Francois-Garrison parameters");
    }
    const double f = frequencyHz / 1000.0;
    const double temperature = parameters.temperatureCelsius;
    const double salinity = parameters.salinityPsu;
    const double depth = parameters.meanDepthMetres;
    const double c = 1412.0 + 3.21 * temperature + 1.19 * salinity +
                     0.0167 * depth;
    const double a1 = 8.86 / c * std::pow(10.0, 0.78 * parameters.ph - 5.0);
    const double f1 = 2.8 * std::sqrt(salinity / 35.0) *
                      std::pow(10.0, 4.0 - 1245.0 / (temperature + 273.0));
    const double a2 = 21.44 * salinity / c * (1.0 + 0.025 * temperature);
    const double p2 = 1.0 - 1.37e-4 * depth + 6.2e-9 * depth * depth;
    const double f2 =
        8.17 * std::pow(10.0, 8.0 - 1990.0 / (temperature + 273.0)) /
        (1.0 + 0.0018 * (salinity - 35.0));
    const double p3 = 1.0 - 3.83e-5 * depth + 4.9e-10 * depth * depth;
    const double a3 = temperature < 20.0
                          ? 4.937e-4 - 2.59e-5 * temperature +
                                9.11e-7 * temperature * temperature -
                                1.5e-8 * temperature * temperature * temperature
                          : 3.964e-4 - 1.146e-5 * temperature +
                                1.45e-7 * temperature * temperature -
                                6.5e-10 * temperature * temperature * temperature;
    const double frequencySquared = f * f;
    const double decibelsPerKilometre =
        a1 * f1 * frequencySquared / (f1 * f1 + frequencySquared) +
        a2 * p2 * f2 * frequencySquared / (f2 * f2 + frequencySquared) +
        a3 * p3 * frequencySquared;
    return decibelsPerKilometre / 8685.8896;
}

double biologicalNepersPerMetre(
    double depthMetres, double frequencyHz,
    const VolumeAbsorptionParameters &parameters)
{
    if (!(frequencyHz > 0.0))
    {
        throw std::invalid_argument("frequency must be positive");
    }
    double result = 0.0;
    for (const BiologicalAbsorptionLayer &layer : parameters.biologicalLayers)
    {
        if (layer.bottomDepthMetres < layer.topDepthMetres ||
            !(layer.resonanceFrequencyHz > 0.0) ||
            !(layer.qualityFactor > 0.0) ||
            layer.peakAttenuationDbPerKm < 0.0)
        {
            throw std::invalid_argument("invalid biological absorption layer");
        }
        if (depthMetres >= layer.topDepthMetres &&
            depthMetres <= layer.bottomDepthMetres)
        {
            const double detuning =
                1.0 - std::pow(layer.resonanceFrequencyHz / frequencyHz, 2.0);
            result += layer.peakAttenuationDbPerKm /
                      (detuning * detuning +
                       1.0 / (layer.qualityFactor * layer.qualityFactor)) /
                      8685.8896;
        }
    }
    return result;
}

std::complex<double> complexSoundSpeed(double depthMetres,
                                       double soundSpeed,
                                       double attenuation,
                                       const AttenuationContext &context)
{
    if (!(soundSpeed > 0.0))
    {
        throw std::invalid_argument("sound speed must be positive");
    }
    if (!(context.frequency > 0.0))
    {
        throw std::invalid_argument("frequency must be positive");
    }
    if (attenuation < 0.0)
    {
        throw std::invalid_argument("attenuation must be non-negative");
    }

    const double omega = 2.0 * pi * context.frequency;
    double nepersPerMetre = 0.0;
    switch (context.unit)
    {
    case 'N':
        nepersPerMetre = attenuation;
        break;
    case 'M':
        nepersPerMetre = attenuation / 8.6858896;
        break;
    case 'm':
        if (!(context.referenceFrequency > 0.0) ||
            !(context.transitionFrequency > 0.0))
        {
            throw std::invalid_argument(
                "lowercase-m attenuation requires positive reference and transition frequencies");
        }
        nepersPerMetre = attenuation / 8.6858896;
        if (context.frequency < context.transitionFrequency)
        {
            nepersPerMetre *= std::pow(
                context.frequency / context.referenceFrequency, context.beta);
        }
        else
        {
            nepersPerMetre *=
                (context.frequency / context.referenceFrequency) *
                std::pow(context.transitionFrequency / context.referenceFrequency,
                         context.beta - 1.0);
        }
        break;
    case 'F':
        nepersPerMetre = attenuation * context.frequency / 8685.8896;
        break;
    case 'W':
        nepersPerMetre = attenuation * context.frequency /
                         (8.6858896 * soundSpeed);
        break;
    case 'Q':
        if (attenuation > 0.0)
        {
            nepersPerMetre = omega / (2.0 * soundSpeed * attenuation);
        }
        break;
    case 'L':
        nepersPerMetre = attenuation * omega / soundSpeed;
        break;
    default:
        throw std::invalid_argument("unsupported attenuation unit");
    }

    switch (context.model)
    {
    case OceanAbsorptionModel::None:
        break;
    case OceanAbsorptionModel::Thorpe:
        nepersPerMetre += thorpNepersPerMetre(context.frequency);
        break;
    case OceanAbsorptionModel::FrancGarr:
        nepersPerMetre += francoisGarrisonNepersPerMetre(
            context.frequency, context.volume);
        break;
    case OceanAbsorptionModel::Biological:
        nepersPerMetre += biologicalNepersPerMetre(
            depthMetres, context.frequency, context.volume);
        break;
    }

    const double imaginarySpeed = nepersPerMetre * soundSpeed * soundSpeed / omega;
    if (!std::isfinite(imaginarySpeed) || imaginarySpeed > soundSpeed)
    {
        throw std::invalid_argument("attenuation produces a non-physical complex sound speed");
    }
    return {soundSpeed, imaginarySpeed};
}
}
