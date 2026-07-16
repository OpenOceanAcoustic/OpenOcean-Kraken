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
    if (!(soundSpeed > 0.0))
    {
        throw std::invalid_argument("sound speed must be positive");
    }
    if (!(frequency > 0.0))
    {
        throw std::invalid_argument("frequency must be positive");
    }
    if (attenuation < 0.0)
    {
        throw std::invalid_argument("attenuation must be non-negative");
    }

    const double omega = 2.0 * pi * frequency;
    double nepersPerMetre = 0.0;
    switch (attenuationUnit)
    {
    case 'N':
        nepersPerMetre = attenuation;
        break;
    case 'M':
    case 'm':
        nepersPerMetre = attenuation / 8.6858896;
        break;
    case 'F':
        nepersPerMetre = attenuation * frequency / 8685.8896;
        break;
    case 'W':
        nepersPerMetre = attenuation * frequency / (8.6858896 * soundSpeed);
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

    const double imaginarySpeed = nepersPerMetre * soundSpeed * soundSpeed / omega;
    if (!std::isfinite(imaginarySpeed) || imaginarySpeed > soundSpeed)
    {
        throw std::invalid_argument("attenuation produces a non-physical complex sound speed");
    }
    return {soundSpeed, imaginarySpeed};
}
}
