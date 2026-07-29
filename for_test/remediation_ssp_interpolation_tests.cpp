#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexNumerics.h"
#include "algorithm/KrakencSolver.h"
#include "algorithm/SspInterpolation.h"

#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>

using namespace OpenOceanKrakenc;

namespace
{
void require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

AcousticCase twoNodeElastic(AcousticInterpolation interpolation)
{
    AcousticCase input;
    input.frequency = 100.0;
    input.referenceFrequency = 100.0;
    input.attenuationUnit = 'W';
    input.interpolation = interpolation;
    AcousticLayer layer;
    layer.baseMesh = 2;
    layer.topDepth = 0.0;
    layer.bottomDepth = 100.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.2, 0.1},
        {100.0, 1800.0, 900.0, 2.0, 0.8, 0.5},
    };
    input.layers.push_back(layer);
    return input;
}

bool nearComplex(std::complex<double> actual,
                 std::complex<double> reference)
{
    return std::abs(actual - reference) <=
           1.0e-11 + 5.0e-13 * std::abs(reference);
}

AttenuationContext lossPerWavelengthContext()
{
    AttenuationContext context;
    context.unit = 'L';
    context.frequency = 100.0;
    context.referenceFrequency = 100.0;
    return context;
}

AcousticLayer threeNodeElasticLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 2;
    layer.topDepth = 5.0;
    layer.bottomDepth = 15.0;
    layer.samples = {
        {5.0, 1000.0, 500.0, 1.0, 10.0 / 1000.0, 5.0 / 500.0},
        {10.0, 1200.0, 600.0, 1.5, 30.0 / 1200.0, 15.0 / 600.0},
        {15.0, 1600.0, 800.0, 2.0, 70.0 / 1600.0, 35.0 / 800.0},
    };
    return layer;
}

AcousticLayer twoNodeFluidLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 2;
    layer.topDepth = 0.0;
    layer.bottomDepth = 10.0;
    layer.samples = {
        {0.0, 1000.0, 0.0, 1.0, 10.0 / 1000.0, 0.0},
        {10.0, 2000.0, 0.0, 2.0, 20.0 / 2000.0, 0.0},
    };
    return layer;
}

AcousticLayer twoNodeSplineLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 2;
    layer.topDepth = 0.0;
    layer.bottomDepth = 1.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {1.0, 1600.0, 800.0, 2.0, 2.0 / 1600.0, 4.0 / 800.0},
    };
    return layer;
}

AcousticLayer threeNodeSplineLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 2;
    layer.topDepth = 0.0;
    layer.bottomDepth = 2.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {1.0, 1501.0, 701.0, 2.0, 1.0 / 1501.0, 2.0 / 701.0},
        {2.0, 1504.0, 704.0, 5.0, 4.0 / 1504.0, 8.0 / 704.0},
    };
    return layer;
}

AcousticLayer fourNodeSplineLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 3;
    layer.topDepth = 0.0;
    layer.bottomDepth = 3.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {1.0, 1501.0, 701.0, 2.0, 1.0 / 1501.0, 2.0 / 701.0},
        {2.0, 1508.0, 708.0, 9.0, 8.0 / 1508.0, 16.0 / 708.0},
        {3.0, 1527.0, 727.0, 28.0, 27.0 / 1527.0, 54.0 / 727.0},
    };
    return layer;
}

AcousticLayer nonUniformFourNodeSplineLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 4;
    layer.topDepth = 0.0;
    layer.bottomDepth = 4.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {0.5, 1500.125, 700.125, 1.125,
         0.125 / 1500.125, 0.25 / 700.125},
        {2.0, 1508.0, 708.0, 9.0,
         8.0 / 1508.0, 16.0 / 708.0},
        {4.0, 1564.0, 764.0, 65.0,
         64.0 / 1564.0, 128.0 / 764.0},
    };
    return layer;
}

AcousticLayer threeNodePchipLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 4;
    layer.topDepth = 0.0;
    layer.bottomDepth = 2.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {1.0, 1600.0, 800.0, 2.0, 0.0, 0.0},
        {2.0, 1600.0, 800.0, 2.0,
         1.0 / 1600.0, 2.0 / 800.0},
    };
    return layer;
}

AcousticLayer fourNodePchipLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 6;
    layer.topDepth = 0.0;
    layer.bottomDepth = 3.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {1.0, 1501.0, 701.0, 2.0,
         1.0 / 1501.0, 2.0 / 701.0},
        {2.0, 1504.0, 704.0, 5.0, 0.0, 0.0},
        {3.0, 1509.0, 709.0, 10.0,
         1.0 / 1509.0, 2.0 / 709.0},
    };
    return layer;
}

AcousticLayer threeNodePchipEndpointCapLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 4;
    layer.topDepth = 0.0;
    layer.bottomDepth = 2.0;
    layer.samples = {
        {0.0, 1500.0, 719.0, 1.0,
         10.0 / 1500.0, 19.0 / 719.0},
        {1.0, 1501.0, 709.0, 2.0,
         11.0 / 1501.0, 9.0 / 709.0},
        {2.0, 1491.0, 710.0, 1.0,
         1.0 / 1491.0, 10.0 / 710.0},
    };
    return layer;
}

AcousticLayer nonUniformFourNodePchipProjectionLayer()
{
    AcousticLayer layer;
    layer.baseMesh = 8;
    layer.topDepth = 0.0;
    layer.bottomDepth = 4.0;
    layer.samples = {
        {0.0, 1500.0, 700.0, 1.0, 0.0, 0.0},
        {1.0, 1501.0, 701.0, 2.0,
         1.0 / 1501.0, 2.0 / 701.0},
        {3.0, 1503.0, 703.0, 4.0,
         3.0 / 1503.0, 6.0 / 703.0},
        {4.0, 1603.0, 803.0, 104.0,
         103.0 / 1603.0, 206.0 / 803.0},
    };
    return layer;
}

AcousticCase matrixCase(AcousticLayer layer,
                        AcousticInterpolation interpolation)
{
    AcousticCase input;
    input.frequency = 100.0;
    input.referenceFrequency = 100.0;
    input.attenuationUnit = 'L';
    input.interpolation = interpolation;
    input.layers.push_back(std::move(layer));
    return input;
}

void requireInvalidArgument(const std::function<void()> &operation,
                            const std::string &message)
{
    try
    {
        operation();
    }
    catch (const std::invalid_argument &)
    {
        return;
    }
    throw std::runtime_error(message);
}

std::string invalidArgumentMessage(
    const std::function<void()> &operation,
    const std::string &message)
{
    try
    {
        operation();
    }
    catch (const std::invalid_argument &error)
    {
        return error.what();
    }
    throw std::runtime_error(message);
}

void requireContains(const std::string &actual,
                     const std::string &fragment,
                     const std::string &message)
{
    require(actual.find(fragment) != std::string::npos, message);
}

void testMatrixRetainsDuplicateInterfaceNodes()
{
    AcousticCase input;
    input.frequency = 100.0;
    input.referenceFrequency = 100.0;
    input.attenuationUnit = 'L';
    input.interpolation = AcousticInterpolation::CLinear;

    AcousticLayer upper;
    upper.baseMesh = 1;
    upper.topDepth = 0.0;
    upper.bottomDepth = 100.0;
    upper.samples = {
        {0.0, 1400.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 1.0 / 1500.0, 0.0},
    };
    AcousticLayer lower;
    lower.baseMesh = 1;
    lower.topDepth = 100.0;
    lower.bottomDepth = 200.0;
    lower.samples = {
        {100.0, 2500.0, 0.0, 2.0, 10.0 / 2500.0, 0.0},
        {200.0, 2600.0, 0.0, 2.0, 12.0 / 2600.0, 0.0},
    };
    input.layers = {upper, lower};

    const AcousticMatrix matrix = buildAcousticMatrix(input, 1);
    require(matrix.depth ==
                std::vector<double>{0.0, 100.0, 100.0, 200.0},
            "ALG-022 matrix must retain both 100 m interface nodes");
    require(nearComplex(matrix.cp[1], {1500.0, 1.0}),
            "ALG-022 upper 100 m node must remain upper-medium data");
    require(nearComplex(matrix.cp[2], {2500.0, 10.0}),
            "ALG-022 lower 100 m node must remain lower-medium data");
}

AcousticLayer localPolynomialLayer(double topDepth, bool upper)
{
    AcousticLayer layer;
    layer.baseMesh = 2;
    layer.topDepth = topDepth;
    layer.bottomDepth = topDepth + 100.0;
    for (int index = 0; index < 4; ++index)
    {
        const double coordinate = static_cast<double>(index) / 3.0;
        const double coordinate2 = coordinate * coordinate;
        const double coordinate3 = coordinate2 * coordinate;
        AcousticSample sample;
        sample.depth = topDepth + 100.0 * coordinate;
        if (upper)
        {
            sample.cp = 1400.0 + 100.0 * coordinate3;
            const double cpImaginary = coordinate2;
            sample.cs = 600.0 + 100.0 * coordinate2;
            const double csImaginary = 2.0 * coordinate3;
            sample.alphaP = cpImaginary / sample.cp;
            sample.alphaS = csImaginary / sample.cs;
            sample.rho = 1.0;
        }
        else
        {
            sample.cp = 2500.0 + 100.0 * coordinate2;
            const double cpImaginary = 10.0 + 2.0 * coordinate3;
            sample.cs = 800.0 + 100.0 * coordinate3;
            const double csImaginary = 20.0 + 4.0 * coordinate2;
            sample.alphaP = cpImaginary / sample.cp;
            sample.alphaS = csImaginary / sample.cs;
            sample.rho = 2.0;
        }
        layer.samples.push_back(sample);
    }
    return layer;
}

void testPreparedLayersUseOnlyLocalPolynomialNodes()
{
    const AttenuationContext context = lossPerWavelengthContext();
    const PreparedComplexSspLayer upper(
        localPolynomialLayer(0.0, true),
        AcousticInterpolation::CubicSpline, context, 1);
    const PreparedComplexSspLayer lower(
        localPolynomialLayer(100.0, false),
        AcousticInterpolation::CubicSpline, context, 2);

    const InterpolatedMaterial upperMidpoint = upper.evaluate(50.0);
    require(nearComplex(upperMidpoint.cp, {1412.5, 0.25}) &&
                nearComplex(upperMidpoint.cs, {625.0, 0.25}),
            "ALG-022 upper-medium P/S midpoint must use upper nodes only");
    const InterpolatedMaterial lowerMidpoint = lower.evaluate(150.0);
    require(nearComplex(lowerMidpoint.cp, {2525.0, 10.25}) &&
                nearComplex(lowerMidpoint.cs, {812.5, 21.0}),
            "ALG-022 lower-medium P/S midpoint must use lower nodes only");
}

void testSampleEndpointToleranceIsStrict()
{
    const double tolerance =
        100.0 * std::numeric_limits<float>::epsilon();
    AcousticLayer layer = twoNodeFluidLayer();
    layer.samples.front().depth = tolerance;
    requireInvalidArgument(
        [&layer]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::CLinear,
                lossPerWavelengthContext(), 1));
        },
        "ALG-022 endpoint delta equal to tolerance must be rejected");

    layer.samples.front().depth =
        std::nextafter(tolerance, 0.0);
    static_cast<void>(PreparedComplexSspLayer(
        layer, AcousticInterpolation::CLinear,
        lossPerWavelengthContext(), 1));

    layer = twoNodeFluidLayer();
    layer.topDepth = -10.0;
    layer.bottomDepth = 0.0;
    layer.samples.front().depth = -10.0;
    layer.samples.back().depth = -tolerance;
    requireInvalidArgument(
        [&layer]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::CLinear,
                lossPerWavelengthContext(), 1));
        },
        "ALG-022 bottom delta equal to tolerance must be rejected");

    layer.samples.back().depth =
        -std::nextafter(tolerance, 0.0);
    static_cast<void>(PreparedComplexSspLayer(
        layer, AcousticInterpolation::CLinear,
        lossPerWavelengthContext(), 1));
}

void testNearEndpointQueryUsesTrueDepth()
{
    AcousticLayer layer = twoNodeFluidLayer();
    const PreparedComplexSspLayer prepared(
        layer, AcousticInterpolation::CLinear,
        lossPerWavelengthContext(), 1);
    const double query =
        std::nextafter(
            100.0 * std::numeric_limits<float>::epsilon(), 0.0);
    const InterpolatedMaterial actual = prepared.evaluate(query);
    const double expected = 1000.0 + 100.0 * query;
    require(actual.cp.real() > 1000.0 &&
                std::abs(actual.cp.real() - expected) <= 1.0e-12,
            "ALG-022 near-endpoint query must interpolate at true depth");
}

void testInvalidSampleMessagesIdentifyComponentSampleAndDepth()
{
    const auto construct =
        [](const AcousticLayer &layer, std::size_t mediumNumber)
    {
        static_cast<void>(PreparedComplexSspLayer(
            layer, AcousticInterpolation::CLinear,
            lossPerWavelengthContext(), mediumNumber));
    };

    AcousticLayer layer =
        twoNodeElastic(AcousticInterpolation::CLinear).layers.front();
    layer.samples[1].cp = 0.0;
    require(
        invalidArgumentMessage(
            [&] { construct(layer, 1); },
            "ALG-022 zero cp must throw") ==
            "medium 1 cp sample 2 depth 100: sound speed must be positive",
        "ALG-022 zero cp error context mismatch");

    layer =
        twoNodeElastic(AcousticInterpolation::CLinear).layers.front();
    layer.samples[1].rho = 0.0;
    require(
        invalidArgumentMessage(
            [&] { construct(layer, 3); },
            "ALG-022 zero rho must throw") ==
            "medium 3 rho sample 2 depth 100: density must be positive",
        "ALG-022 zero rho error context mismatch");

    layer =
        twoNodeElastic(AcousticInterpolation::CLinear).layers.front();
    layer.samples[1].alphaP = -0.25;
    require(
        invalidArgumentMessage(
            [&] { construct(layer, 4); },
            "ALG-022 negative P attenuation must throw") ==
            "medium 4 cp sample 2 depth 100: attenuation must be non-negative",
        "ALG-022 negative P attenuation error context mismatch");

    layer =
        twoNodeElastic(AcousticInterpolation::CLinear).layers.front();
    layer.samples[1].alphaS = -0.5;
    require(
        invalidArgumentMessage(
            [&] { construct(layer, 5); },
            "ALG-022 negative S attenuation must throw") ==
            "medium 5 cs sample 2 depth 100: attenuation must be non-negative",
        "ALG-022 negative S attenuation error context mismatch");
}

void testMixedShearSpeedMessageIdentifiesOffendingSample()
{
    AcousticLayer layer =
        twoNodeElastic(AcousticInterpolation::CLinear).layers.front();
    layer.samples[0].cs = 0.0;
    layer.samples[0].alphaS = 0.0;
    layer.samples[1].cs = 800.0;
    layer.samples[1].alphaS = 0.0;
    const std::string message = invalidArgumentMessage(
        [&]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::CLinear,
                lossPerWavelengthContext(), 1));
        },
        "ALG-022 mixed shear speeds must throw");
    requireContains(message, "medium 1",
                    "ALG-022 mixed cs error must identify medium");
    requireContains(message, "cs",
                    "ALG-022 mixed cs error must identify component");
    requireContains(message, "sample 2",
                    "ALG-022 mixed cs error must identify sample");
}

void testNonFiniteSplineCoefficientIdentifiesInterval()
{
    AcousticLayer layer = fourNodePchipLayer();
    layer.samples[0].rho = 1.0;
    layer.samples[1].rho = 1.0;
    layer.samples[2].rho = 1.0;
    layer.samples[3].rho = std::numeric_limits<double>::max();
    const std::string message = invalidArgumentMessage(
        [&]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::Pchip,
                lossPerWavelengthContext(), 2));
        },
        "ALG-022 non-finite rho coefficient must throw");
    require(
        message ==
            "medium 2 rho interval 3: non-finite spline coefficient",
        "ALG-022 non-finite spline coefficient error mismatch");
}

void testNonFiniteCpSplineCoefficientIdentifiesInterval()
{
    AcousticLayer layer = fourNodePchipLayer();
    for (AcousticSample &sample : layer.samples)
    {
        sample.cp = 1500.0;
        sample.alphaP = 0.0;
    }
    layer.samples[3].cp = std::numeric_limits<double>::max();
    const std::string message = invalidArgumentMessage(
        [&]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::Pchip,
                lossPerWavelengthContext(), 6));
        },
        "ALG-022 non-finite cp coefficient must throw");
    require(
        message ==
            "medium 6 cp interval 3: non-finite spline coefficient",
        "ALG-022 non-finite cp coefficient error mismatch: " + message);
}

void testNonFiniteCsSplineCoefficientIdentifiesInterval()
{
    AcousticLayer layer = fourNodePchipLayer();
    for (AcousticSample &sample : layer.samples)
    {
        sample.cs = 700.0;
        sample.alphaS = 0.0;
    }
    layer.samples[3].cs = std::numeric_limits<double>::max();
    const std::string message = invalidArgumentMessage(
        [&]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::Pchip,
                lossPerWavelengthContext(), 7));
        },
        "ALG-022 non-finite cs coefficient must throw");
    require(
        message ==
            "medium 7 cs interval 3: non-finite spline coefficient",
        "ALG-022 non-finite cs coefficient error mismatch: " + message);
}

void testDegenerateFiniteSplineSystemIdentifiesNodeRange()
{
    const double huge = std::numeric_limits<double>::max();
    AcousticLayer layer;
    layer.baseMesh = 4;
    layer.topDepth = 0.0;
    layer.bottomDepth = huge;
    layer.samples = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {0.25, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {0.5, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {huge, 1500.0, 0.0, 1.0, 0.0, 0.0},
    };
    const std::string message = invalidArgumentMessage(
        [&]
        {
            static_cast<void>(PreparedComplexSspLayer(
                layer, AcousticInterpolation::CubicSpline,
                lossPerWavelengthContext(), 2));
        },
        "ALG-022 degenerate finite spline system must throw");
    require(
        message ==
            "medium 2 cp nodes 1-4: singular spline system",
        "ALG-022 singular spline system error mismatch: " + message);
}

void testPreparedLayerRebuildsComplexNodesForEachFrequency()
{
    const AcousticLayer layer = twoNodeFluidLayer();
    AttenuationContext first;
    first.unit = 'N';
    first.frequency = 100.0;
    first.referenceFrequency = 100.0;
    AttenuationContext second = first;
    second.frequency = 200.0;

    const PreparedComplexSspLayer firstPrepared(
        layer, AcousticInterpolation::CLinear, first, 1);
    const PreparedComplexSspLayer secondPrepared(
        layer, AcousticInterpolation::CLinear, second, 1);
    const std::complex<double> firstExpected =
        complexSoundSpeed(
            layer.samples[0].depth,
            layer.samples[0].cp,
            layer.samples[0].alphaP, first);
    const std::complex<double> secondExpected =
        complexSoundSpeed(
            layer.samples[0].depth,
            layer.samples[0].cp,
            layer.samples[0].alphaP, second);

    require(nearComplex(firstPrepared.top().cp, firstExpected) &&
                nearComplex(secondPrepared.top().cp, secondExpected),
            "ALG-022 each prepared layer must complexify at its frequency");
    require(std::abs(
                secondPrepared.top().cp.imag() -
                0.5 * firstPrepared.top().cp.imag()) <= 1.0e-13,
            "ALG-022 doubled frequency must halve N-unit imaginary speed");
}

void testPreparedComplexCLinear()
{
    const AcousticLayer layer = threeNodeElasticLayer();
    const AttenuationContext context = lossPerWavelengthContext();
    const PreparedComplexSspLayer prepared(
        layer, AcousticInterpolation::CLinear, context, 1);

    const InterpolatedMaterial left = prepared.evaluate(7.5);
    require(nearComplex(left.cp, {1100.0, 20.0}),
            "ALG-022 C left-segment cp mismatch");
    require(nearComplex(left.cs, {550.0, 10.0}),
            "ALG-022 C left-segment cs mismatch");
    require(std::abs(left.rho - 1.25) <= 1.0e-13,
            "ALG-022 C left-segment density mismatch");

    const InterpolatedMaterial right = prepared.evaluate(12.5);
    require(nearComplex(right.cp, {1400.0, 50.0}),
            "ALG-022 C right-segment cp mismatch");
    require(nearComplex(right.cs, {700.0, 25.0}),
            "ALG-022 C right-segment cs mismatch");
    require(std::abs(right.rho - 1.75) <= 1.0e-13,
            "ALG-022 C right-segment density mismatch");

    const InterpolatedMaterial exactNode = prepared.evaluate(10.0);
    require(exactNode.cp == complexSoundSpeed(
                                10.0, 1200.0, 30.0 / 1200.0,
                                context),
            "ALG-022 exact C node must return stored complex cp");
    require(exactNode.cs == complexSoundSpeed(
                                10.0, 600.0, 15.0 / 600.0,
                                context),
            "ALG-022 exact C node must return stored complex cs");
    require(prepared.top().cp ==
                    complexSoundSpeed(
                        5.0, 1000.0, 10.0 / 1000.0, context) &&
                prepared.bottom().cp == complexSoundSpeed(
                                            15.0, 1600.0,
                                            70.0 / 1600.0,
                                            context),
            "ALG-022 C endpoint nodes must remain exact");
    require(prepared.elastic(), "ALG-022 elastic layer classification mismatch");
}

void testPreparedComplexN2LinearAndZeroShear()
{
    const AcousticLayer layer = twoNodeFluidLayer();
    const PreparedComplexSspLayer prepared(
        layer, AcousticInterpolation::N2Linear,
        lossPerWavelengthContext(), 2);

    const InterpolatedMaterial midpoint = prepared.evaluate(5.0);
    require(nearComplex(
                midpoint.cp,
                {1264.9110640673518, 12.649110640673518}),
            "ALG-022 N midpoint must interpolate complex inverse squared speed");
    require(std::abs(midpoint.rho - 1.5) <= 1.0e-13,
            "ALG-022 N midpoint density mismatch");
    require(prepared.top().cs == std::complex<double>{} &&
                midpoint.cs == std::complex<double>{} &&
                prepared.bottom().cs == std::complex<double>{},
            "ALG-022 all-zero cs must remain exact complex zero");
    const AttenuationContext context =
        lossPerWavelengthContext();
    require(prepared.evaluate(0.0).cp == complexSoundSpeed(
                                               0.0, 1000.0,
                                               10.0 / 1000.0,
                                               context) &&
                prepared.evaluate(10.0).cp == complexSoundSpeed(
                                                  10.0, 2000.0,
                                                  20.0 / 2000.0,
                                                  context),
            "ALG-022 exact N nodes must return stored complex cp");
    require(!prepared.elastic(), "ALG-022 fluid layer classification mismatch");
}

void testPreparedComplexCubicSplineSpecialCases()
{
    const AttenuationContext context = lossPerWavelengthContext();

    const PreparedComplexSspLayer twoNode(
        twoNodeSplineLayer(), AcousticInterpolation::CubicSpline,
        context, 1);
    const InterpolatedMaterial twoNodeMidpoint =
        twoNode.evaluate(0.5);
    require(nearComplex(twoNodeMidpoint.cp, {1550.0, 1.0}),
            "ALG-022 S two-node cp mismatch");
    require(nearComplex(twoNodeMidpoint.cs, {750.0, 2.0}),
            "ALG-022 S two-node cs mismatch");
    require(std::abs(twoNodeMidpoint.rho - 1.5) <= 1.0e-13,
            "ALG-022 S two-node density mismatch");

    const PreparedComplexSspLayer threeNode(
        threeNodeSplineLayer(), AcousticInterpolation::CubicSpline,
        context, 1);
    const InterpolatedMaterial threeNodeLeft =
        threeNode.evaluate(0.5);
    require(nearComplex(threeNodeLeft.cp, {1500.25, 0.25}),
            "ALG-022 S three-node left cp mismatch");
    require(nearComplex(threeNodeLeft.cs, {700.25, 0.5}),
            "ALG-022 S three-node left cs mismatch");
    require(std::abs(threeNodeLeft.rho - 1.25) <= 1.0e-13,
            "ALG-022 S three-node left density mismatch");
    const InterpolatedMaterial threeNodeRight =
        threeNode.evaluate(1.5);
    require(nearComplex(threeNodeRight.cp, {1502.25, 2.25}),
            "ALG-022 S three-node right cp mismatch");
    require(nearComplex(threeNodeRight.cs, {702.25, 4.5}),
            "ALG-022 S three-node right cs mismatch");
    require(std::abs(threeNodeRight.rho - 3.25) <= 1.0e-13,
            "ALG-022 S three-node right density mismatch");

    const PreparedComplexSspLayer fourNode(
        fourNodeSplineLayer(), AcousticInterpolation::CubicSpline,
        context, 1);
    const std::array<double, 3> queries = {0.5, 1.5, 2.5};
    const std::array<std::complex<double>, 3> expectedCp = {
        std::complex<double>{1500.125, 0.125},
        std::complex<double>{1503.375, 3.375},
        std::complex<double>{1515.625, 15.625},
    };
    const std::array<std::complex<double>, 3> expectedCs = {
        std::complex<double>{700.125, 0.25},
        std::complex<double>{703.375, 6.75},
        std::complex<double>{715.625, 31.25},
    };
    const std::array<double, 3> expectedRho = {
        1.125, 4.375, 16.625};
    for (std::size_t index = 0; index < queries.size(); ++index)
    {
        const InterpolatedMaterial actual =
            fourNode.evaluate(queries[index]);
        require(nearComplex(actual.cp, expectedCp[index]),
                "ALG-022 S four-node cp mismatch");
        require(nearComplex(actual.cs, expectedCs[index]),
                "ALG-022 S four-node cs mismatch");
        require(std::abs(actual.rho - expectedRho[index]) <=
                    1.0e-13 +
                        5.0e-13 * std::abs(expectedRho[index]),
                "ALG-022 S four-node density mismatch");
    }

    const PreparedComplexSspLayer nonUniform(
        nonUniformFourNodeSplineLayer(),
        AcousticInterpolation::CubicSpline, context, 1);
    const InterpolatedMaterial nonUniformValue =
        nonUniform.evaluate(3.0);
    require(nearComplex(nonUniformValue.cp, {1527.0, 27.0}),
            "ALG-022 S non-uniform cp mismatch");
    require(nearComplex(nonUniformValue.cs, {727.0, 54.0}),
            "ALG-022 S non-uniform cs mismatch");
    require(std::abs(nonUniformValue.rho - 28.0) <=
                1.0e-13 + 5.0e-13 * 28.0,
            "ALG-022 S non-uniform density mismatch");
}

void requireSplineContinuity(const AcousticLayer &layer,
                             const std::vector<double> &internalKnots)
{
    const AttenuationContext context = lossPerWavelengthContext();
    const PreparedComplexSspLayer prepared(
        layer, AcousticInterpolation::CubicSpline, context, 1);
    for (const double knot : internalKnots)
    {
        const auto sample = std::find_if(
            layer.samples.begin(), layer.samples.end(),
            [knot](const AcousticSample &candidate)
            {
                return candidate.depth == knot;
            });
        require(sample != layer.samples.end(),
                "ALG-022 S continuity fixture knot missing");
        const InterpolatedMaterial exact = prepared.evaluate(knot);
        require(exact.cp == complexSoundSpeed(
                                knot, sample->cp, sample->alphaP,
                                context) &&
                    exact.cs == complexSoundSpeed(
                                    knot, sample->cs, sample->alphaS,
                                    context) &&
                    exact.rho == sample->rho,
                "ALG-022 S exact knot must return stored node");

        const double leftDepth = std::nextafter(
            knot, -std::numeric_limits<double>::infinity());
        const double rightDepth = std::nextafter(
            knot, std::numeric_limits<double>::infinity());
        const InterpolatedMaterial left =
            prepared.evaluate(leftDepth);
        const InterpolatedMaterial right =
            prepared.evaluate(rightDepth);
        require(std::isfinite(left.cp.real()) &&
                    std::isfinite(left.cp.imag()) &&
                    std::isfinite(left.cs.real()) &&
                    std::isfinite(left.cs.imag()) &&
                    std::isfinite(left.rho) &&
                    std::isfinite(right.cp.real()) &&
                    std::isfinite(right.cp.imag()) &&
                    std::isfinite(right.cs.real()) &&
                    std::isfinite(right.cs.imag()) &&
                    std::isfinite(right.rho),
                "ALG-022 S knot neighbours must remain finite");
        require(nearComplex(left.cp, exact.cp) &&
                    nearComplex(right.cp, exact.cp) &&
                    nearComplex(left.cs, exact.cs) &&
                    nearComplex(right.cs, exact.cs) &&
                    std::abs(left.rho - exact.rho) <=
                        1.0e-13 +
                            5.0e-13 * std::abs(exact.rho) &&
                    std::abs(right.rho - exact.rho) <=
                        1.0e-13 +
                            5.0e-13 * std::abs(exact.rho),
                "ALG-022 S spline must be continuous at knots");
    }
}

void testPreparedComplexCubicSplineNodeContinuity()
{
    requireSplineContinuity(threeNodeSplineLayer(), {1.0});
    requireSplineContinuity(fourNodeSplineLayer(), {1.0, 2.0});
    requireSplineContinuity(
        nonUniformFourNodeSplineLayer(), {0.5, 2.0});
}

void testPreparedComplexPchipTwoNodeIsLinear()
{
    const AttenuationContext context = lossPerWavelengthContext();
    const PreparedComplexSspLayer prepared(
        twoNodeSplineLayer(), AcousticInterpolation::Pchip,
        context, 1);
    const InterpolatedMaterial midpoint = prepared.evaluate(0.5);
    require(nearComplex(midpoint.cp, {1550.0, 1.0}),
            "ALG-022 P two-node cp must be linear");
    require(nearComplex(midpoint.cs, {750.0, 2.0}),
            "ALG-022 P two-node cs must be linear");
    require(std::abs(midpoint.rho - 1.5) <= 1.0e-13,
            "ALG-022 P two-node density must be linear");

    const PreparedComplexSspLayer fluid(
        twoNodeFluidLayer(), AcousticInterpolation::Pchip,
        context, 2);
    require(fluid.top().cs == std::complex<double>{} &&
                fluid.evaluate(5.0).cs == std::complex<double>{} &&
                fluid.bottom().cs == std::complex<double>{},
            "ALG-022 P fluid cs must remain exact complex zero");
}

void testPreparedComplexPchipThreeNodeEndpointLimits()
{
    const AttenuationContext context =
        lossPerWavelengthContext();
    const PreparedComplexSspLayer prepared(
        threeNodePchipLayer(), AcousticInterpolation::Pchip,
        context, 1);

    const InterpolatedMaterial left = prepared.evaluate(0.5);
    require(nearComplex(left.cp, {1568.75, 0.0}),
            "ALG-022 P three-node left cp endpoint limit mismatch");
    require(nearComplex(left.cs, {768.75, 0.0}),
            "ALG-022 P three-node left cs endpoint limit mismatch");
    require(std::abs(left.rho - 1.6875) <= 1.0e-13,
            "ALG-022 P three-node left density endpoint limit mismatch");

    const InterpolatedMaterial right = prepared.evaluate(1.5);
    require(nearComplex(right.cp, {1600.0, 0.3125}),
            "ALG-022 P three-node right cp endpoint limit mismatch");
    require(nearComplex(right.cs, {800.0, 0.625}),
            "ALG-022 P three-node right cs endpoint limit mismatch");
    require(std::abs(right.rho - 2.0) <= 1.0e-13,
            "ALG-022 P three-node right density endpoint limit mismatch");

    require(prepared.evaluate(1.0).cp ==
                    complexSoundSpeed(
                        1.0, 1600.0, 0.0, context) &&
                prepared.evaluate(2.0).cp ==
                    complexSoundSpeed(
                        2.0, 1600.0, 1.0 / 1600.0,
                        context),
            "ALG-022 P exact nodes must return stored complex values");
}

void testMatrixUsesFourNodePchipProjection()
{
    const AcousticMatrix matrix = buildAcousticMatrix(
        matrixCase(
            fourNodePchipLayer(),
            AcousticInterpolation::Pchip),
        1);
    const std::array<std::size_t, 3> indices = {1, 3, 5};
    const std::array<std::complex<double>, 3> expectedCp = {
        std::complex<double>{1500.25, 0.75},
        std::complex<double>{1502.25, 0.5},
        std::complex<double>{1506.25, 0.25},
    };
    const std::array<std::complex<double>, 3> expectedCs = {
        std::complex<double>{700.25, 1.5},
        std::complex<double>{702.25, 1.0},
        std::complex<double>{706.25, 0.5},
    };
    const std::array<double, 3> expectedRho = {
        1.25, 3.25, 7.25};
    for (std::size_t index = 0; index < indices.size(); ++index)
    {
        const std::size_t matrixIndex = indices[index];
        require(nearComplex(matrix.cp[matrixIndex], expectedCp[index]),
                "ALG-022 P four-node matrix cp projection mismatch");
        require(nearComplex(matrix.cs[matrixIndex], expectedCs[index]),
                "ALG-022 P four-node matrix cs projection mismatch");
        require(std::abs(matrix.rho[matrixIndex] - expectedRho[index]) <=
                    1.0e-13,
                "ALG-022 P four-node matrix density projection mismatch");
    }
}

void testPreparedComplexPchipEndpointThreeTimesCap()
{
    const PreparedComplexSspLayer prepared(
        threeNodePchipEndpointCapLayer(),
        AcousticInterpolation::Pchip,
        lossPerWavelengthContext(), 1);

    const InterpolatedMaterial left = prepared.evaluate(0.5);
    require(nearComplex(left.cp, {1500.875, 10.875}),
            "ALG-022 P left endpoint three-times cap mismatch");

    const InterpolatedMaterial right = prepared.evaluate(1.5);
    require(nearComplex(right.cs, {709.125, 9.125}),
            "ALG-022 P right endpoint three-times cap mismatch");
}

void testPreparedComplexPchipNonUniformInteriorProjection()
{
    const PreparedComplexSspLayer prepared(
        nonUniformFourNodePchipProjectionLayer(),
        AcousticInterpolation::Pchip,
        lossPerWavelengthContext(), 1);
    const std::array<double, 3> queries = {0.5, 2.0, 3.5};
    const std::array<std::complex<double>, 3> expectedCp = {
        std::complex<double>{1500.625, 0.625},
        std::complex<double>{1501.25, 1.25},
        std::complex<double>{1536.75, 36.75},
    };
    const std::array<std::complex<double>, 3> expectedCs = {
        std::complex<double>{700.625, 1.25},
        std::complex<double>{701.25, 2.5},
        std::complex<double>{736.75, 73.5},
    };
    const std::array<double, 3> expectedRho = {
        1.625, 2.25, 37.75};
    for (std::size_t index = 0; index < queries.size(); ++index)
    {
        const InterpolatedMaterial actual =
            prepared.evaluate(queries[index]);
        require(nearComplex(actual.cp, expectedCp[index]),
                "ALG-022 P non-uniform cp projection mismatch");
        require(nearComplex(actual.cs, expectedCs[index]),
                "ALG-022 P non-uniform cs projection mismatch");
        require(std::abs(actual.rho - expectedRho[index]) <=
                    1.0e-13,
                "ALG-022 P non-uniform density projection mismatch");
    }
}

void testRealCpCAndNAndSAndPInterpolation()
{
    require(std::abs(interpolateRealCp(
                         threeNodeElasticLayer(),
                         AcousticInterpolation::CLinear, 7.5, 1) -
                     1100.0) <= 1.0e-13,
            "ALG-022 real C interpolation mismatch");
    require(std::abs(interpolateRealCp(
                         twoNodeFluidLayer(),
                         AcousticInterpolation::N2Linear, 5.0, 2) -
                     1264.9110640673518) <= 1.0e-10,
            "ALG-022 real N interpolation mismatch");
    require(std::abs(interpolateRealCp(
                         threeNodeSplineLayer(),
                         AcousticInterpolation::CubicSpline,
                         0.5, 3) -
                     1500.25) <= 1.0e-13,
            "ALG-022 real S interpolation mismatch");
    require(std::abs(interpolateRealCp(
                         nonUniformFourNodeSplineLayer(),
                         AcousticInterpolation::CubicSpline,
                         3.0, 3) -
                     1527.0) <= 1.0e-13,
            "ALG-022 real non-uniform S interpolation mismatch");
    require(std::abs(interpolateRealCp(
                         threeNodePchipLayer(),
                         AcousticInterpolation::Pchip,
                         0.5, 3) -
                     1568.75) <= 1.0e-13,
            "ALG-022 real P interpolation mismatch");
    requireInvalidArgument(
        []
        {
            static_cast<void>(interpolateRealCp(
                twoNodeFluidLayer(), AcousticInterpolation::CLinear,
                -1.0, 2));
        },
        "ALG-022 real cp interpolation must reject out-of-layer queries");
}

void testSoundSpeedAtUsesSourceAlignedInterpolation()
{
    AcousticCase spline;
    spline.interpolation = AcousticInterpolation::CubicSpline;
    spline.layers.push_back(threeNodeSplineLayer());
    require(std::abs(soundSpeedAt(spline, 0.5) - 1500.25) <= 1.0e-13,
            "ALG-022 soundSpeedAt S midpoint mismatch");

    AcousticCase pchip;
    pchip.interpolation = AcousticInterpolation::Pchip;
    pchip.layers.push_back(threeNodePchipLayer());
    require(std::abs(soundSpeedAt(pchip, 0.5) - 1568.75) <= 1.0e-13,
            "ALG-022 soundSpeedAt P midpoint mismatch");

    AcousticCase twoNode = twoNodeElastic(AcousticInterpolation::CLinear);
    require(std::abs(soundSpeedAt(twoNode, 50.0) - 1650.0) <= 1.0e-13,
            "ALG-022 soundSpeedAt C midpoint mismatch");

    twoNode.interpolation = AcousticInterpolation::N2Linear;
    require(std::abs(soundSpeedAt(twoNode, 50.0) -
                     1629.6434287653337) <= 1.0e-10,
            "ALG-022 soundSpeedAt N midpoint mismatch");

    AcousticCase interfaces;
    interfaces.interpolation = AcousticInterpolation::CLinear;
    interfaces.top.cp = 1400.0;
    interfaces.bottom.cp = 2200.0;
    AcousticLayer upper;
    upper.baseMesh = 2;
    upper.topDepth = 0.0;
    upper.bottomDepth = 10.0;
    upper.samples = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {10.0, 1510.0, 0.0, 1.0, 0.0, 0.0},
    };
    AcousticLayer lower;
    lower.baseMesh = 2;
    lower.topDepth = 10.0;
    lower.bottomDepth = 20.0;
    lower.samples = {
        {10.0, 1700.0, 700.0, 1.5, 0.0, 0.0},
        {20.0, 1800.0, 800.0, 1.6, 0.0, 0.0},
    };
    interfaces.layers = {upper, lower};
    require(soundSpeedAt(interfaces, 0.0) == 1400.0,
            "ALG-022 global top halfspace must have priority");
    require(soundSpeedAt(interfaces, 20.0) == 2200.0,
            "ALG-022 global bottom halfspace must have priority");
    require(soundSpeedAt(interfaces, 10.0) == 1510.0,
            "ALG-022 internal interface must use upper medium endpoint");
}

void testMatrixUsesPreparedThreeNodePaths()
{
    AcousticLayer layer = threeNodeElasticLayer();
    layer.baseMesh = 4;
    const AttenuationContext context =
        lossPerWavelengthContext();

    const AcousticMatrix c = buildAcousticMatrix(
        matrixCase(layer, AcousticInterpolation::CLinear), 1);
    require(c.depth.size() == 5 &&
                c.depth[0] == 5.0 &&
                c.depth[2] == 10.0 &&
                c.depth[4] == 15.0,
            "ALG-022 aligned C matrix grid mismatch");
    require(c.cp[0] == complexSoundSpeed(
                           5.0, 1000.0, 10.0 / 1000.0,
                           context) &&
                c.cp[2] == complexSoundSpeed(
                               10.0, 1200.0, 30.0 / 1200.0,
                               context) &&
                c.cp[4] == complexSoundSpeed(
                               15.0, 1600.0, 70.0 / 1600.0,
                               context),
            "ALG-022 C matrix top/exact-node/bottom wiring mismatch");
    require(nearComplex(c.cp[1], {1100.0, 20.0}) &&
                nearComplex(c.cp[3], {1400.0, 50.0}),
            "ALG-022 C matrix evaluate wiring mismatch");
    require(nearComplex(c.cs[1], {550.0, 10.0}) &&
                nearComplex(c.cs[3], {700.0, 25.0}) &&
                std::abs(c.rho[1] - 1.25) <= 1.0e-13 &&
                std::abs(c.rho[3] - 1.75) <= 1.0e-13,
            "ALG-022 C matrix material wiring mismatch");

    const AcousticMatrix n = buildAcousticMatrix(
        matrixCase(layer, AcousticInterpolation::N2Linear), 1);
    require(n.cp[0] == c.cp[0] &&
                n.cp[2] == c.cp[2] &&
                n.cp[4] == c.cp[4],
            "ALG-022 N matrix top/exact-node/bottom wiring mismatch");
    require(nearComplex(
                n.cp[1],
                {1086.5176073853672, 17.541404629783205}) &&
                nearComplex(
                    n.cp[3],
                    {1357.8097346599843, 43.098846066928019}),
            "ALG-022 N matrix cp evaluate wiring mismatch");
    require(nearComplex(
                n.cs[1],
                {543.25880369268361, 8.7707023148916026}) &&
                nearComplex(
                    n.cs[3],
                    {678.90486732999216, 21.549423033464009}) &&
                std::abs(n.rho[1] - 1.25) <= 1.0e-13 &&
                std::abs(n.rho[3] - 1.75) <= 1.0e-13,
            "ALG-022 N matrix material wiring mismatch");
}

void testMatrixFluidShearRemainsExactZero()
{
    const AcousticMatrix fluid = buildAcousticMatrix(
        matrixCase(
            twoNodeFluidLayer(),
            AcousticInterpolation::N2Linear),
        1);
    require(fluid.cs.size() == 3 &&
                fluid.cs[0] == std::complex<double>{} &&
                fluid.cs[1] == std::complex<double>{} &&
                fluid.cs[2] == std::complex<double>{},
            "ALG-022 matrix fluid cs must be exact zero at top/internal/bottom");
}

using OracleJson = nlohmann::json;

std::filesystem::path fixturePath(const std::string &name)
{
    return std::filesystem::path(OPENOCEANKRAKENC_SOURCE_DIR) /
           "for_test" / "fixtures" / name;
}

OracleJson readSspInterpolationOracle()
{
    const std::filesystem::path path =
        fixturePath("ssp_interp_fortran_oracle.json");
    std::ifstream stream(path);
    require(stream.is_open(),
            "ALG-022 Fortran SSP oracle is missing: " + path.string());
    OracleJson oracle;
    stream >> oracle;
    return oracle;
}

bool isSha256(const OracleJson &value)
{
    static const std::regex sha256("^[0-9a-f]{64}$");
    return value.is_string() &&
           std::regex_match(value.get<std::string>(), sha256);
}

void requireFiniteNumber(const OracleJson &object,
                         const std::string &key,
                         const std::string &location)
{
    require(object.contains(key) && object.at(key).is_number(),
            location + "." + key + " must be numeric");
    require(std::isfinite(object.at(key).get<double>()),
            location + "." + key + " must be finite");
}

void requireProvenanceFile(const OracleJson &record,
                           const std::string &location)
{
    require(record.is_object(), location + " must be an object");
    require(record.contains("path") && record.at("path").is_string() &&
                !record.at("path").get<std::string>().empty(),
            location + ".path must be a non-empty string");
    require(std::filesystem::path(
                record.at("path").get<std::string>()).is_absolute(),
            location + ".path must be absolute");
    require(record.contains("sha256") && isSha256(record.at("sha256")),
            location + ".sha256 must be 64 lowercase hexadecimal digits");
}

void requireCommandRecord(const OracleJson &record,
                          const std::string &location)
{
    require(record.is_object(), location + " must be an object");
    require(record.contains("command") && record.at("command").is_array() &&
                !record.at("command").empty(),
            location + ".command must be a non-empty array");
    for (const OracleJson &argument : record.at("command"))
    {
        require(argument.is_string() &&
                    !argument.get<std::string>().empty(),
                location + ".command arguments must be non-empty strings");
    }
    require(record.contains("returncode") &&
                record.at("returncode").is_number_integer() &&
                record.at("returncode").get<int>() == 0,
            location + ".returncode must be integer zero");
    require(record.contains("cwd") && record.at("cwd").is_string() &&
                std::filesystem::path(
                    record.at("cwd").get<std::string>()).is_absolute(),
            location + ".cwd must be absolute");
}

void validateOracleContract(const OracleJson &oracle)
{
    require(oracle.is_object(), "ALG-022 oracle root must be an object");
    require(oracle.at("schema") ==
                "OpenOcean-Krakenc.ssp-interpolation-oracle",
            "ALG-022 oracle schema mismatch");
    require(oracle.at("schema_version").is_number_integer() &&
                oracle.at("schema_version").get<int>() == 1,
            "ALG-022 oracle schema version mismatch");
    require(oracle.at("generated_date") == "2026-07-27",
            "ALG-022 oracle generated date mismatch");

    const OracleJson &provenance = oracle.at("provenance");
    requireProvenanceFile(provenance.at("krakenc"),
                          "provenance.krakenc");
    requireProvenanceFile(provenance.at("libMisc"),
                          "provenance.libMisc");
    for (const std::string groupName : {"sources", "inputs"})
    {
        const OracleJson &group = provenance.at(groupName);
        require(group.is_object() && !group.empty(),
                "provenance." + groupName + " must be non-empty");
        for (auto entry = group.begin(); entry != group.end(); ++entry)
        {
            requireProvenanceFile(
                entry.value(),
                "provenance." + groupName + "." + entry.key());
        }
    }
    const OracleJson &driver = provenance.at("driver");
    requireProvenanceFile(driver, "provenance.driver");
    require(driver.at("compile_command").is_array() &&
                !driver.at("compile_command").empty(),
            "provenance.driver.compile_command must be non-empty");
    for (const OracleJson &argument : driver.at("compile_command"))
    {
        require(argument.is_string() &&
                    !argument.get<std::string>().empty(),
                "driver compile arguments must be non-empty strings");
    }
    require(driver.at("compile_returncode").is_number_integer() &&
                driver.at("compile_returncode").get<int>() == 0,
            "provenance.driver.compile_returncode must be integer zero");
    require(driver.at("compile_cwd").is_string() &&
                std::filesystem::path(
                    driver.at("compile_cwd").get<std::string>()).is_absolute(),
            "provenance.driver.compile_cwd must be absolute");

    const OracleJson &commands = oracle.at("commands");
    for (const std::string profile : {"N", "C", "P", "S"})
    {
        requireCommandRecord(commands.at("profile_runs").at(profile),
                             "commands.profile_runs." + profile);
        requireCommandRecord(commands.at("krakenc_runs").at(profile),
                             "commands.krakenc_runs." + profile);
    }
    requireCommandRecord(
        commands.at("krakenc_runs").at("lossy_gradient_interpolation"),
        "commands.krakenc_runs.lossy_gradient_interpolation");
}

void requireFiniteWavenumbers(const OracleJson &profile,
                              const std::string &location)
{
    require(profile.at("mode_count").is_number_integer() &&
                profile.at("mode_count").get<int>() > 0,
            location + ".mode_count must be a positive integer");
    const OracleJson &wavenumbers = profile.at("wavenumbers");
    require(wavenumbers.is_array() && !wavenumbers.empty(),
            location + ".wavenumbers must be non-empty");
    require(static_cast<std::size_t>(
                profile.at("mode_count").get<int>()) == wavenumbers.size(),
            location + " mode count must equal wavenumber count");
    for (std::size_t index = 0; index < wavenumbers.size(); ++index)
    {
        const std::string item =
            location + ".wavenumbers[" + std::to_string(index) + "]";
        requireFiniteNumber(wavenumbers.at(index), "real", item);
        requireFiniteNumber(wavenumbers.at(index), "imag", item);
    }
}

constexpr double oracleComplexAbsoluteTolerance = 1.0e-11;
constexpr double oracleComplexRelativeTolerance = 5.0e-13;
constexpr double oracleDensityAbsoluteTolerance = 1.0e-13;
constexpr double oracleDensityRelativeTolerance = 5.0e-13;
constexpr double oracleRelativeErrorDenominatorFloor = 1.0e-30;

bool oracleComplexNear(std::complex<double> actual,
                       std::complex<double> reference)
{
    return std::abs(actual - reference) <=
           oracleComplexAbsoluteTolerance +
               oracleComplexRelativeTolerance * std::abs(reference);
}

bool oracleDensityNear(double actual, double reference)
{
    return std::abs(actual - reference) <=
           oracleDensityAbsoluteTolerance +
               oracleDensityRelativeTolerance * std::abs(reference);
}

double oracleRelativeError(std::complex<double> actual,
                           std::complex<double> reference)
{
    return std::abs(actual - reference) /
           std::max(std::abs(reference),
                    oracleRelativeErrorDenominatorFloor);
}

bool oracleWavenumberNear(std::complex<double> actual,
                          std::complex<double> reference,
                          double maximumAllowedRelativeError)
{
    return oracleRelativeError(actual, reference) <=
           maximumAllowedRelativeError;
}

std::pair<double, double> toleranceBoundaryValues(
    double reference,
    double tolerance)
{
    double accepted = reference + tolerance;
    while (std::abs(accepted - reference) > tolerance)
    {
        accepted = std::nextafter(accepted, reference);
    }
    double rejected =
        std::nextafter(accepted, std::numeric_limits<double>::infinity());
    while (std::abs(rejected - reference) <= tolerance)
    {
        rejected = std::nextafter(
            rejected, std::numeric_limits<double>::infinity());
    }
    return {accepted, rejected};
}

void testOracleAcceptanceToleranceContract()
{
    constexpr double requiredComplexAbsoluteTolerance = 1.0e-11;
    constexpr double requiredComplexRelativeTolerance = 5.0e-13;
    constexpr double requiredDensityAbsoluteTolerance = 1.0e-13;
    constexpr double requiredDensityRelativeTolerance = 5.0e-13;
    constexpr double requiredWavenumberRelativeTolerance = 1.0e-6;

    std::string failures;
    const std::array<std::complex<double>, 2> complexReferences = {
        std::complex<double>{0.0, 0.0},
        std::complex<double>{768.0, 1024.0},
    };
    for (std::size_t index = 0;
         index < complexReferences.size(); ++index)
    {
        const std::complex<double> reference = complexReferences[index];
        const double tolerance =
            requiredComplexAbsoluteTolerance +
            requiredComplexRelativeTolerance * std::abs(reference);
        const auto boundary =
            toleranceBoundaryValues(reference.real(), tolerance);
        const std::complex<double> accepted{
            boundary.first, reference.imag()};
        const std::complex<double> rejected{
            boundary.second, reference.imag()};
        require(std::abs(accepted - reference) <= tolerance &&
                    std::abs(rejected - reference) > tolerance,
                "ALG-022 complex oracle boundary fixture is not "
                "discriminating");
        if (!oracleComplexNear(accepted, reference))
        {
            failures += " complex boundary rejected at reference " +
                        std::to_string(index) + ";";
        }
        if (oracleComplexNear(rejected, reference))
        {
            failures +=
                " complex value just above boundary accepted at reference " +
                std::to_string(index) + ";";
        }
    }

    const std::array<double, 2> densityReferences = {0.0, 4.0};
    for (std::size_t index = 0;
         index < densityReferences.size(); ++index)
    {
        const double reference = densityReferences[index];
        const double tolerance =
            requiredDensityAbsoluteTolerance +
            requiredDensityRelativeTolerance * std::abs(reference);
        const auto boundary =
            toleranceBoundaryValues(reference, tolerance);
        require(std::abs(boundary.first - reference) <= tolerance &&
                    std::abs(boundary.second - reference) > tolerance,
                "ALG-022 density oracle boundary fixture is not "
                "discriminating");
        if (!oracleDensityNear(boundary.first, reference))
        {
            failures += " density boundary rejected at reference " +
                        std::to_string(index) + ";";
        }
        if (oracleDensityNear(boundary.second, reference))
        {
            failures +=
                " density value just above boundary accepted at reference " +
                std::to_string(index) + ";";
        }
    }

    const std::complex<double> wavenumberReference{1.0, 0.0};
    const auto wavenumberBoundary = toleranceBoundaryValues(
        wavenumberReference.real(),
        requiredWavenumberRelativeTolerance);
    require(oracleRelativeError(
                {wavenumberBoundary.first, 0.0},
                wavenumberReference) <=
                    requiredWavenumberRelativeTolerance &&
                std::abs(wavenumberBoundary.second -
                         wavenumberReference.real()) >
                    requiredWavenumberRelativeTolerance,
            "ALG-022 wavenumber oracle boundary fixture is not discriminating");
    if (!oracleWavenumberNear(
            {wavenumberBoundary.first, 0.0},
            wavenumberReference,
            requiredWavenumberRelativeTolerance))
    {
        failures += " wavenumber boundary rejected;";
    }
    if (oracleWavenumberNear(
            {wavenumberBoundary.second, 0.0},
            wavenumberReference,
            requiredWavenumberRelativeTolerance))
    {
        failures += " wavenumber value just above boundary accepted;";
    }

    const std::complex<double> tinyReference{1.0e-40, 0.0};
    const std::complex<double> tinyActual{2.0001e-36, 0.0};
    const double requiredTinyRelativeError =
        std::abs(tinyActual - tinyReference) / 1.0e-30;
    require(requiredTinyRelativeError >
                    requiredWavenumberRelativeTolerance &&
                requiredTinyRelativeError <
                    3.0 * requiredWavenumberRelativeTolerance,
            "ALG-022 tiny-reference fixture is not discriminating");
    if (oracleWavenumberNear(
            tinyActual, tinyReference,
            requiredWavenumberRelativeTolerance))
    {
        failures += " tiny-reference wavenumber error hidden by denominator;";
    }

    require(failures.empty(),
            "ALG-022 oracle tolerance contract violations:" + failures);
}

struct OracleErrorStatistics
{
    double maximumAbsolute = 0.0;
    double maximumRelative = 0.0;
};

void updateErrorStatistics(OracleErrorStatistics &statistics,
                           std::complex<double> actual,
                           std::complex<double> expected)
{
    const double absoluteError = std::abs(actual - expected);
    statistics.maximumAbsolute =
        std::max(statistics.maximumAbsolute, absoluteError);
    statistics.maximumRelative = std::max(
        statistics.maximumRelative,
        oracleRelativeError(actual, expected));
}

OracleErrorStatistics comparePreparedLayerWithOracle(
    const AcousticCase &input,
    const OracleJson &profile,
    const std::string &profileType)
{
    require(input.layers.size() == 2,
            "ALG-022 " + profileType + " ENV must contain two layers");
    const AcousticLayer &layer = input.layers.at(1);
    const PreparedComplexSspLayer prepared(
        layer, input.interpolation,
        attenuationContext(input, layer.attenuationPower,
                           layer.transitionFrequency),
        2);
    const OracleJson &samples = profile.at("samples");
    require(samples.is_array() && samples.size() == 5,
            "profiles." + profileType +
                ".samples must contain exactly five samples");
    OracleErrorStatistics statistics;
    for (std::size_t index = 0; index < samples.size(); ++index)
    {
        const OracleJson &sample = samples.at(index);
        const std::string location =
            "profiles." + profileType + ".samples[" +
            std::to_string(index) + "]";
        for (const std::string key :
             {"depth", "cp_re", "cp_im", "cs_re", "cs_im", "rho"})
        {
            requireFiniteNumber(sample, key, location);
        }
        const double depth = sample.at("depth").get<double>();
        require(std::abs(depth - (50.0 + 12.5 * index)) <= 1.0e-12,
                location + ".depth is not on the controlled five-point grid");
        const InterpolatedMaterial actual = prepared.evaluate(depth);
        const std::complex<double> expectedCp{
            sample.at("cp_re").get<double>(),
            sample.at("cp_im").get<double>()};
        const std::complex<double> expectedCs{
            sample.at("cs_re").get<double>(),
            sample.at("cs_im").get<double>()};
        const double expectedRho = sample.at("rho").get<double>();
        require(std::isfinite(actual.cp.real()) &&
                    std::isfinite(actual.cp.imag()) &&
                    std::isfinite(actual.cs.real()) &&
                    std::isfinite(actual.cs.imag()) &&
                    std::isfinite(actual.rho),
                "ALG-022 " + profileType +
                    " prepared interpolation produced a non-finite value");
        require(oracleComplexNear(actual.cp, expectedCp),
                "ALG-022 " + profileType + " prepared cp mismatch at " +
                    std::to_string(depth) + " m");
        require(oracleComplexNear(actual.cs, expectedCs),
                "ALG-022 " + profileType + " prepared cs mismatch at " +
                    std::to_string(depth) + " m");
        require(oracleDensityNear(actual.rho, expectedRho),
                "ALG-022 " + profileType + " prepared rho mismatch at " +
                    std::to_string(depth) + " m");
        updateErrorStatistics(statistics, actual.cp, expectedCp);
        updateErrorStatistics(statistics, actual.cs, expectedCs);
        updateErrorStatistics(
            statistics, {actual.rho, 0.0}, {expectedRho, 0.0});
    }
    return statistics;
}

double compareSolvedWavenumbers(const AcousticCase &input,
                                const OracleJson &profile,
                                const std::string &profileType,
                                double maximumAllowedRelativeError)
{
    const AcousticSolveResult result = solveAcousticModes(input);
    const OracleJson &expected = profile.at("wavenumbers");
    require(result.modes.size() == expected.size(),
            "ALG-022 " + profileType +
                " mode count differs from the Fortran oracle");
    double maximumRelativeError = 0.0;
    for (std::size_t index = 0; index < result.modes.size(); ++index)
    {
        const std::complex<double> actual =
            result.modes.at(index).wavenumber;
        const std::complex<double> reference{
            expected.at(index).at("real").get<double>(),
            expected.at(index).at("imag").get<double>()};
        require(std::isfinite(actual.real()) &&
                    std::isfinite(actual.imag()),
                "ALG-022 " + profileType +
                    " solver produced a non-finite wavenumber");
        const double relativeError =
            oracleRelativeError(actual, reference);
        maximumRelativeError =
            std::max(maximumRelativeError, relativeError);
        require(oracleWavenumberNear(
                    actual, reference, maximumAllowedRelativeError),
                "ALG-022 " + profileType +
                    " ordered wavenumber relative error exceeds " +
                    std::to_string(maximumAllowedRelativeError) +
                    " at mode " +
                    std::to_string(index + 1));
    }
    return maximumRelativeError;
}

void testFrozenFortranSspOracle()
{
    const OracleJson oracle = readSspInterpolationOracle();
    validateOracleContract(oracle);
    OracleErrorStatistics interpolationErrors;
    double maximumProfileWavenumberError = 0.0;
    for (const std::pair<std::string, std::string> fixture : {
             std::pair<std::string, std::string>{"N", "ssp_interp_n.env"},
             {"C", "ssp_interp_c.env"},
             {"P", "ssp_interp_p.env"},
             {"S", "ssp_interp_s.env"}})
    {
        const OracleJson &profile =
            oracle.at("profiles").at(fixture.first);
        requireFiniteWavenumbers(
            profile, "profiles." + fixture.first);
        const AcousticCase input =
            readAcousticEnv(fixturePath(fixture.second));
        const OracleErrorStatistics profileInterpolationErrors =
            comparePreparedLayerWithOracle(
            input, profile, fixture.first);
        interpolationErrors.maximumAbsolute = std::max(
            interpolationErrors.maximumAbsolute,
            profileInterpolationErrors.maximumAbsolute);
        interpolationErrors.maximumRelative = std::max(
            interpolationErrors.maximumRelative,
            profileInterpolationErrors.maximumRelative);
        maximumProfileWavenumberError = std::max(
            maximumProfileWavenumberError,
            compareSolvedWavenumbers(
                input, profile, fixture.first, 1.0e-6));
    }

    const OracleJson &lossy =
        oracle.at("lossy_gradient_interpolation");
    requireFiniteWavenumbers(
        lossy, "lossy_gradient_interpolation");
    const double lossyMaximumRelativeError =
        compareSolvedWavenumbers(
            readAcousticEnv(
                fixturePath("lossy_gradient_interpolation.env")),
            lossy, "lossy_gradient_interpolation", 2.0e-5);
    require(lossyMaximumRelativeError <= 2.0e-5,
            "ALG-022 lossy maximum relative error exceeds 2e-5");
    std::cout << "ALG-022 Fortran oracle maximum prepared interpolation "
                 "absolute error: "
              << interpolationErrors.maximumAbsolute << '\n'
              << "ALG-022 Fortran oracle maximum prepared interpolation "
                 "relative error: "
              << interpolationErrors.maximumRelative << '\n'
              << "ALG-022 Fortran oracle maximum profile wavenumber "
                 "relative error: "
              << maximumProfileWavenumberError << '\n'
              << "ALG-022 Fortran oracle lossy maximum relative error: "
              << lossyMaximumRelativeError << '\n';
}

void testAlg022()
{
    testOracleAcceptanceToleranceContract();
    testFrozenFortranSspOracle();
    testMatrixRetainsDuplicateInterfaceNodes();
    testPreparedLayersUseOnlyLocalPolynomialNodes();
    testSampleEndpointToleranceIsStrict();
    testNearEndpointQueryUsesTrueDepth();
    testInvalidSampleMessagesIdentifyComponentSampleAndDepth();
    testMixedShearSpeedMessageIdentifiesOffendingSample();
    testNonFiniteSplineCoefficientIdentifiesInterval();
    testNonFiniteCpSplineCoefficientIdentifiesInterval();
    testNonFiniteCsSplineCoefficientIdentifiesInterval();
    testDegenerateFiniteSplineSystemIdentifiesNodeRange();
    testPreparedLayerRebuildsComplexNodesForEachFrequency();
    testPreparedComplexCLinear();
    testPreparedComplexN2LinearAndZeroShear();
    testPreparedComplexCubicSplineSpecialCases();
    testPreparedComplexCubicSplineNodeContinuity();
    testPreparedComplexPchipTwoNodeIsLinear();
    testPreparedComplexPchipThreeNodeEndpointLimits();
    testMatrixUsesFourNodePchipProjection();
    testPreparedComplexPchipEndpointThreeTimesCap();
    testPreparedComplexPchipNonUniformInteriorProjection();
    testRealCpCAndNAndSAndPInterpolation();
    testSoundSpeedAtUsesSourceAlignedInterpolation();
    testMatrixUsesPreparedThreeNodePaths();
    testMatrixFluidShearRemainsExactZero();

    const AcousticMatrix c = buildAcousticMatrix(
        twoNodeElastic(AcousticInterpolation::CLinear), 1);
    require(nearComplex(c.cp[1], {1650.0, 15.941349345488913}),
            "ALG-022 C must complexify nodes before interpolation");
    require(nearComplex(c.cs[1], {800.0, 4.764081413594388}),
            "ALG-022 C shear order mismatch");
    require(std::abs(c.rho[1] - 1.5) <= 1.0e-13,
            "ALG-022 C density mismatch");

    const AcousticMatrix n = buildAcousticMatrix(
        twoNodeElastic(AcousticInterpolation::N2Linear), 1);
    require(nearComplex(n.cp[1],
                        {1629.7148845925744, 13.314114275026862}),
            "ALG-022 N must interpolate complex inverse squared speed");
    require(nearComplex(n.cs[1],
                        {781.4336150195090, 3.5904710085476124}),
            "ALG-022 N shear order mismatch");

    AcousticCase bio;
    bio.frequency = 2000.0;
    bio.referenceFrequency = 2000.0;
    bio.attenuationUnit = 'W';
    bio.absorptionModel = OceanAbsorptionModel::Biological;
    bio.volumeAbsorption.biologicalLayers.push_back(
        {0.0, 10.0, 1000.0, 5.0, 2.0});
    bio.interpolation = AcousticInterpolation::CLinear;
    AcousticLayer bioLayer;
    bioLayer.baseMesh = 2;
    bioLayer.topDepth = 0.0;
    bioLayer.bottomDepth = 100.0;
    bioLayer.samples = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
    };
    bio.layers.push_back(bioLayer);
    const AcousticMatrix biological = buildAcousticMatrix(bio, 1);
    require(nearComplex(
                biological.cp[1], {1500.0, 0.0342137988528164}),
            "ALG-022 biological loss must be applied at raw nodes");
}
}

int main(int argc, char **argv)
{
    if (argc != 2 || std::string(argv[1]) != "ALG-022")
    {
        std::cerr << "exactly one remediation scenario ALG-022 is required\n";
        return 2;
    }
    try
    {
        testAlg022();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
