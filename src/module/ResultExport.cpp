#include "module/ResultExport.h"

#include "algorithm/ComplexNumerics.h"
#include "algorithm/ModeFileWriter.h"
#include "module/ParameterAdapters.h"

#include <complex>
#include <stdexcept>
#include <vector>

namespace OpenOceanKrakenc
{
namespace
{
char attenuationCode(AttenuationUnit unit)
{
    switch (unit)
    {
    case AttenuationUnit::MODE_F_dB_per_m_kHz: return 'F';
    case AttenuationUnit::MODE_L_params_lose: return 'L';
    case AttenuationUnit::MODE_M_dB_per_m: return 'M';
    case AttenuationUnit::MODE_m_dB_per_m: return 'm';
    case AttenuationUnit::MODE_N_Nepers_per_m: return 'N';
    case AttenuationUnit::MODE_Q_Quality_Factor: return 'Q';
    case AttenuationUnit::MODE_W_db_per_lambda: return 'W';
    }
    return 'W';
}

char boundaryCode(AcousticBoundaryType type)
{
    switch (type)
    {
    case AcousticBoundaryType::Vacuum: return 'V';
    case AcousticBoundaryType::Rigid: return 'R';
    case AcousticBoundaryType::HalfSpace: return 'A';
    case AcousticBoundaryType::ReflectionCoefficient: return 'F';
    case AcousticBoundaryType::InternalReflection: return 'P';
    }
    return 'V';
}

ModeBoundaryData makeBoundary(const AcousticBoundary &source,
                              const AcousticCase &input,
                              char attenuationUnit)
{
    AttenuationContext context = attenuationContext(
        input, source.attenuationPower, source.transitionFrequency);
    context.unit = attenuationUnit;
    ModeBoundaryData result;
    result.type = boundaryCode(source.type);
    result.cp = source.cp > 0.0
                    ? complexSoundSpeed(source.depth, source.cp,
                                        source.alphaP, context)
                    : std::complex<double>{};
    result.cs = source.cs > 0.0
                    ? complexSoundSpeed(source.depth, source.cs,
                                        source.alphaS, context)
                    : std::complex<double>{};
    result.rho = source.rho;
    result.depth = source.depth;
    return result;
}

ModeProfileData makeProfile(const AcousticCase &input,
                            const EigenParams &eigen,
                            char attenuationUnit)
{
    if (eigen.M < 1 || eigen.k.size() != eigen.M ||
        eigen.ModeZ.size() < 2 || eigen.PhiMode.rows() != eigen.M ||
        eigen.PhiMode.cols() != eigen.ModeZ.size())
    {
        throw std::invalid_argument("eigen output is incomplete for MOD export");
    }
    ModeProfileData profile;
    profile.depths.assign(eigen.ModeZ.data(), eigen.ModeZ.data() + eigen.ModeZ.size());
    profile.wavenumbers.assign(eigen.k.data(), eigen.k.data() + eigen.k.size());
    profile.modes.resize(static_cast<std::size_t>(eigen.M));
    for (Eigen::Index mode = 0; mode < eigen.PhiMode.rows(); ++mode)
    {
        auto &values = profile.modes[static_cast<std::size_t>(mode)];
        values.resize(static_cast<std::size_t>(eigen.PhiMode.cols()));
        for (Eigen::Index depth = 0; depth < eigen.PhiMode.cols(); ++depth)
        {
            values[static_cast<std::size_t>(depth)] = eigen.PhiMode(mode, depth);
        }
    }
    std::size_t firstAcoustic = input.layers.size();
    std::size_t lastAcoustic = input.layers.size();
    for (std::size_t medium = 0; medium < input.layers.size(); ++medium)
    {
        if (!input.layers[medium].samples.empty() &&
            input.layers[medium].samples.front().cs == 0.0)
        {
            if (firstAcoustic == input.layers.size()) firstAcoustic = medium;
            lastAcoustic = medium;
        }
    }
    if (firstAcoustic == input.layers.size())
    {
        throw std::invalid_argument("MOD export requires an acoustic layer");
    }
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        const AcousticLayer &layer = input.layers[medium];
        profile.meshCounts.push_back(layer.baseMesh);
        profile.materials.push_back("ACOUSTIC");
        profile.mediumDepths.push_back(layer.topDepth);
        profile.mediumDensities.push_back(
            layer.samples.empty() ? 0.0 : layer.samples.front().rho);
    }
    profile.top = makeBoundary(input.top, input, attenuationUnit);
    profile.bottom = makeBoundary(input.bottom, input, attenuationUnit);
    profile.top.depth = input.layers.front().topDepth;
    profile.bottom.depth = input.layers.back().bottomDepth;
    return profile;
}

std::vector<std::complex<double>> widen(
    const std::vector<std::complex<float>> &source)
{
    std::vector<std::complex<double>> result;
    result.reserve(source.size());
    for (std::complex<float> value : source)
    {
        result.emplace_back(value.real(), value.imag());
    }
    return result;
}

std::vector<double> realVector(const Eigen::VectorXd &source)
{
    return std::vector<double>(source.data(), source.data() + source.size());
}
}

std::filesystem::path resultPath(const std::string &root, const char *extension)
{
    if (root.empty()) throw std::invalid_argument("result root must not be empty");
    std::filesystem::path path = root;
    if (path.extension() != extension) path += extension;
    return path;
}

void exportModeResult(const OOKC_parameters &params,
                      const OOKC_output &output,
                      const std::filesystem::path &path)
{
    const std::vector<AcousticCase> inputs = toAcousticCases(params);
    if (output.eigen.size() != inputs.size())
    {
        throw std::logic_error("eigen results are unavailable for every profile");
    }
    ModeFileData modes;
    modes.title = params.Title;
    modes.frequency = params.freqinfo.freqvec.size() == 1
                          ? params.freqinfo.freqvec[0]
                          : params.freqinfo.freq;
    const char unit = attenuationCode(params.AttenUnit.attnUnit);
    static_cast<ModeProfileData &>(modes) = makeProfile(inputs.front(), output.eigen.front(), unit);
    for (std::size_t profile = 1; profile < inputs.size(); ++profile)
    {
        modes.additionalProfiles.push_back(makeProfile(inputs[profile], output.eigen[profile], unit));
    }
    writeModeFile(modes, path);
}

void exportShadeResult(const OOKC_parameters &params,
                       const OOKC_output &output,
                       const std::filesystem::path &path,
                       ShadeDataType type)
{
    const std::size_t expected = output.sourceCount *
                                 output.receiverDepthCount * output.rangeCount;
    if (expected == 0 || output.sourceCount != static_cast<std::size_t>(params.Pos.Sz.size()) ||
        output.receiverDepthCount != static_cast<std::size_t>(params.Pos.Rz.size()) ||
        output.rangeCount != static_cast<std::size_t>(params.Pos.Rr.size()))
    {
        throw std::logic_error("field result shape is unavailable for SHD export");
    }
    PressureField field;
    field.title = params.Title;
    field.frequency = params.freqinfo.freqvec.size() == 1
                          ? params.freqinfo.freqvec[0]
                          : params.freqinfo.freq;
    field.sourceDepths = realVector(params.Pos.Sz);
    field.receiverDepths = realVector(params.Pos.Rz);
    field.rangesMetres = realVector(params.Pos.Rr);
    field.values = widen(output.pressure);
    field.verticalValues = widen(output.verticalVelocity);
    field.horizontalValues = widen(output.horizontalVelocity);
    const std::vector<std::complex<double>> *selected = nullptr;
    switch (type)
    {
    case ShadeDataType::Pressure: selected = &field.values; break;
    case ShadeDataType::VerticalVelocity: selected = &field.verticalValues; break;
    case ShadeDataType::HorizontalVelocity: selected = &field.horizontalValues; break;
    }
    if (!selected || selected->size() != expected)
    {
        throw std::logic_error("selected field vector is unavailable for SHD export");
    }
    writeShadeFile(field, path, type, params.Pos.GridType);
}
}
