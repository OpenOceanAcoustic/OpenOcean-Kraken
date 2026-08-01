#include "OpenOceanKrakenKernelInterface.h"
#include "ThreadPool.h"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstring>
#include <stdexcept>

namespace py = pybind11;

namespace
{
    template <typename T>
    py::array_t<T> vector_array(const std::vector<T> &values)
    {
        py::array_t<T> result(values.size());
        if (!values.empty())
        {
            std::memcpy(result.mutable_data(), values.data(), values.size() * sizeof(T));
        }
        return result;
    }

    py::array_t<std::complex<float>> field_array(const OpenOceanKraken::FieldSnapshot &snapshot)
    {
        py::array_t<std::complex<float>> result({
            static_cast<py::ssize_t>(snapshot.source_count),
            static_cast<py::ssize_t>(snapshot.range_count),
            static_cast<py::ssize_t>(snapshot.depth_count)});
        if (static_cast<std::size_t>(result.size()) != snapshot.values.size())
        {
            throw std::runtime_error("Field snapshot dimensions do not match its value count.");
        }
        if (!snapshot.values.empty())
        {
            std::memcpy(
                result.mutable_data(),
                snapshot.values.data(),
                snapshot.values.size() * sizeof(std::complex<float>));
        }
        return result;
    }

    py::array_t<std::complex<double>> matrix_array(const Eigen::MatrixXcd &values)
    {
        py::array_t<std::complex<double>> result({
            static_cast<py::ssize_t>(values.rows()),
            static_cast<py::ssize_t>(values.cols())});
        auto output = result.mutable_unchecked<2>();
        for (Eigen::Index row = 0; row < values.rows(); ++row)
        {
            for (Eigen::Index column = 0; column < values.cols(); ++column)
            {
                output(row, column) = values(row, column);
            }
        }
        return result;
    }
}

PYBIND11_MODULE(OpenOceanKraken, module)
{
    using namespace OpenOceanKraken;
    module.doc() = "Python bindings for OpenOcean-Kraken";

    py::class_<ThreadPool>(module, "ThreadPool")
        .def(py::init<std::size_t>(), py::arg("threads"));

    py::enum_<SSP_Mode>(module, "SSP_Mode")
        .value("MODE_N_n2Linear", SSP_Mode::MODE_N_n2Linear)
        .value("MODE_C_cLinear", SSP_Mode::MODE_C_cLinear)
        .value("MODE_P_cPCHIP", SSP_Mode::MODE_P_cPCHIP)
        .value("MODE_S_cCubic", SSP_Mode::MODE_S_cCubic)
        .value("MODE_A_Analytic", SSP_Mode::MODE_A_Analytic)
        .export_values();
    py::enum_<Media_Mode>(module, "Media_Mode")
        .value("MODE_A_Acoustic", Media_Mode::MODE_A_Acoustic)
        .value("MODE_E_Elastic", Media_Mode::MODE_E_Elastic)
        .export_values();
    py::enum_<AttenuationUnit>(module, "AttenuationUnit")
        .value("MODE_F_dB_per_m_kHz", AttenuationUnit::MODE_F_dB_per_m_kHz)
        .value("MODE_L_params_lose", AttenuationUnit::MODE_L_params_lose)
        .value("MODE_M_dB_per_m", AttenuationUnit::MODE_M_dB_per_m)
        .value("MODE_m_dB_per_m", AttenuationUnit::MODE_m_dB_per_m)
        .value("MODE_N_Nepers_per_m", AttenuationUnit::MODE_N_Nepers_per_m)
        .value("MODE_Q_Quality_Factor", AttenuationUnit::MODE_Q_Quality_Factor)
        .value("MODE_W_db_per_lambda", AttenuationUnit::MODE_W_db_per_lambda)
        .export_values();
    py::enum_<OceanAbsorptionModel>(module, "OceanAbsorptionModel")
        .value("None_", OceanAbsorptionModel::None)
        .value("Thorpe", OceanAbsorptionModel::Thorpe)
        .value("FrancGarr", OceanAbsorptionModel::FrancGarr)
        .value("Biological", OceanAbsorptionModel::Biological)
        .export_values();
    py::enum_<BC_Mode>(module, "BC_Mode")
        .value("MODE_R_Rigid", BC_Mode::MODE_R_Rigid)
        .value("MODE_V_Vacuum", BC_Mode::MODE_V_Vacuum)
        .value("MODE_F_File", BC_Mode::MODE_F_File)
        .value("MODE_A_Half_space", BC_Mode::MODE_A_Half_space)
        .value("MODE_G_Grain", BC_Mode::MODE_G_Grain)
        .value("MODE_P_Precomputed", BC_Mode::MODE_P_Precomputed)
        .export_values();
    py::enum_<Source_Mode>(module, "Source_Mode")
        .value("MODE_R_Point", Source_Mode::MODE_R_Point)
        .value("MODE_X_Line", Source_Mode::MODE_X_Line)
        .export_values();
    py::enum_<Grid_Mode>(module, "Grid_Mode")
        .value("MODE_R_Rectangular", Grid_Mode::MODE_R_Rectangular)
        .value("MODE_I_Irregular", Grid_Mode::MODE_I_Irregular)
        .export_values();
    py::enum_<Run_Mode>(module, "Run_Mode")
        .value("MODE_M_Modes", Run_Mode::MODE_M_Modes)
        .value("MODE_F_Field", Run_Mode::MODE_F_Field)
        .value("MODE_B_Both", Run_Mode::MODE_B_Both)
        .export_values();
    py::enum_<CoherenceType>(module, "CoherenceType")
        .value("Coherent", CoherenceType::Coherent)
        .value("Incoherent", CoherenceType::Incoherent)
        .export_values();
    py::enum_<ModeType>(module, "ModeType")
        .value("Adiabatic", ModeType::Adiabatic)
        .value("Couple", ModeType::Couple)
        .export_values();

    py::class_<BiologicalAttenuationLayer>(
        module, "BiologicalAttenuationLayer")
        .def(py::init<>())
        .def_readwrite("Z1", &BiologicalAttenuationLayer::Z1)
        .def_readwrite("Z2", &BiologicalAttenuationLayer::Z2)
        .def_readwrite("f0", &BiologicalAttenuationLayer::f0)
        .def_readwrite("Q", &BiologicalAttenuationLayer::Q)
        .def_readwrite("a0", &BiologicalAttenuationLayer::a0);
    py::class_<Atten_Mode>(module, "Atten_Mode")
        .def(py::init<>())
        .def_readwrite("attnUnit", &Atten_Mode::attnUnit)
        .def_readwrite("absModel", &Atten_Mode::absModel)
        .def_readwrite(
            "biologicalLayers",
            &Atten_Mode::biologicalLayers);
    py::class_<HSInfo>(module, "HSInfo")
        .def(py::init<>())
        .def_readwrite("alphaR", &HSInfo::alphaR)
        .def_readwrite("alphaI", &HSInfo::alphaI)
        .def_readwrite("betaR", &HSInfo::betaR)
        .def_readwrite("betaI", &HSInfo::betaI)
        .def_readwrite("beta", &HSInfo::beta)
        .def_readwrite("ft", &HSInfo::ft)
        .def_readwrite("cp", &HSInfo::cp)
        .def_readwrite("cs", &HSInfo::cs)
        .def_readwrite("rho", &HSInfo::rho)
        .def_readwrite("Depth", &HSInfo::Depth)
        .def_readwrite("BC", &HSInfo::BC);
    py::class_<ssp::SSPLayer>(module, "SSPLayer")
        .def(py::init<>())
        .def_readwrite("npts", &ssp::SSPLayer::npts)
        .def_readwrite("nmesh", &ssp::SSPLayer::nmesh)
        .def_readwrite("beta", &ssp::SSPLayer::beta)
        .def_readwrite("ft", &ssp::SSPLayer::ft)
        .def_readwrite("sigma", &ssp::SSPLayer::sigma)
        .def_readwrite("Material", &ssp::SSPLayer::Material)
        .def_readwrite("z", &ssp::SSPLayer::z)
        .def_readwrite("rho", &ssp::SSPLayer::rho)
        .def_readwrite("alphaR", &ssp::SSPLayer::alphaR)
        .def_readwrite("alphaI", &ssp::SSPLayer::alphaI)
        .def_readwrite("betaR", &ssp::SSPLayer::betaR)
        .def_readwrite("betaI", &ssp::SSPLayer::betaI);
    py::class_<ssp::Range_Independent_Area>(module, "Range_Independent_Area")
        .def(py::init<>())
        .def_readwrite("SSPType", &ssp::Range_Independent_Area::SSPType)
        .def_readwrite("layers", &ssp::Range_Independent_Area::layers)
        .def_readwrite("HSTop", &ssp::Range_Independent_Area::HSTop)
        .def_readwrite("HSBot", &ssp::Range_Independent_Area::HSBot)
        .def_readwrite("Range", &ssp::Range_Independent_Area::Range)
        .def("addLayer", &ssp::Range_Independent_Area::addLayer)
        .def("clearLayer", &ssp::Range_Independent_Area::clearLayer)
        .def("NMedia", &ssp::Range_Independent_Area::NMedia)
        .def("set_Bottom_Line", &ssp::Range_Independent_Area::set_Bottom_Line)
        .def("set_Top_Line", &ssp::Range_Independent_Area::set_Top_Line)
        .def("set_Bottom_type", &ssp::Range_Independent_Area::set_Bottom_type)
        .def("set_Top_type", &ssp::Range_Independent_Area::set_Top_type);
    py::class_<ReflectionCoef>(module, "ReflectionCoef")
        .def(py::init<>())
        .def(py::init([](double theta, double magnitude, double phase) {
            return ReflectionCoef{theta, magnitude, phase};
        }))
        .def_readwrite("theta", &ReflectionCoef::theta)
        .def_readwrite("R", &ReflectionCoef::R)
        .def_readwrite("phi", &ReflectionCoef::phi);

    py::class_<FieldSnapshot>(module, "FieldSnapshot")
        .def_property_readonly("title", [](const FieldSnapshot &value) { return value.title; })
        .def_property_readonly("frequency", [](const FieldSnapshot &value) { return value.frequency; })
        .def_property_readonly("grid_type", [](const FieldSnapshot &value) { return value.grid_type; })
        .def_property_readonly("source_count", [](const FieldSnapshot &value) { return value.source_count; })
        .def_property_readonly("range_count", [](const FieldSnapshot &value) { return value.range_count; })
        .def_property_readonly("depth_count", [](const FieldSnapshot &value) { return value.depth_count; })
        .def_property_readonly("source_depths", [](const FieldSnapshot &value) { return vector_array(value.source_depths); })
        .def_property_readonly("receiver_ranges", [](const FieldSnapshot &value) { return vector_array(value.receiver_ranges); })
        .def_property_readonly("receiver_depths", [](const FieldSnapshot &value) { return vector_array(value.receiver_depths); })
        .def_property_readonly("values", &field_array);
    py::class_<ModeProfileSnapshot>(module, "ModeProfileSnapshot")
        .def_property_readonly("profile_range", [](const ModeProfileSnapshot &value) { return value.profile_range; })
        .def_property_readonly("depth", [](const ModeProfileSnapshot &value) { return vector_array(value.depth); })
        .def_property_readonly("wavenumbers", [](const ModeProfileSnapshot &value) { return vector_array(value.wavenumbers); })
        .def_property_readonly("group_velocity", [](const ModeProfileSnapshot &value) { return vector_array(value.group_velocity); })
        .def_property_readonly("mode_shapes", [](const ModeProfileSnapshot &value) { return matrix_array(value.mode_shapes); });

    py::class_<KernelInterface>(module, "KernelInterface")
        .def(py::init<>())
        .def(py::init<ThreadPool &>(), py::arg("pool"), py::keep_alive<1, 2>())
        .def("setNumThreads", &KernelInterface::setNumThreads)
        .def("setThreadPool", &KernelInterface::setThreadPool, py::keep_alive<1, 2>())
        .def("getNumThreads", &KernelInterface::getNumThreads)
        .def("getHardwareThreads", &KernelInterface::getHardwareThreads)
        .def("from_env", &KernelInterface::from_env)
        .def("from_json", &KernelInterface::from_json)
        .def("to_json", &KernelInterface::to_json)
        .def("to_json_string", &KernelInterface::to_json_string)
        .def("set_Title", [](KernelInterface &self, const std::string &title) {
            std::string copy = title;
            self.set_Title(copy);
        })
        .def("set_Freq", &KernelInterface::set_Freq)
        .def("set_freqvec", &KernelInterface::set_freqvec)
        .def("set_SSP", &KernelInterface::set_SSP)
        .def("set_AttenUnit", &KernelInterface::set_AttenUnit)
        .def("set_Sz", [](KernelInterface &self, const Eigen::VectorXd &values) { self.set_Sz(values); })
        .def("set_Sz", [](KernelInterface &self, double start, double end, int count) { self.set_Sz(start, end, count); })
        .def("set_Rr", [](KernelInterface &self, const Eigen::VectorXd &values) { self.set_Rr(values); })
        .def("set_Rr", [](KernelInterface &self, double start, double end, int count) { self.set_Rr(start, end, count); })
        .def("set_Rz", [](KernelInterface &self, const Eigen::VectorXd &values) { self.set_Rz(values); })
        .def("set_Rz", [](KernelInterface &self, double start, double end, int count) { self.set_Rz(start, end, count); })
        .def("set_Ro", [](KernelInterface &self, const Eigen::VectorXd &values) { self.set_Ro(values); })
        .def("set_Ro", [](KernelInterface &self, double start, double end, int count) { self.set_Ro(start, end, count); })
        .def("set_RProf", [](KernelInterface &self, const Eigen::VectorXd &values) { self.set_RProf(values); })
        .def("set_RProf", [](KernelInterface &self, double start, double end, int count) { self.set_RProf(start, end, count); })
        .def("set_MLimit", &KernelInterface::set_MLimit)
        .def("set_cPhase", &KernelInterface::set_cPhase)
        .def("set_GridType", &KernelInterface::set_GridType)
        .def("set_Rmax", &KernelInterface::set_Rmax)
        .def("set_SourceType", &KernelInterface::set_SourceType)
        .def("set_RunMode", &KernelInterface::set_RunMode)
        .def("set_CoherenceType", &KernelInterface::set_CoherenceType)
        .def("set_ModeType", &KernelInterface::set_ModeType)
        .def("set_Velocity_enable", &KernelInterface::set_Velocity_enable)
        .def("set_ReflCoef_Top", &KernelInterface::set_ReflCoef_Top)
        .def("set_ReflCoef_Bottom", &KernelInterface::set_ReflCoef_Bottom)
        .def("set_SBP", &KernelInterface::set_SBP)
        .def("run", &KernelInterface::run, py::call_guard<py::gil_scoped_release>())
        .def("runEigen", &KernelInterface::runEigen, py::call_guard<py::gil_scoped_release>())
        .def("runField", &KernelInterface::runField, py::call_guard<py::gil_scoped_release>())
        .def("clearResults", &KernelInterface::clearResults)
        .def("free", &KernelInterface::free)
        .def("export_result", &KernelInterface::export_result)
        .def("export_mod", &KernelInterface::export_mod)
        .def("export_shd", &KernelInterface::export_shd)
        .def("get_pressure", [](const KernelInterface &self) { return field_array(self.getPressureCopy()); })
        .def("get_vertical_velocity", [](const KernelInterface &self) { return field_array(self.getVerticalVelocityCopy()); })
        .def("get_horizontal_velocity", [](const KernelInterface &self) { return field_array(self.getHorizontalVelocityCopy()); })
        .def("get_pressure_snapshot", &KernelInterface::getPressureCopy)
        .def("get_vertical_velocity_snapshot", &KernelInterface::getVerticalVelocityCopy)
        .def("get_horizontal_velocity_snapshot", &KernelInterface::getHorizontalVelocityCopy)
        .def("get_modes", &KernelInterface::getModesCopy);

    module.attr("Interface") = module.attr("KernelInterface");
}
