#ifndef JSON_EIGEN_HPP
#define JSON_EIGEN_HPP
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <stdexcept>

namespace nlohmann {
// 辅助 trait：判断是否为 Eigen 稠密矩阵/向量
namespace details {
    template<typename T>
    struct is_eigen_matrix : std::false_type {};

    template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    struct is_eigen_matrix<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> : std::true_type {};
} // namespace detail

// 自由函数模板：适用于任何 BasicJsonType（包括 json 和 ordered_json）
template<typename BasicJsonType, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void to_json(BasicJsonType& j,
             const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& mat) {
    using MatrixType = Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>;
    if constexpr (MatrixType::RowsAtCompileTime == 1 || MatrixType::ColsAtCompileTime == 1) {
        std::vector<Scalar> vec(mat.data(), mat.data() + mat.size());
        j = std::move(vec);
    } else {
        std::vector<std::vector<Scalar>> mat2d;
        mat2d.reserve(mat.rows());
        for (int i = 0; i < mat.rows(); ++i) {
            std::vector<Scalar> row(mat.cols());
            for (int k = 0; k < mat.cols(); ++k) {
                row[k] = mat(i, k);
            }
            mat2d.push_back(std::move(row));
        }
        j = std::move(mat2d);
    }
}

template<typename BasicJsonType, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void from_json(const BasicJsonType& j,
               Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& mat) {
    using MatrixType = Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>;
    if (j.is_array()) {
        if (!j.empty() && j.front().is_array()) {
            if constexpr (MatrixType::IsVectorAtCompileTime) {
                throw std::runtime_error("Expected 1D JSON array for Eigen vector");
            } else {
                auto rows = static_cast<int>(j.size());
                auto cols = static_cast<int>(j[0].size());

                if constexpr (MatrixType::RowsAtCompileTime != Eigen::Dynamic) {
                    if (rows != MatrixType::RowsAtCompileTime)
                        throw std::runtime_error("Row count mismatch");
                }
                if constexpr (MatrixType::ColsAtCompileTime != Eigen::Dynamic) {
                    if (cols != MatrixType::ColsAtCompileTime)
                        throw std::runtime_error("Column count mismatch");
                }

                mat.resize(rows, cols);
                for (int i = 0; i < rows; ++i) {
                    if (j[i].size() != static_cast<size_t>(cols))
                        throw std::runtime_error("Inconsistent column size");
                    for (int k = 0; k < cols; ++k) {
                        mat(i, k) = j[i][k].template get<Scalar>();
                    }
                }
            }
        } else {
            if constexpr (!MatrixType::IsVectorAtCompileTime) {
                throw std::runtime_error("Expected 2D JSON array for Eigen matrix");
            } else {
                auto size = static_cast<int>(j.size());
                if constexpr (MatrixType::SizeAtCompileTime != Eigen::Dynamic) {
                    if (size != MatrixType::SizeAtCompileTime)
                        throw std::runtime_error("Vector size mismatch");
                }
                mat.resize(size);
                for (int i = 0; i < size; ++i) {
                    mat(i) = j[i].template get<Scalar>();
                }
            }
        }
    } else {
        throw std::runtime_error("Expected JSON array for Eigen matrix/vector");
    }
}
}
#endif // JSON_EIGEN_HPP
