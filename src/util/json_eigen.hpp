#ifndef JSON_EIGEN_HPP
#define JSON_EIGEN_HPP
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include <stdexcept>

// 辅助 trait：判断是否为 Eigen 稠密矩阵/向量
namespace details {
    template<typename T>
    struct is_eigen_matrix : std::false_type {};

    template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
    struct is_eigen_matrix<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> : std::true_type {};
} // namespace detail

// 通用 adl_serializer 特化
namespace nlohmann {

template<typename T>
struct adl_serializer<T, std::enable_if_t<details::is_eigen_matrix<T>::value>> {
    using Scalar = typename T::Scalar;
    
    // === to_json: 支持向量（1D）和矩阵（2D）===
    static void to_json(json& j, const T& mat) {
        if constexpr (T::RowsAtCompileTime == 1 || T::ColsAtCompileTime == 1) {
            // 向量：转为一维数组 [a, b, c]
            std::vector<Scalar> vec(mat.data(), mat.data() + mat.size());
            j = std::move(vec);
        } else {
            // 矩阵：转为二维数组 [[a,b],[c,d]]
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

    // === from_json: 从 JSON 恢复 Eigen 对象 ===
    static void from_json(const json& j, T& mat) {
        if (j.is_array()) {
            if (!j.empty() && j.front().is_array()) {
                // 二维数组 → 矩阵
                auto rows = static_cast<int>(j.size());
                auto cols = static_cast<int>(j[0].size());

                // 检查固定大小是否匹配
                if constexpr (T::RowsAtCompileTime != Eigen::Dynamic) {
                    if (rows != T::RowsAtCompileTime)
                        throw std::runtime_error("Row count mismatch");
                }
                if constexpr (T::ColsAtCompileTime != Eigen::Dynamic) {
                    if (cols != T::ColsAtCompileTime)
                        throw std::runtime_error("Column count mismatch");
                }

                mat.resize(rows, cols);
                for (int i = 0; i < rows; ++i) {
                    if (j[i].size() != static_cast<size_t>(cols))
                        throw std::runtime_error("Inconsistent column size");
                    for (int k = 0; k < cols; ++k) {
                        mat(i, k) = j[i][k].get<Scalar>();
                    }
                }
            } else {
                // 一维数组 → 向量
                auto size = static_cast<int>(j.size());
                if constexpr (T::SizeAtCompileTime != Eigen::Dynamic) {
                    if (size != T::SizeAtCompileTime)
                        throw std::runtime_error("Vector size mismatch");
                }
                mat.resize(size); // 对向量，resize(单参数) 设置总元素数
                for (int i = 0; i < size; ++i) {
                    mat(i) = j[i].get<Scalar>();
                }
            }
        } else {
            throw std::runtime_error("Expected JSON array for Eigen matrix/vector");
        }
    }
};

} // namespace nlohmann

#endif // JSON_EIGEN_HPP