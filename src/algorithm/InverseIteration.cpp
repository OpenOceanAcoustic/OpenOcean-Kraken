#include "InverseIteration.h"



void InverseIterationD(int N, 
                       VectorXd& D, 
                       VectorXd& E, 
                       int& IERR, 
                       VectorXd& PhiVector) {
    IERR = 0;

    // 计算矩阵的无穷范数
    double norm = D.array().abs().sum() + E.segment(1, N-1).array().abs().sum();

    // 计算小量参数（保持原逻辑的100倍系数）
    double eps3 = 100.0 * std::numeric_limits<double>::epsilon() * norm;
    double uk = N;
    double eps4 = uk * eps3;
    uk = eps4 / std::sqrt(uk);

    // 局部变量和临时向量
    double u, uk_val, v, xu;
    Eigen::VectorXd RV1(N), RV2(N), RV3(N), RV4(N);

    // 消去过程（带交换）
    xu = 1.0;
    u = D(0);  // 0-based索引，对应Fortran的D(1)
    v = E(1);  // E(2) in Fortran

    for (int i = 1; i < N; ++i) {  // i对应Fortran的2~N（0-based转换）
        if (std::abs(E(i)) >= std::abs(u)) {
            xu = u / E(i);
            RV4(i) = xu;
            RV1(i-1) = E(i);
            RV2(i-1) = D(i);
            RV3(i-1) = E(i+1);
            u = v - xu * RV2(i-1);
            v = -xu * RV3(i-1);
        } else {
            xu = E(i) / u;
            RV4(i) = xu;
            RV1(i-1) = u;
            RV2(i-1) = v;
            RV3(i-1) = 0.0;
            u = D(i) - xu * v;
            v = E(i+1);
        }
    }
    if (u == 0.0) u = eps3;

    // 初始化剩余的临时向量元素
    if (N >= 2) RV3(N-2) = 0.0;  // 对应Fortran的RV3(N-1)
    RV1(N-1) = u;
    RV2(N-1) = 0.0;
    RV3(N-1) = 0.0;

    // 初始化特征向量
    PhiVector.setConstant(uk);

    // 反迭代主循环
    for (int iter = 0; iter < MAXIT; ++iter) {
        // 回代过程
        u = 0.0;
        v = 0.0;
        for (int i = N-1; i >= 0; --i) {  // 从最后一个元素向前
            double temp = PhiVector(i) - u * RV2(i) - v * RV3(i);
            PhiVector(i) = temp / RV1(i);
            v = u;
            u = PhiVector(i);
        }

        // 计算向量范数并检查收敛
        norm = PhiVector.array().abs().sum();
        if (norm >= 1.0) {
            return;  // 收敛，返回结果
        }

        // 缩放向量         
        xu = eps4 / norm;
        PhiVector *= xu;

        // 前向消去
        for (int i = 1; i < N; ++i) {  // i对应Fortran的2~N
            double u_val = PhiVector(i);

            // 检查是否在三角化过程中交换过行
            if (RV1(i-1) == E(i)) {
                u_val = PhiVector(i-1);
                PhiVector(i-1) = PhiVector(i);
            }
            PhiVector(i) = u_val - RV4(i) * PhiVector(i-1);
        }
    }

    // 若循环结束仍未收敛，设置错误码
    IERR = -1;
}