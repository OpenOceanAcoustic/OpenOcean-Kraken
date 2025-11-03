
#include "MergeVectors.h"

void MergeVectors(VectorXd& x, VectorXd& y, VectorXd& z, int& NzTab) {
    int ix = 0, iy = 0, iz = 0;       // x的索引（0基）
    int Nx = x.size();  // x的长度
    int Ny = y.size();  // y的长度

    // 临时存储合并结果（动态扩容）
    VectorXd temp;
    temp.resize(Nx + Ny);  // 预分配内存，提升效率

    // 循环合并两个向量
    while (ix < Nx || iy < Ny) {
        double current;

        if (iy >= Ny) {
            // y已遍历完，取x的元素
            current = x(ix++);
        } else if (ix >= Nx) {
            // x已遍历完，取y的元素
            current = y(iy++);
        } else if (x(ix) <= y(iy)) {
            // x的当前元素更小，取x的元素
            current = x(ix++);
        } else {
            // y的当前元素更小，取y的元素
            current = y(iy++);
        }

        // 检查是否为近似重复元素
        if (iz > 0) {
            double last = temp(iz - 1);
            // 与Fortran的EPSILON(z)对应，使用double的机器epsilon
            if (std::abs(current - last) < 100.0 * std::numeric_limits<double>::epsilon()) {
                continue;  // 重复元素，不加入临时向量
            }
        }

        temp(iz++) = current;           // 添加元素
    }
    NzTab = iz;
    z.resize(iz);
    z = temp.head(iz);  // 赋值给输出向量
}


void Weight_dble(VectorXd& x, int Nx, 
                 VectorXd& xTab, int NxTab, 
                 VectorXd& w, VectorXi& Ix) {
    // 快速返回：如果插值点只有一个
    if (Nx == 1) {
        w(0) = 0.0;
        Ix(0) = 0;  // 0-based索引
        return;
    }

    int L = 0;  // 初始索引（0-based）

    // 为每个需要计算权重的点循环
    for (int IxTab = 0; IxTab < NxTab; ++IxTab) {
        // 搜索满足[x(L), x(L+1)]包含xTab(IxTab)的索引L
        // 循环条件：当前点大于右邻点，且未到达倒数第二个点（避免L+1越界）
        while (xTab(IxTab) > x(L + 1) && L < Nx - 2) {
            L++;
        }

        // 记录0-based索引和插值权重
        Ix(IxTab) = L;
        // 权重计算：(目标点 - 左端点) / (右端点 - 左端点)
        w(IxTab) = (xTab(IxTab) - x(L)) / (x(L + 1) - x(L));

        // 原注释中的特殊处理代码（底跟踪接收器）
        // if (w(IxTab) != 0.0) {
        //     Ix(IxTab) = std::min(1, NxTab - 1);  // 0-based的第二个元素索引为1
        //     w(IxTab) = 0.0;
        // }
    }
}