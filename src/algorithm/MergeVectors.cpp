
#include "MergeVectors.h"

void MergeVectors(VectorXd& x, VectorXd& y, VectorXd& z) {
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
    z.resize(iz);
    z = temp.head(iz);  // 赋值给输出向量
}