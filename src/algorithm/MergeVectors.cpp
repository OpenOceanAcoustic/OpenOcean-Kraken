
#include "MergeVectors.h"

namespace OpenOceanKraken
{
    void MergeVectors(const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &z, int &NzTab, Eigen::VectorXi &Ix, Eigen::VectorXi &Iy)
    {
        int ix = 0, iy = 0, iz = 0; // x的索引（0基）
        int Nx = x.size();          // x的长度
        int Ny = y.size();          // y的长度

        // 临时存储合并结果（动态扩容）
        Eigen::VectorXd temp;
        Eigen::VectorXi tempIx, tempIy;
        temp.resize(Nx + Ny);   // 预分配内存，提升效率
        tempIx.resize(Nx + Ny); // 存储x中元素在合并后向量中的索引
        tempIy.resize(Nx + Ny); // 存储y中元素在合并后向量中的索引

        // 初始化索引数组为-1，表示未分配
        tempIx.setConstant(-1);
        tempIy.setConstant(-1);

        // 循环合并两个向量
        while (ix < Nx || iy < Ny)
        {
            double current;
            int source;                     // 0表示来自x，1表示来自y
            int orig_ix = ix, orig_iy = iy; // 记录原始索引

            if (iy >= Ny)
            {
                // y已遍历完，取x的元素
                current = x(ix++);
                source = 0;
            }
            else if (ix >= Nx)
            {
                // x已遍历完，取y的元素
                current = y(iy++);
                source = 1;
            }
            else if (x(ix) <= y(iy))
            {
                // x的当前元素更小，取x的元素
                current = x(ix++);
                source = 0;
            }
            else
            {
                // y的当前元素更小，取y的元素
                current = y(iy++);
                source = 1;
            }

            // 检查是否为近似重复元素
            if (iz > 0)
            {
                double last = temp(iz - 1);
                // 与Fortran的EPSILON(z)对应，使用double的机器epsilon
                if (std::abs(current - last) < 100.0 * std::numeric_limits<double>::epsilon())
                {
                    // 重复元素，不加入临时向量，但仍需更新索引
                    if (source == 0)
                    {
                        // 如果是x中的元素且与前一个元素重复，则将其索引指向同一个位置
                        tempIx(orig_ix) = iz - 1;
                    }
                    else
                    {
                        // 如果是y中的元素且与前一个元素重复，则将其索引指向同一个位置
                        tempIy(orig_iy) = iz - 1;
                    }
                    continue;
                }
            }

            temp(iz) = current; // 添加元素
            if (source == 0)
            {
                tempIx(orig_ix) = iz; // 记录x中元素在合并后向量中的索引
            }
            else
            {
                tempIy(orig_iy) = iz; // 记录y中元素在合并后向量中的索引
            }
            iz++;
        }
        NzTab = iz;
        z.resize(iz);
        z = temp.head(iz); // 赋值给输出向量

        // 调整索引向量的大小并赋值
        Ix.resize(Nx);
        Iy.resize(Ny);
        Ix = tempIx.head(Nx);
        Iy = tempIy.head(Ny);
    }

    void Weight_dble(const Eigen::VectorXd &x, const int Nx,
                     const Eigen::VectorXd &xTab, int NxTab,
                     Eigen::VectorXd &w, Eigen::VectorXi &Ix)
    {
        // 快速返回：如果插值点只有一个
        if (Nx == 1)
        {
            for (int IxTab = 0; IxTab < NxTab; ++IxTab)
            {
                w(IxTab) = 0.0;
                Ix(IxTab) = 0; // 0-based索引
            }
            return;
        }

        int L = 0; // 初始索引（0-based）

        // 为每个需要计算权重的点循环
        for (int IxTab = 0; IxTab < NxTab; ++IxTab)
        {
            // 搜索满足[x(L), x(L+1)]包含xTab(IxTab)的索引L
            // 循环条件：当前点大于右邻点，且未到达倒数第二个点（避免L+1越界）
            while (xTab(IxTab) > x(L + 1) && L < Nx - 2)
            {
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
}
