#include "RootFinderBrent.h"


void ZBRENTX(double &x, double &a, double &b, const double t, 
        int& iset, int &mode, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, bool coutmodes, int &modeCount,
    std::string &errorMessage, FunctType funct) {
    int iExpA, iExpB, iExpC;
    const double MACHEP = 1.0E-16;
    const double TEN = 10.0;
    
    double fa, fb, fc, F1, F2, C, D, E, M, P, Q, R, S, TOL;
    
    errorMessage = "";
    
    // 计算区间端点的函数值
    funct(iset, mode, a, fa, iExpA, kramtrx, params, EVMat, false, modeCount);
    funct(iset, mode, b, fb, iExpB, kramtrx, params, EVMat, false, modeCount);
    
    // 检查区间端点函数值是否异号
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        errorMessage = "Function sign is the same at the interval endpoints";
        return;
    }
    
    // 内部根初始化
    C = a;
    fc = fa;
    iExpC = iExpA;
    E = b - a;
    D = E;
    // 调整值，确保b是当前最佳近似点
    if (iExpA < iExpB) {
        F1 = fc * std::pow(TEN, iExpC - iExpB);
        F2 = fb;
    } else {
        F1 = fc;
        F2 = fb * std::pow(TEN, iExpB - iExpC);
    }

    while (true) {

        
        if (std::abs(F1) < std::abs(F2)) {
            // 交换a, b, c和相应的函数值
            a = b;
            b = C;
            C = a;
            
            fa = fb;
            iExpA = iExpB;
            
            fb = fc;
            iExpB = iExpC;
            
            fc = fa;
            iExpC = iExpA;
            

        }
        
        // 计算容差和中点
        TOL = 2.0 * MACHEP * std::abs(b) + t;
        M = 0.5 * (C - b);
        
        // 检查是否收敛
        if (std::abs(M) <= TOL || fb == 0.0) {
            break;
        }
        
        // 检查是否需要强制二分法
        if (iExpA < iExpB) {
            F1 = fa * std::pow(TEN, iExpA - iExpB);
            F2 = fb;
        } else {
            F1 = fa;
            F2 = fb * std::pow(TEN, iExpB - iExpA);
        }
        
        if (std::abs(E) < TOL || std::abs(F1) <= std::abs(F2)) {
            // 使用二分法
            E = M;
            D = E;
        } else {
            // 计算插值步长
            S = fb / fa * std::pow(TEN, iExpB - iExpA);
            
            if (a == C) {
                // 线性插值
                P = 2.0 * M * S;
                Q = 1.0 - S;
            } else {
                // 逆二次插值
                Q = fa / fc * std::pow(TEN, iExpA - iExpC);
                R = fb / fc * std::pow(TEN, iExpB - iExpC);
                P = S * (2.0 * M * Q * (Q - R) - (b - a) * (R - 1.0));
                Q = (Q - 1.0) * (R - 1.0) * (S - 1.0);
            }
            
            // 保证Q为正
            if (P > 0.0) {
                Q = -Q;
            } else {
                P = -P;
            }
            
            S = E;
            E = D;
            
            // 检查步长是否可接受
            if (2.0 * P < 3.0 * M * Q - std::abs(TOL * Q) && P < std::abs(0.5 * S * Q)) {
                D = P / Q;
            } else {
                E = M;
                D = E;
            }
        }
        
        // 更新a和b
        a = b;
        fa = fb;
        iExpA = iExpB;
        
        if (std::abs(D) > TOL) {
            b += D;
        } else {
            if (M > 0.0) {
                b += TOL;
            } else {
                b -= TOL;
            }
        }
        
        // 计算新的b点函数值
        funct(iset, mode, b, fb, iExpB, kramtrx, params, EVMat, false, modeCount);
        
        // 更新c点
        if ((fb > 0.0) == (fc > 0.0)) {
            C = a;
            fc = fa;
            iExpC = iExpA;
            E = b - a;
            D = E;

            // 调整值，确保b是当前最佳近似点
            if (iExpA < iExpB) {
                F1 = fc * std::pow(TEN, iExpC - iExpB);
                F2 = fb;
            } else {
                F1 = fc;
                F2 = fb * std::pow(TEN, iExpB - iExpC);
            }
        }
    }
    
    // 返回找到的零点
    x = b;
}