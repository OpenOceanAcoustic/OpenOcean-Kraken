#include "json_eigen.hpp"
#include "OpenOceanKrakenParams.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <Eigen/Dense>

namespace OpenOceanKraken
{
    bool read_sbp_file(const std::string &envPath, OOK_parameters &params)
    {
        std::string sbpPath = envPath.substr(0, envPath.length() - 4) + ".sbp"; // 打开同名sbp文件
        std::ifstream sbpFile(sbpPath);
        if (!sbpFile.is_open())
        {
            std::ostringstream oss;
            oss << "打开 sbp 文件 " << sbpPath << " 失败" << std::endl;
            throw std::runtime_error(oss.str());
            return false;
        }

        sbpFile >> params.SBP.NSBPPts;
        std::vector<double> theta(params.SBP.NSBPPts);
        std::vector<double> pat(params.SBP.NSBPPts);
        for (int i = 0; i < params.SBP.NSBPPts; i++)
        {
            sbpFile >> theta[i] >> pat[i];
        }
        params.SBP.theta = Eigen::Map<Eigen::VectorXd>(theta.data(), theta.size());
        params.SBP.pat = Eigen::Map<Eigen::VectorXd>(pat.data(), pat.size());
        params.SBP.isSet = true;
        sbpFile.close();
        return true;
    }
    bool read_flp_file(const std::string &envPath, OOK_parameters &params)
    {
        // 生成对应的 flp 文件路径
        std::string flpPath;
        size_t dotPos = envPath.find_last_of('.'); // 找到最后一个.的位置

        if (dotPos != std::string::npos && envPath.substr(dotPos) == ".env")
        {
            flpPath = envPath.substr(0, dotPos) + ".flp";
        }
        else
        {
            // 无.env 后缀，直接追加.flp
            flpPath = envPath + ".flp";
        }

        // 读取flp文件
        std::ifstream flpFile(flpPath);
        if (!flpFile.is_open())
        {
            std::cerr << "打开 flp 文件 " << flpPath << " 失败" << std::endl;
            return false;
        }

        std::string flp_line;
        int flp_line_num = 0;
        std::vector<std::string> flp_lines;

        // 读取所有行并预处理
        while (std::getline(flpFile, flp_line))
        {
            flp_line_num++;
            // 移除注释（通常flp文件以!开头为注释）
            size_t comment_pos = flp_line.find('!');
            if (comment_pos != std::string::npos)
            {
                flp_line = flp_line.substr(0, comment_pos);
            }

            // 去除首尾空白字符
            flp_line.erase(0, flp_line.find_first_not_of(" \t\r\n"));
            flp_line.erase(flp_line.find_last_not_of(" \t\r\n") + 1);

            // 去掉单引号 '
            flp_line.erase(std::remove(flp_line.begin(), flp_line.end(), '\''), flp_line.end());

            if (!flp_line.empty())
            {
                flp_lines.push_back(flp_line);
            }
        }

        flpFile.close();
        if (flp_lines.empty())
        {
            std::cerr << "Flp 文件为空或只包含注释: " << flpPath << std::endl;
            return false;
        }

        // 解析flp文件
        try
        {
            size_t line_idx = 0;
            // env的第1行：标题行无需读取
            line_idx++;
            // 第2行：SOURCE OPTION 声源选项
            if (line_idx < flp_lines.size())
            {
                char option1, option2, option3, option4;
                std::istringstream iss(flp_lines[line_idx++]);
                if (iss >> option1)
                {
                    if (option1 == 'R')
                    {
                        params.SourceType = Source_Mode::MODE_R_Point;
                    }
                    else if (option1 == 'X')
                    {
                        params.SourceType = Source_Mode::MODE_X_Line;
                    }
                }
                if (iss >> option2)
                {
                    if (option2 == 'A')
                    {
                        params.modeType = ModeType::Adiabatic;
                    }
                    else if (option2 == 'C')
                    {
                        params.modeType = ModeType::Couple;
                    }
                }
                if (iss >> option3)
                {
                    if (option3 == '*')
                    {
                        read_sbp_file(envPath, params);
                        if (!read_sbp_file(envPath, params))
                        {
                            std::cerr << "警告: sbp 文件解析失败，使用已解析的 env 参数继续" << std::endl;
                        }
                    }
                    else if (option3 == 'O')
                    {
                        params.SBP.isSet = false;
                    }
                }
                else
                {
                    params.SBP.isSet = false;
                }
                if (iss >> option4)
                {
                    if (option4 == 'C')
                    {
                        params.coherenceType = CoherenceType::Coherent;
                    }
                    else if (option4 == 'I')
                    {
                        params.coherenceType = CoherenceType::Incoherent;
                    }
                }
                else
                {
                    params.coherenceType = CoherenceType::Coherent;
                }
            }
            // 第3~5行：模态数量跳过，声速剖面个数和范围在这里kraken也不解析，而且在env中已经赋值过了
            for (int i = 0; i < 3 && line_idx < flp_lines.size(); ++i)
                line_idx++;

            // 第6行：接收器水平个数
            if (line_idx < flp_lines.size())
            {
                std::istringstream iss(flp_lines[line_idx++]);
                iss >> params.Pos.NRr;
            }

            // 第7行：接收器距离范围
            Eigen::VectorXd Rr;
            std::vector<double> temp_rr; // 临时存所有接收器距离范围
            if (line_idx < flp_lines.size())
            {
                double val;
                std::istringstream iss(flp_lines[line_idx++]);
                while (iss >> val)
                {
                    temp_rr.push_back(val); // 先把一行所有值读完
                }
            }
            Rr = Eigen::Map<Eigen::VectorXd>(temp_rr.data(), temp_rr.size());
            if (Rr.size() != params.Pos.NRr)
            {
                params.Pos.is_Linspace_Rr = true;
                params.Pos.Rr = Rr * 1e3;
            }
            else
            {
                params.Pos.is_Linspace_Rr = false;
                params.Pos.Rr = Rr;
            }

            // 第8~11行:NSz,Sz,NRz,Rz，均已经赋值过此处跳过
            for (int i = 0; i < 4 && line_idx < flp_lines.size(); ++i)
                line_idx++;

            // 第12行：接收器垂直个数
            if (line_idx < flp_lines.size())
            {
                std::istringstream iss(flp_lines[line_idx++]);
                iss >> params.Pos.NRo;
            }

            // 第13行：接收器在水平方向的偏移（倾斜阵）
            Eigen::VectorXd Ro;
            std::vector<double> temp_ro; // 临时存所有接收器距离范围
            if (line_idx < flp_lines.size())
            {
                double val;
                std::istringstream iss(flp_lines[line_idx++]);
                while (iss >> val)
                {
                    temp_ro.push_back(val); // 先把一行所有值读完
                }
            }
            Ro = Eigen::Map<Eigen::VectorXd>(temp_ro.data(), temp_ro.size());
            if (Ro.size() == 1 && Ro(0) == 0.0)
            {
                Ro.resize(2);
                Ro << 0.0, 0.0;
                params.Pos.Ro = Ro;
                params.Pos.is_Linspace_Ro = true;
            }
            else if (Ro.size() != params.Pos.NRo)
            {
                params.Pos.is_Linspace_Ro = true;
                params.Pos.Ro = Ro;
            }
            else
            {
                params.Pos.is_Linspace_Ro = false;
                params.Pos.Ro = Ro;
            }

            params.Pos.GridType = Grid_Mode::MODE_R_Rectangular;
            params.runMode = Run_Mode::MODE_B_Both;
        }
        catch (const std::exception &e)
        {
            std::cerr << "解析 flp 文件时出错: " << e.what() << std::endl;
            return false;
        }

        return true;
    }

    bool read_refCoef_file(const std::string &envPath, OOK_parameters &params, std::string pattern)
    {
        int trc_err_line = 0;                                                        // 行数计数器，出错时抛出行数
        std::string refCoefPath = envPath.substr(0, envPath.length() - 4) + pattern; // 打开同名trc文件
        std::ifstream refCoefFile(refCoefPath);
        if (!refCoefFile.is_open())
        {
            std::ostringstream oss;
            oss << "无法打开" << pattern << "文件: " << refCoefPath;
            throw std::runtime_error(oss.str());
        }
        try
        {
            int Npts;
            refCoefFile >> Npts;
            std::vector<double> thetas(Npts), Rs(Npts), phis(Npts);
            if (pattern == ".trc")
            {
                // .trc 文件是顶部反射系数文件
                params.ReflectionCoef.RTop.resize(Npts);
                for (int i = 0; i < Npts; i++)
                {
                    refCoefFile >> params.ReflectionCoef.RTop(i).theta;
                    refCoefFile >> params.ReflectionCoef.RTop(i).R;
                    refCoefFile >> params.ReflectionCoef.RTop(i).phi;
                }
            }
            else if (pattern == ".brc")
            {
                // .brc 文件是底部反射系数文件
                params.ReflectionCoef.RBot.resize(Npts);
                for (int i = 0; i < Npts; i++)
                {
                    refCoefFile >> params.ReflectionCoef.RBot(i).theta;
                    refCoefFile >> params.ReflectionCoef.RBot(i).R;
                    refCoefFile >> params.ReflectionCoef.RBot(i).phi;
                }
            }
            else
            {
                std::cerr << "未知的文件类型: " << pattern << std::endl;
                return false;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "反射系数文件时格式错误:" << refCoefPath << std::endl;
            std::cerr << e.what() << '\n';
        }
        refCoefFile.close();
        return true;
    }

    // 解析单个环境的所有行
    static bool parse_single_env(const std::vector<std::string>& lines,
                                 size_t& line_idx,
                                 OOK_parameters& params,
                                 const std::string& envPath)
    {
        try
        {
            // 初始化 sspInput[0] 和其内部的 layers[0]
            params.sspInput.resize(1);
            params.sspInput[0].Range = 0.0;
            params.sspInput[0].layers.resize(1);

            // 第1行: 标题
            if (line_idx < lines.size())
            {
                params.Title = lines[line_idx++];
            }

            // 第2行: 频率
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                iss >> params.freqinfo.freq;
                params.freqinfo.Nfreq = 1;
            }

            // 第3行: 介质层数
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                int val;
                if (iss >> val)
                {
                    params.NMediaMax = static_cast<int>(val);
                }
                else
                {
                    std::cerr << "第3行: 介质层数格式错误" << std::endl;
                    return false;
                }
            }

            // 第4行：海面选项（TOP OPTION）
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                char option1, option2, option3, option4;

                if (iss >> option1 >> option2 >> option3)
                {
                    // 1. 第1个字母：声速剖面插值类型
                    switch (option1)
                    {
                    case 'C':
                        params.sspInput[0].SSPType = SSP_Mode::MODE_C_cLinear;
                        break;
                    case 'N':
                        params.sspInput[0].SSPType = SSP_Mode::MODE_N_n2Linear;
                        break;
                    case 'P':
                        params.sspInput[0].SSPType = SSP_Mode::MODE_P_cPCHIP;
                        break;
                    case 'S':
                        params.sspInput[0].SSPType = SSP_Mode::MODE_S_cCubic;
                        break;
                    default:
                        params.sspInput[0].SSPType = SSP_Mode::MODE_C_cLinear;
                        break;
                    }

                    // 2. 第2个字母：海面边界条件 TOP
                    switch (option2)
                    {
                    case 'V':
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_V_Vacuum;
                        break;
                    case 'R':
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_R_Rigid;
                        break;
                    case 'A':
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_A_Half_space;
                        break;
                    case 'F':
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_F_File;
                        read_refCoef_file(envPath, params, ".trc");
                        break;
                    case 'G':
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_G_Grain;
                        break;
                    case 'P':
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_P_Precomputed;
                        break;
                    default:
                        params.sspInput[0].HSTop.BC = BC_Mode::MODE_V_Vacuum;
                        break;
                    }

                    // 3. 第3个字母：衰减系数单位
                    switch (option3)
                    {
                    case 'W':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_W_db_per_lambda;
                        break;
                    case 'M':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_M_dB_per_m;
                        break;
                    case 'N':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_N_Nepers_per_m;
                        break;
                    case 'Q':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_Q_Quality_Factor;
                        break;
                    case 'L':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_L_params_lose;
                        break;
                    case 'F':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_F_dB_per_m_kHz;
                        break;
                    case 'm':
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_m_dB_per_m;
                        break;
                    default:
                        params.AttenUnit.attnUnit = AttenuationUnit::MODE_W_db_per_lambda;
                        break;
                    }

                    // 4. 第4个字母：附加体积衰减选项或计算选项(可选)
                    if (iss >> option4)
                    {
                        switch (option4)
                        {
                        case 'T':
                            params.AttenUnit.absModel = OceanAbsorptionModel::Thorpe;
                            break;
                        case 'F':
                            params.AttenUnit.absModel = OceanAbsorptionModel::FrancGarr;
                            break;
                        default:
                            params.AttenUnit.absModel = OceanAbsorptionModel::None;
                            break;
                        }
                    }
                    else
                    {
                        params.AttenUnit.absModel = OceanAbsorptionModel::None;
                    }
                }
                else
                {
                    std::cerr << "第4行: 海面选项格式错误" << std::endl;
                    return false;
                }
            }

            // 第5行：海水竖直网格层个数，界面粗糙度，海水深度
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                double nmesh, roughness, depth;
                if (iss >> nmesh >> roughness >> depth)
                {
                    params.sspInput[0].layers[0].nmesh = static_cast<int>(nmesh);
                    params.sspInput[0].layers[0].sigma = roughness;
                    params.sspInput[0].HSBot.Depth = depth;
                }
                else
                {
                    std::cerr << "第5行: 网格层个数、粗糙度、深度格式错误" << std::endl;
                    return false;
                }
            }

            // 声速剖面
            struct Point
            {
                double z = 0, alphaR = 0, betaR = 0, rho = 1, alphaI = 0, betaI = 0;
            };
            std::vector<Point> points;
            Point last_p;

            while (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                Point p = last_p;
                int count = 0;
                if (iss >> p.z)
                    count++;
                if (iss >> p.alphaR)
                    count++;
                if (iss >> p.betaR)
                    count++;
                if (iss >> p.rho)
                    count++;
                if (iss >> p.alphaI)
                    count++;
                if (iss >> p.betaI)
                    count++;

                if (count < 2)
                {
                    if (!iss)
                        line_idx--;
                    break;
                }
                points.push_back(p);
                last_p = p;
            }
            const int npts = static_cast<int>(points.size());

            params.sspInput[0].layers[0].z.resize(npts);
            params.sspInput[0].layers[0].alphaR.resize(npts);
            params.sspInput[0].layers[0].betaR.resize(npts);
            params.sspInput[0].layers[0].rho.resize(npts);
            params.sspInput[0].layers[0].alphaI.resize(npts);
            params.sspInput[0].layers[0].betaI.resize(npts);

            for (int i = 0; i < npts; ++i)
            {
                params.sspInput[0].layers[0].z[i] = points[i].z;
                params.sspInput[0].layers[0].alphaR[i] = points[i].alphaR;
                params.sspInput[0].layers[0].betaR[i] = points[i].betaR;
                params.sspInput[0].layers[0].rho[i] = points[i].rho;
                params.sspInput[0].layers[0].alphaI[i] = points[i].alphaI;
                params.sspInput[0].layers[0].betaI[i] = points[i].betaI;
            }

            params.sspInput[0].layers[0].npts = npts;
            if (params.sspInput[0].layers[0].nmesh == 0)
            {
                params.sspInput[0].layers[0].nmesh = static_cast<int>(points.size());
            }

            // 海底半空间
            char bottomType;
            double dummy;
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                iss >> bottomType;
                switch (bottomType)
                {
                case 'A':
                    params.sspInput[0].HSBot.BC = BC_Mode::MODE_A_Half_space;
                    break;
                case 'V':
                    params.sspInput[0].HSBot.BC = BC_Mode::MODE_V_Vacuum;
                    break;
                case 'R':
                    params.sspInput[0].HSBot.BC = BC_Mode::MODE_R_Rigid;
                    break;
                case 'F':
                    params.sspInput[0].HSBot.BC = BC_Mode::MODE_F_File;
                    read_refCoef_file(envPath, params, ".brc");
                    break;
                case 'G':
                    params.sspInput[0].HSBot.BC = BC_Mode::MODE_G_Grain;
                    break;
                case 'P':
                    params.sspInput[0].HSBot.BC = BC_Mode::MODE_P_Precomputed;
                    break;
                }
            }

            if (params.sspInput[0].HSBot.BC != BC_Mode::MODE_F_File)
            {
                Point bottom_p;
                if (line_idx < lines.size())
                {
                    std::istringstream iss(lines[line_idx++]);
                    iss >> params.sspInput[0].HSBot.Depth;
                    iss >> params.sspInput[0].HSBot.alphaR;
                    iss >> params.sspInput[0].HSBot.alphaI;
                    iss >> params.sspInput[0].HSBot.betaR;
                    iss >> params.sspInput[0].HSBot.betaI;
                    iss >> params.sspInput[0].HSBot.rho;
                }
            }

            // 相速度
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                iss >> params.cLow >> params.cHigh;
            }

            // 最大距离
            double Rmax;
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                iss >> Rmax;
                params.Rmax = Rmax * 1e3;
            }

            // 声源个数
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                iss >> params.Pos.NSz;
            }

            // 声源深度
            std::vector<double> temp_sz;
            if (line_idx < lines.size())
            {
                double val;
                std::istringstream iss(lines[line_idx++]);
                while (iss >> val)
                {
                    temp_sz.push_back(val);
                }
            }
            params.Pos.Sz = Eigen::Map<Eigen::VectorXd>(temp_sz.data(), temp_sz.size());
            if (params.Pos.Sz.size() != params.Pos.NSz)
            {
                params.Pos.is_Linspace_Sz = true;
            }
            else
            {
                params.Pos.is_Linspace_Sz = false;
            }

            // 接收器深度个数
            if (line_idx < lines.size())
            {
                std::istringstream iss(lines[line_idx++]);
                iss >> params.Pos.NRz;
            }

            // 接收器深度
            std::vector<double> temp_rz;
            if (line_idx < lines.size())
            {
                double val;
                std::istringstream iss(lines[line_idx++]);
                while (iss >> val)
                {
                    temp_rz.push_back(val);
                }
            }

            params.Pos.Rz = Eigen::Map<Eigen::VectorXd>(temp_rz.data(), temp_rz.size());
            if (params.Pos.Rz.size() != params.Pos.NRz)
            {
                params.Pos.is_Linspace_Rz = true;
            }
            else
            {
                params.Pos.is_Linspace_Rz = false;
            }

            read_flp_file(envPath, params);
            if (!read_flp_file(envPath, params))
            {
                std::cerr << "警告: flp 文件解析失败，使用已解析的 env 参数继续" << std::endl;
            }

            return true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "解析 env 文件时出错，行号 " << line_idx << ": " << e.what() << std::endl;
            return false;
        }
    }

    bool read_env_file(const std::string &envPath, OOK_parameters &params)
    {
        std::ifstream file(envPath);
        if (!file.is_open())
        {
            std::cerr << "无法打开 env 文件: " << envPath << std::endl;
            return false;
        }

        std::string line;
        int line_num = 0;
        std::vector<std::string> lines;

        // 读取所有行并预处理
        while (std::getline(file, line))
        {
            line_num++;

            // 移除注释（通常env文件以!开头为注释）
            size_t comment_pos = line.find('!');
            if (comment_pos != std::string::npos)
            {
                line = line.substr(0, comment_pos);
            }

            // 去除首尾空白字符
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            // 去掉单引号 '
            line.erase(std::remove(line.begin(), line.end(), '\''), line.end());

            if (!line.empty())
            {
                lines.push_back(line);
            }
        }

        file.close();

        if (lines.empty())
        {
            std::cerr << "Env 文件为空或只包含注释: " << envPath << std::endl;
            return false;
        }

        // 解析env文件
        size_t line_idx = 0;
        params.sspInput.clear();

        while (line_idx < lines.size())
        {
            OOK_parameters single_params;
            if (!parse_single_env(lines, line_idx, single_params, envPath))
            {
                return false;
            }

            // 第一个环境：拷贝所有字段
            if (params.sspInput.empty())
            {
                params = single_params;
                params.sspInput.clear();
            }

            // 每个环境追加 sspInput
            params.sspInput.push_back(single_params.sspInput[0]);
        }

        return true;
    }
}