# include "OpenOceanKrakenInterface.h"
// #include "src/algorithm/run.h"
int main()
{   

    #ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
    #endif
    OpenOceanKraken::Interface kraken_interface;
    OpenOceanKraken::ssp::SSPLayer ssp_layer;
    OpenOceanKraken::ssp::Range_Independent_Area area;
    std::vector<OpenOceanKraken::ssp::Range_Independent_Area> areas;

    Eigen::VectorXd z(2);
    // z << 0,5000;

    // Eigen::VectorXd alphaR(2);
    // alphaR << 1500,1500;
    // Eigen::VectorXd alphaI(2);
    // alphaI << 0.8,0.8;
    // Eigen::VectorXd betaR(2);
    // betaR << 0,0;
    // Eigen::VectorXd betaI(2);
    // betaI << 0,0;
    // Eigen::VectorXd rho(2);
    // rho << 1,1;
    // ssp_layer.alphaR = alphaR;
    // ssp_layer.alphaI = alphaI;
    // ssp_layer.betaR = betaR;
    // ssp_layer.betaI = betaI;
    // ssp_layer.rho = rho;

    // area.addLayer(ssp_layer);
    // areas.push_back(area);
    // kraken_interface.set_SSP(areas);




    int numThreads = kraken_interface.getHardwareThreads();
    // numThreads = 1;
    kraken_interface.setNumThreads(numThreads);
    ThreadPool threadPool(numThreads);
    kraken_interface.setThreadPool(threadPool);
    kraken_interface.run();
    kraken_interface.export_result("result");
    bool success = kraken_interface.to_json("result.json"); 
    if (success)
    {
        std::cout << "JSON export successful" << std::endl;
    }
    else
    {
        std::cout << "JSON export failed" << std::endl;
    }

    // #ifdef _WIN32
    // // 设置控制台输出为 UTF-8
    // system("chcp 65001 >nul");
    // #endif
    // std::cout << "OpenOcean-Kraken" << std::endl;
    // kkc_interface kkc;
    // kkc.run();
    // std::string baseFilename = "result";
    // kkc.export_result(baseFilename);

    return 0;
}