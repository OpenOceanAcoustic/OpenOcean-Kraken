# include "kkc_interface.h"
// #include "src/algorithm/run.h"

int main()
{
    #ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
    #endif
    std::cout << "OpenOcean-Kraken" << std::endl;
    kkc_interface kkc;
    kkc.run();
    std::string baseFilename = "result";
    kkc.export_result(baseFilename);

    return 0;
}