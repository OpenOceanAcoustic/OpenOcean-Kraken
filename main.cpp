# include "kkc_interface.h"
// #include "src/algorithm/run.h"

int main()
{
    #ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
    #endif
    // run();
    kkc_interface kkc;
    kkc.run();

    return 0;
}