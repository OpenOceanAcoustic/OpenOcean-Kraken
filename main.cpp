#include "src/algorithm/run.h"

int main()
{
    #ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
    #endif
    run();
    return 0;
}