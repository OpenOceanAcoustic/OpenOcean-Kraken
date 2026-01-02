#include "kkc_interface.h"
// #include "src/algorithm/run.h"

int main(int argc, char* argv[])
{
#ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
#endif
    // run();
    // 解析可选参数
    int numThreads = -1;  
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "-t" && i + 1 < argc)
        {
            try
            {
                numThreads = std::stoi(argv[++i]);
                if (numThreads <= 0)
                {
                    std::cerr << "错误: 线程数必须大于 0" << std::endl;
                    return 1;
                }
            }
            catch (...)
            {
                std::cerr << "错误: 线程数必须是有效整数" << std::endl;
                return 1;
            }
        }
        else
        {
            std::cerr << "未知参数: " << arg << std::endl;
            return 1;
        }
    }

    kkc_interface kkc;
    if (numThreads <= 0) {
        numThreads = kkc.getHardwareThreads();
    }
    kkc.setNumThreads(numThreads);
    ThreadPool threadPool(numThreads);
    kkc.setThreadPool(threadPool);
    kkc.run();
    std::string baseFilename = "result";
    kkc.export_result(baseFilename);

    return 0;
}