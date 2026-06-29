#include <iostream>
#include "OpenOceanKrakenInterface.h"
#include "ThreadPool.h"

int main()
{

    #ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
    #endif

    try {
        OpenOceanKraken::Interface kraken_interface;

        // 测试env文件转换功能
        std::string envName = "arcticK";
        std::string envFilePath = "E:\\my_project\\env_for_OOKtest\\" + envName + ".env";
        std::string resultSavePath = "E:\\my_project\\test_in_matlab\\env_test_OOK\\";
        std::cout << "正在加载 env 文件: " << envFilePath << std::endl;
        bool success = kraken_interface.from_env(envFilePath);
        if (!success) {
            std::cerr << "env 加载失败" << std::endl;
            return 1;
        }

        // 将 env 参数转为 JSON 并保存
        kraken_interface.to_json(resultSavePath + envName + ".json");
        std::cout << "JSON export successful: " << resultSavePath + envName + ".json" << std::endl;

        // 设置线程池和线程数
        int numThreads = kraken_interface.getHardwareThreads();
        kraken_interface.setNumThreads(numThreads);
        ThreadPool threadPool(numThreads);
        kraken_interface.setThreadPool(threadPool);
        // std::cout << "向量：\n" << kraken_interface.getParams().Pos.Rr << std::endl;

        // 运行声场计算
        kraken_interface.run();
        std::cout << "环境运行完成" << std::endl;
        // 导出声场
        kraken_interface.export_shd(resultSavePath + envName, 1);
        std::cout << "声场结果导出完成" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "运行出错: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}