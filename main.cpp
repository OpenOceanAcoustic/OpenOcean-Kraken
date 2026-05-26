#include <iostream>
#include "OpenOceanKrakenInterface.h"
int main()
{

    #ifdef _WIN32
    // 设置控制台输出为 UTF-8
    system("chcp 65001 >nul");
    #endif

    try {
        OpenOceanKraken::Interface kraken_interface;

        // 测试env文件转换功能
        std::string envFilePath = "E:\\my_project\\env_for_OOKtest\\stepK_rd.env";
        std::string jsonSavePath = "E:\\my_project\\OOK\\default_params.json";
        std::cout << "正在加载 env 文件: " << envFilePath << std::endl;
        bool success = kraken_interface.from_env(envFilePath);
        if (!success) {
            std::cerr << "env 加载失败" << std::endl;
            return 1;
        }

        // 将 env 参数转为 JSON 并保存
        kraken_interface.to_json(jsonSavePath);

    } catch (const std::exception& e) {
        std::cerr << "运行出错: " << e.what() << std::endl;
        return 1;
    }


    // kraken_interface.run();
    // kraken_interface.export_result("result");
    // bool success = kraken_interface.to_json("result.json"); 
    // if (success)
    // {
    //     std::cout << "JSON export successful" << std::endl;
    // }
    // else
    // {
    //     std::cout << "JSON export failed" << std::endl;
    // }

    //    创建接口实例
    // OpenOceanKraken::Interface kraken_interface;


    return 0;
}