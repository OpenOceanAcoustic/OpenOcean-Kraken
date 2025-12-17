#ifndef KKC_LOG_H
#define KKC_LOG_H
#include <iostream>
#include <thread>
#include <mutex>
#include <fstream>
#include <sstream>
class kkc_Log
{
public:
    // 枚举定义
    enum class LogLevel
    {
        Debug,
        Info,
        Warning,
        Error
    };

    enum class ErrorCode
    {
        Success,

        Blank_input,
        Wrong_value,
        Wrong_param_type,

        UnknownError
    };

    // 获取单例实例
    static kkc_Log &get_instance();

    // 流插入操作符重载，支持各种数据类型
    template <typename T>
    kkc_Log &operator<<(const T &content);
    // 特殊处理字符串字面量
    kkc_Log &operator<<(const char *content);
    // 换行操作
    kkc_Log &operator<<(std::ostream &(*manip)(std::ostream &));

    kkc_Log &operator<<(LogLevel level);
    kkc_Log &operator<<(ErrorCode code);

    void append(const std::string &content); // 添加内容到缓冲区
    void append(LogLevel level);
    void append(ErrorCode code);

    // 打印到控制台
    void print_to_Console();

    // 写入文件
    void print_to_file(const std::string &filename, bool append = true);

    // 刷新缓冲区
    void flush();

    // 设置默认文件名
    void set_default_file(const std::string &filename);

    // 获取当前缓冲区内容
    std::string get_buffer();
    bool raise_err_WrongType(std::string params_description, const std::string &filePath,
                             int line_num, bool raise_err); // 用于输出“空白输入/错误变量类型”的错误

    void raise_runtime_err_and_exit(std::string file_type, const std::string &filePath); // 用于抛出runtime err，但是在抛出前确保log打印到控制台上

private:
    kkc_Log();
    ~kkc_Log();
    std::stringstream ss_;
    std::mutex mtx_;
    std::string filename_;

    static std::string toString(LogLevel level);
    static std::string toString(ErrorCode code);
};

#endif // KKC_LOG_H
