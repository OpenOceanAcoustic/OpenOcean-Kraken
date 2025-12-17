#include "kkc_log.h"

// 获取单例实例
kkc_Log::kkc_Log() : filename_("default.log") {}

kkc_Log::~kkc_Log()
{
    flush();
}

kkc_Log &kkc_Log::get_instance()
{
    static kkc_Log instance;
    return instance;
}

template <typename T>
kkc_Log &kkc_Log::operator<<(const T &content)
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_ << content;
    return *this;
}

kkc_Log &kkc_Log::operator<<(const char *content)
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_ << content;
    return *this;
}

kkc_Log &kkc_Log::operator<<(std::ostream &(*manip)(std::ostream &))
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_ << manip;
    return *this;
}
kkc_Log &kkc_Log::operator<<(LogLevel level)
{
    append(level); // 复用 append 的逻辑
    return *this;
}

kkc_Log &kkc_Log::operator<<(ErrorCode code)
{
    append(code); // 复用 append 的逻辑
    return *this;
}

void kkc_Log::append(const std::string &content)
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_ << content;
}

void kkc_Log::append(LogLevel level)
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_ << toString(level);
}
void kkc_Log::append(ErrorCode code)
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_ << toString(code);
}

void kkc_Log::print_to_Console()
{
    std::lock_guard<std::mutex> lock(mtx_);
    std::cout << ss_.str();
}

void kkc_Log::print_to_file(const std::string &filename, bool append)
{
    std::lock_guard<std::mutex> lock(mtx_);
    if (filename.empty())
    {
        return; // 不处理空文件名
    }
    std::ofstream file;
    if (append)
    {
        file.open(filename, std::ios::app);
    }
    else
    {
        file.open(filename, std::ios::trunc);
    }
    if (file.is_open())
    {
        file << ss_.str();
        file.close();
    }
}

void kkc_Log::flush()
{
    std::lock_guard<std::mutex> lock(mtx_);
    ss_.str(""); // 清空缓冲区
    ss_.clear(); // 清除错误状态
}

void kkc_Log::set_default_file(const std::string &filename)
{
    std::lock_guard<std::mutex> lock(mtx_);
    filename_ = filename;
}

std::string kkc_Log::get_buffer()
{
    std::lock_guard<std::mutex> lock(mtx_);
    return ss_.str();
}

std::string kkc_Log::toString(LogLevel level)
{
    switch (level)
    {
    case LogLevel::Debug:
        return "[DEBUG] ";
    case LogLevel::Info:
        return "[INFO] ";
    case LogLevel::Warning:
        return "[WARNING] ";
    case LogLevel::Error:
        return "[ERROR] ";
    default:
        return "<unknown log level> ";
    }
}
std::string kkc_Log::toString(ErrorCode code)
{
    switch (code)
    {

    case ErrorCode::Success:
        return "Success";
    case ErrorCode::UnknownError:
        return "Unknown error";
    case ErrorCode::Blank_input:
        return "Blank input";
    case ErrorCode::Wrong_value:
        return "Wrong value range";
    case ErrorCode::Wrong_param_type:
        return "Wrong param type";

    default:
        return "<unknown error code>\n";
    }
}

/**
 * @brief 用于输出“空白输入/错误变量类型”的错误
 * @param log 向该log输出信息
 * @param params_description 对参数类型的说明，如AAA：int
 * @param filePath 文件路径
 * @param line_num 行数，应该由对应文件的行数指示器决定
 * @param raise_err 该值为T时抛出runtime err，终止程序，为F时仅向控制台输出错误信息
 * @return 永远为true，对接err—flag
 */
bool kkc_Log::raise_err_WrongType(std::string params_description, const std::string &filePath, int line_num, bool raise_err)
{
    *this << kkc_Log::LogLevel::Error << kkc_Log::ErrorCode::Blank_input << "/" << kkc_Log::ErrorCode::Wrong_param_type << "--failed to get/input wrong param type in 1 or more params (" << params_description << ") in file: "
          << filePath << " , line: " << line_num << std::endl;
    if (raise_err)
    {
        this->print_to_Console();
        throw std::runtime_error("----\nError in file: " + filePath + " , program terminated");
    }
    return true;
}

/**
 * @brief 用于抛出runtime err，但是在抛出前确保log打印到控制台上
 * @param log 向该log输出信息
 * @param file_type string类型，如“.env”
 * @param filePath 文件路径
 */
void kkc_Log::raise_runtime_err_and_exit(std::string file_type, const std::string &filePath)
{
    this->print_to_Console();
    std::string temp_text = "----\nError in " + file_type + " file: " + filePath + " , program terminated";
    std::cout << temp_text << std::endl;
    throw std::runtime_error(temp_text);
}

// 模板函数的显式实例化（可选，避免链接错误）
template kkc_Log &kkc_Log::operator<<(const int &content);
template kkc_Log &kkc_Log::operator<<(const double &content);
template kkc_Log &kkc_Log::operator<<(const std::string &content);
