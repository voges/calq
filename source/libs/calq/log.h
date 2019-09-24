#ifndef CALQ_LOG_H_
#define CALQ_LOG_H_

#include <chrono>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>

#define LOG_FATAL(x) { calq::Log &log = calq::Log::instance(); log.out(x, calq::Log::Level::FATAL, std::cerr); }
#define LOG_ERROR(x) { calq::Log &log = calq::Log::instance(); log.out(x, calq::Log::Level::ERROR, std::cerr); }
#define LOG_WARN(x) { calq::Log &log = calq::Log::instance(); log.out(x, calq::Log::Level::WARN, std::cerr); }
#define LOG_INFO(x) { calq::Log &log = calq::Log::instance(); log.out(x, calq::Log::Level::INFO, std::cout); }
#define LOG_DEBUG(x) { calq::Log &log = calq::Log::instance(); log.out(x, calq::Log::Level::DEBUG, std::cout); }
#define LOG_TRACE(x) { calq::Log &log = calq::Log::instance(); log.out(x, calq::Log::Level::TRACE, std::cout); }

namespace calq {

static std::string currentIso8601TimeUtc() {
    auto now = std::chrono::system_clock::now();
    auto itt = std::chrono::system_clock::to_time_t(now);
    std::ostringstream ss;
    ss << std::put_time(gmtime(&itt), "%FT%TZ");  // ISO 8601 format: 2007-04-05T14:30:21Z
    return ss.str();
}

class Log {
 public:
     enum class Level {
         FATAL = 0,
         ERROR = 1,
         WARN = 2,
         INFO = 3,
         DEBUG = 4,
         TRACE = 5
     };

    static Log &instance() {
        static Log log;
        return log;
    }

    static void setLevel(const Level level = Level::INFO) {
        std::lock_guard<std::mutex> lock(instance().mtx_);
        instance().level_ = level;
    }

    void out(const std::string& msg, const Level level = Level::INFO, std::ostream& os = std::cout) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (static_cast<int>(level) <= static_cast<int>(level_)) {
            std::ostringstream ss;
            ss << "[" << currentIso8601TimeUtc() << "] [" << levelStr(level) << "] "  << msg << std::endl;
            os << ss.str();
            os.flush();
        }
    }

    static std::string levelStr(const Level level) {
        switch (level) {
            case Level::FATAL: return "fatal";
            case Level::ERROR: return "error";
            case Level::WARN: return "warn";
            case Level::INFO: return "info";
            case Level::DEBUG: return "debug";
            case Level::TRACE: return "trace";
            default: return "invalid";
        }
    }

    Level level_;
    std::mutex mtx_;

 private:
    Log() { level_ = Level::INFO; }
    Log(const Log&) = delete;
    Log& operator=(const Log&) = delete;
    Log(Log &&) = delete;
    Log &operator=(Log &&) = delete;
    ~Log() = default;
};

inline bool isValidLogLevel(const int l) {
    if (l < static_cast<int>(Log::Level::FATAL) || l > static_cast<int>(Log::Level::TRACE)) {
        return false;
    }
    return true;
}

}  // namespace calq

#endif  // CALQ_LOG_H_
