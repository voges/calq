#ifndef CALQAPP_ERROR_REPORTER_H_
#define CALQAPP_ERROR_REPORTER_H_

// -----------------------------------------------------------------------------

#include <exception>
#include <iostream>
#include <string>
#include <utility>

// -----------------------------------------------------------------------------

#include "calq/exceptions.h"


// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

inline void throwErrorException(const std::string& msg){
    std::cout.flush();
    throw calq::ErrorException(msg);
}

// -----------------------------------------------------------------------------

class ErrorExceptionReporter
{
    // -------------------------------------------------------------------------

 public:
    ErrorExceptionReporter(std::string file,
                           std::string function,
                           const int& line
    )
            : file_(std::move(file)),
            function_(std::move(function)),
            line_(line){
    }

    // -------------------------------------------------------------------------

    void operator()(const std::string& msg){
        // std::cerr << file << ":" << function << ":" << line << ": ";
        std::string tmp = file_.substr(file_.find_last_of("/\\") + 1) +
                          ":" +
                          function_ +
                          ":" +
                          std::to_string(line_) +
                          ": " +
                          msg;
        // Can use the original name here, as it is still defined
        throwErrorException(tmp);
    }

    // -------------------------------------------------------------------------

 private:
    std::string file_;
    std::string function_;
    int line_;
};

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------

// Remove the symbol for the function, then define a new version that instead
// creates a stack temporary instance of ErrorExceptionReporter initialized
// with the caller.
#undef throwErrorException
#define throwErrorException calqapp::ErrorExceptionReporter(__FILE__, \
                                                            __FUNCTION__, \
                                                            __LINE__ \
)

// -----------------------------------------------------------------------------

#endif  // CALQAPP_ERROR_REPORTER_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
