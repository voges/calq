/** @file Exceptions.h
 *  @brief This files contains the definitions of some custom exception 
 *         classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <exception>
#include <iostream>
#include <string>

/** @brief Class: Exception
 *
 *  Specific exception classes inherit from this class.
 */
class Exception : public std::exception {
public:
    Exception(const std::string &msg);
    virtual ~Exception(void) throw();
    virtual std::string getMessage(void) const;
    virtual const char * what(void) const throw();

protected:
    std::string msg;
};

/** @brief Class: UserException
 *
 *  This exception is thrown by calls to the throwUserException function.
 *  If a UserException is thrown at any point in the range of the try block
 *  (including in functions called from that code), control will jump
 *  immediately to the error handler.
 */
class UserException : public Exception {
public:
    UserException(const std::string &msg): Exception(msg) {}
};

/** @brief Global function: throwUserException
 *
 *  Signals an error condition in a program by throwing a
 *  UserException with the specified message.
 *
 *  @param msg Error message
 *  @return Void
 */
inline void throwUserException(const std::string &msg)
{
    throw UserException(msg);
}

/** @brief Class: ErrorException
 *
 *  This exception is thrown by calls to the >throwErrorException function.
 *  If an ErrorException is thrown at any point in the range of the try block
 *  (including in functions called from that code), control will jump
 *  immediately to the error handler.
 */
class ErrorException : public Exception {
public:
    ErrorException(const std::string &msg): Exception(msg) {}
};

/** @brief Global function: throwErrorException
 *
 *  Signals an error condition in a program by throwing an ErrorException
 *  with the specified message.
 *
 *  @param msg Error message
 *  @return Void
 */
inline void throwErrorException(const std::string &msg)
{
    throw ErrorException(msg);
}

/** @brief Class: ErrorExceptionReporter
 *
 *  Reporter class catching the caller of an ErrorException.
 */
class ErrorExceptionReporter {
public:
    ErrorExceptionReporter(const std::string &function, 
                           const std::string &file,
                           const int line)
        : function(function)
        , file(file)
        , line(line)
    {}

    void operator()(const std::string &msg)
    {
        std::cerr << file << ":" << function << ":" << line << ": ";
        // can use the original name here, as it is still defined
        throwErrorException(msg);
    }
private:
    std::string function;
    std::string file;
    int line;
};

// Remove the symbol for the function, then define a new version that instead
// creates a stack temporary instance of ErrorExceptionReporter initialized
// with the caller.
#undef throwErrorException
#define throwErrorException ErrorExceptionReporter(__FUNCTION__, __FILE__, __LINE__)

#endif // EXCEPTIONS_H

