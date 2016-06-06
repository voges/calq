/** @file ErrorException.h
 *  @brief This files contains the definitions of the ErrorException class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef ERROREXCEPTION_H
#define ERROREXCEPTION_H

#include <sstream>
#include <stdexcept>
#include <string>

/** @brief Class: ErrorException
 *
 *  This class is designed to make it a little easier to throw informative
 *  exceptions. It's a little lame, but I do like to be able to write code
 *  like this errors:
 *
 *  throw ErrorException() << "Game over," << mHealth << " health points!";
 */
class ErrorException : public std::exception {
public :
    /** @brief Constructor: ErrorException
     *
     *  The standard constructor.
     */
    ErrorException() {};

    /** @brief Copy constructor: ErrorException
     *
     *  Need a copy constructor, because the runtime of the compiler is
     *  allowed to insist on it when it throws an object of this type,
     *  even if it doesn't actually make a copy. When I make a copy, I
     *  need to capture whatever is in the stringstream object. Note:
     *  in many cases, attempting to copy an iostream object leads to
     *  errors, so the copy constructor here constructs a brand new
     *  mstream object.
     *
     *  @param that The ErrorException instance to be copied during
     *         construction.
     */
    ErrorException(const ErrorException &that)
    {
        mWhat += that.mStream.str();
    }

    /** @brief Destructor: ErrorException
     *
     *  Virtual destructor not really needed; here it is anyway
     */
    virtual ~ErrorException() throw() {};

    /** @brief Member function: what
     *
     *  When I finally get this object to an exception handler,
     *  (hopefully catching by reference) I want to display the error
     *  message that I've inserted. To do that, I just capture
     *  whatever is in the mWhat string object concatenated
     *  with anything that might be in the stringstring mStream object,
     *  and return it. (Odds are that only one of them will contain
     *  anything, depending on whether or not the copy constructor
     *  was called.
     *
     *  @return Returns the exception message as C-String (a.k.a. const char *)
     */
    virtual const char * what() const throw()
    {
        if (mStream.str().size()) {
            mWhat += mStream.str();
            mStream.str("");
        }
        return mWhat.c_str();
    }

    /** @brief Template member function: operator<<
     *
     *  The template function used to create insertion operators for all
     *  of the various types of objects one might insert into this guy.
     *
     *  @param t Argument to be inserted into the message stream
     *  @return Void
     */
    template<typename T>
    ErrorException & operator<<(const T &t)
    {
        mStream << t;
        return *this;
    }

private:
    mutable std::stringstream mStream;
    mutable std::string mWhat;
};

#endif // ERROREXCEPTION_H

