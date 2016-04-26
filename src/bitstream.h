/** @file bitstream.h
 *  @brief This file defines the ibitstream and obitstream classes.
 *
 *  This file defines the ibitstream and obitstream classes which are basically
 *  the same as the ordinary istream and ostream classes, but add the
 *  functionality to read and write one bit at a time.
 *
 *  The idea is that you can substitute an ibitstream in place of an
 *  istream and use the same operations (get, fail, >>, etc.) along with added 
 *  member functions of readBit, rewind, and size.
 *
 *  Similarly, the obitstream can be used in place of ostream, and has
 *  same operations (put, fail, <<, etc.) along with additional member 
 *  functions writeBit and size.
 *
 *  There are two subclasses of ibitstream: ifbitstream and istringbitstream,
 *  which are similar to the ifstream and istringstream classes. The
 *  obitstream class similarly has ofbitstream and ostringbitstream as
 *  subclasses.
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef BITSTREAM_H
#define BITSTREAM_H

#include <fstream>
#include <iostream>
#include <sstream>

/** @brief Class: ibitstream
 *
 *  Defines a class for reading files with all the functionality of istream
 *  along with an added member function for reading a single bit and convenience
 *  functions for rewinding the stream back to the beginning and getting the 
 *  stream size.
 *  You will probably not create instances of this class directly. Instead, you
 *  will create ifbitstreams or istringbitstreams to read from files or string
 *  buffers.
 */
class ibitstream: public std::istream {
public:

    /** @brief Constructor: ibitstream
     *
     *  Initializes a new ibitstream that is not attached to any source.
     *  You are unlikely to use this function directly.
     */
    ibitstream(void);

    /** @brief Member function: readBit
     *
     *  Reads a single bit. Raises an error if this ibitstream has not
     *  been properly opened.
     *
     *  @return Returns 0 or 1 depending on the bit value. If the stream is
     *          exhausted, EOF (-1) is returned.
     */
    int readBit(void);

    /** @brief Member function: rewind
     *
     *  Rewinds the ibitstream back to the beginning so that subsequent reads
     *  start again from the beginning. Raises an error if this ibitstream
     *  has not been properly opened.
     *
     *  @return Void.
     */
    void rewind(void);

    /** @brief Member function: size
     *
     *  Returns the size in bytes of the data attached to this stream.
     *  Raises an error if this ibitstream has not been properly opened.
     *
     *  @return Returns the size in bytes of the data attached to this stream.
     */
    long size(void);

    /** @brief Member function: is_open
     *
     *  Returns whether or not this ibitstream is opened. This only has
     *  meaning if the ibitstream is a file stream; otherwise it always
     *  returns true.
     *
     *  @return Returns true if this ibitstream is open and false otherwise.
     */
    virtual bool is_open(void);

private:
    std::streampos lastTell;
    int currByte;
    int pos;
};


/** @brief Class: obitstream
 *
 *  Defines a class for writing files with all the functionality of ostream
 *  along with an added member function for writing a single bit and a
 *  convenience function for getting the stream size.
 *  You are unlikely to instantiate this class directly; instead, instantiate
 *  one of the subclasses.
 */
class obitstream: public std::ostream {
public:
    /** @brief Constructor: obitstream
     *
     *  Initializes a new obitstream that is not attached to any file. Use the
     *  open member function from ofstream to attach the stream to a file.
     */
    obitstream(void);

    /** @brief Member function: writeBit
     *
     *  Writes a single bit to the obitstream.
     *  Raises an error if this obitstream has not been properly opened.
     *
     *  @param bit Bit to write
     *  @return Void.
     */
    void writeBit(int bit);

    /** @brief Member function: size
     *
     *  Returns the size in bytes of the file attached to this stream.
     *  Raises an error if this obitstream has not been properly opened.
     * 
     *  @return Returns the size in bytes of the file attached to this stream.
     */
    long size(void);

    /** @brief Member function: is_open
     *
     *  Returns whether or not this obitstream is opened. This only has
     *  meaning if the obitstream is a file stream; otherwise it always
     *  returns true.
     *
     *  @return Returns true if this obitstream is open and false otherwise.
     */
    virtual bool is_open(void);

private:
    std::streampos lastTell;
    int currByte;
    int pos;
};


/** @brief Class: ifbitstream
 *
 *  A class for reading files in all of the usual ways, plus bit-by-bit.
 *  You can treat this class like a normal ifstream, except that there is
 *  extra support for bit-level operations.
 */
class ifbitstream: public ibitstream {
public:
    /** @brief Constructor: ifbitstream
     *
     *  Constructs a new ifbitstream not attached to any file. You can open a 
     *  file for reading using the .open() member function.
     */
    ifbitstream(void);

    /** @brief Constructor: ifbitstream
     *
     *  Constructs a new ifbitstream that reads the specified file, if it 
     *  exists. If not, the stream enters an error state.
     *
     *  @param filename File to read from
     */
    ifbitstream(const char *filename);
    ifbitstream(const std::string &filename);

    /** @brief Member function: open
     *
     *  Opens the specified file for reading. If an error occurs, the stream
     *  enters a failure state, which can be detected by calling .fail().
     *
     *  @param filename File to open
     *  @return Void.
     */
    void open(const char *filename);
    void open(const std::string &filename);

    /** @brief Member function: is_open
     *
     *  Returns whether or not this ifbitstream is connected to a file for
     *  reading.
     *
     *  @return Returns true if ifbitstream is open and false otherwise.
     */
    bool is_open(void);

    /** @brief Member function: close
     *
     *  Closes the currently-opened file, if the stream is open. If the stream 
     *  is not open, puts the stream into a fail state.
     *
     *  @return Void.
     */
    void close(void);

private:
    /* The actual file buffer which does reading and writing. */
    std::filebuf fileBuffer;
};


/** @brief Class: ofbitstream
 *
 *  A class for writing files in all of the usual ways, plus bit-by-bit.
 *  You can treat this class like a normal ofstream, except that there is
 *  extra support for bit-level operations.
 */
class ofbitstream: public obitstream {
public:
    /** @brief Constructor: ofbitstream
     *
     *  Constructs a new ofbitstream not attached to any file. You can open a 
     *  file for writing using the .open() member functions.
     */
    ofbitstream(void);

    /** @brief Constructor: ofbitstream
     *
     *  Constructs a new ofbitstream that writes the specified file, if it 
     *  exists. If not, the stream enters an error state. Read the documentation
     *  on "open" for more details.
     *
     *  @param filename File to write to
     */
    ofbitstream(const char *filename);
    ofbitstream(const std::string &filename);

    /** @brief Member function: open
     *
     *  Opens the specified file for writing. If an error occurs, the stream
     *  enters a failure state, which can be detected by calling .fail(). If an
     *  invalid filename is specified (for example, a source file), reports an
     *  error.
     *
     *  @param filename File to open
     *  @return Void.
     */
    void open(const char *filename);
    void open(const std::string &filename);

    /** @brief Member function: is_open
     *
     *  Returns whether or not this ofbitstream is connected to a file for
     *  reading.
     *
     *  @return Returns true if this ofbitstream is open and false otherwise.
     */
    bool is_open(void);

    /** @brief Member function: close
     *
     *  Closes the currently-opened file, if the stream is open. If the stream
     *  is not open, puts the stream into a fail state.
     *
     *  @return Void.
     */
    void close(void);

private:
    /* The actual file buffer which does reading and writing. */
    std::filebuf fileBuffer;
};


/** @brief Class: istringbitstream
 *
 *  A variant on C++'s istringstream class, which acts as a stream that reads
 *  its data from a string. This is mostly used by testing code.
 */
class istringbitstream: public ibitstream {
public:
    /** @brief Constructor: istringbitstream
     *
     *  Constructs an istringbitstream reading the specified string.
     *
     *  @param s Empty string
     */
    istringbitstream(const std::string &s = "");

    /** @brief Member Function: str
     *
     *  Sets the underlying string of the istringbitstream.
     *
     *  @return Void.
     */
    void str(const std::string &s);
private:
    /// The actual string buffer that does character storage.
    std::stringbuf stringBuffer;
};


/** @brief Class: ostringbitstream
 *
 *  A variant on C++'s ostringstream class, which acts as a stream that writes
 *  its data to a string. This is mostly used by testing code.
 */
class ostringbitstream: public obitstream {
public:
    /** @brief Constructor: ostringbitstream
     *
     *  Constructs an ostringbitstream.
     */
    ostringbitstream(void);

    /** @brief Member function: str
     *
     *  Retrieves the underlying string of the istringbitstream.
     *
     *  @return Returns the string associated to this istringbitstream.
     */
    std::string str(void);

private:
    /// The actual string buffer that does character storage.
    std::stringbuf stringBuffer;
};

#endif // BITSTREAM_H

