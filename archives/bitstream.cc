/** @file bitstream.cc
 *  @brief Implementation of i...bitstream and o...bitstream  classes.
 *
 *  This file contains the implementation of i...bitstream and o...bitstream
 *  classes. These classes are patterned after (and, in fact, inherit from) the
 *  standard ifstream and ofstream classes. Please see bitstream.h for
 *  information about how a client properly uses these classes.
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "bitstream.h"
#include "Common/Exceptions.h"

static const int NUM_BITS_IN_BYTE = 8;

static inline int getNthBit(const int n, const int fromByte) 
{
    return ((fromByte & (1 << n)) != 0);
}

static inline void setNthBit(const int n, int &inByte) 
{
    inByte |= (1 << n);
}

ibitstream::ibitstream(void)
    : std::istream(NULL)
    , lastTell(0)
    , currByte(0)
    , pos(-1)
{
    ///
    /// Each ibitstream tracks 3 integers as private data.
    /// - "lastTell" is streampos of the last byte that was read (this is used
    ///   to detect when other non-readBit activity has changed the tell).
    /// - "currByte" contains contents of byte currently being read.
    /// - "pos" is the bit position within currByte that is next to read.
    /// We set initial states for lastTell and currByte to 0, then pos is set at 8 
    /// so that next readBit will trigger a fresh read.
    ///
}

int ibitstream::readBit(void)
{
    ///
    /// If bits remain in currByte, retrieve next bit and increment pos.
    /// Else if end of currByte (or some other read happened), then read next
    /// byte and start reading from bit position 0 of that byte.
    /// If read byte from file at EOF, return EOF.
    ///
    
    if (!is_open()) {
        throwErrorException("Cannot read a bit from a stream that is not open");
    }
    
    // if just finished bits from currByte or if data read from stream after
    // last readBit
    if (lastTell != tellg() || pos == -1) {
        // read next single byte from file
        if ((currByte = get()) == EOF) {
            return EOF;
        }
        pos = 7; // start reading from first bit of new byte
        lastTell = tellg();
    }
    
    int result = getNthBit(pos, currByte);
    pos--; // advance bit position for next call to readBit
    return result;
}

void ibitstream::reset(){
    currByte=0;
    pos=-1;
}

int ibitstream::readByte(BYTE &byte)
{
    if (!is_open()) {
        throwErrorException("Cannot read a byte from a stream that is not open");
    }

    if ((byte = (BYTE)get()) == EOF) {
        return EOF;
    }

//    fprintf(stderr, "read 0x%x\n", byte);
    return 0;
}

int ibitstream::readUint64(uint64_t &x)
{
    if (!is_open()) {
        throwErrorException("Cannot read a uint64_t from a stream that is not open");
    }

    x = 0x0;
    
    BYTE hhhh = 0x00;
    BYTE hhh = 0x00;
    BYTE hh = 0x00;
    BYTE h = 0x00;
    BYTE l = 0x00;
    BYTE ll = 0x00;
    BYTE lll = 0x00;
    BYTE llll = 0x00;

    if (readByte(hhhh) == EOF) { return EOF; }
    if (readByte(hhh) == EOF) { return EOF; }
    if (readByte(hh) == EOF) { return EOF; }
    if (readByte(h) == EOF) { return EOF; }
    if (readByte(l) == EOF) { return EOF; }
    if (readByte(ll) == EOF) { return EOF; }
    if (readByte(lll) == EOF) { return EOF; }
    if (readByte(llll) == EOF) { return EOF; }

    x ^= (uint64_t)hhhh << 56;
    x ^= (uint64_t)hhh  << 48;
    x ^= (uint64_t)hh   << 40;
    x ^= (uint64_t)h    << 32;
    x ^= (uint64_t)l    << 24;
    x ^= (uint64_t)ll   << 16;
    x ^= (uint64_t)lll  <<  8;
    x ^= (uint64_t)llll;

    return 0;
}

int ibitstream::readBuffer(BYTE &buf, const size_t n)
{
    if (!is_open()) {
        throwErrorException("Cannot read a byte buffer from a stream that is not open");
    }

    read((char*)&buf, n);

    if (eof()) {
        return EOF;
    }

    if (fail()) {
        throwErrorException("Error on stream while reading buffer");
    }

    return 0;
}

void ibitstream::rewind(void)
{
    ///
    /// Simply seeks back to beginning of file, so reading begins again
    /// from start.
    ///

    if (!is_open()) {
        throwErrorException("Cannot rewind a stream that is not open");
    }

    clear();
    seekg(0, std::ios::beg);
}

long ibitstream::size(void) 
{
    ///
    /// Seek to file end and use tell to retrieve position.
    /// In order to not disrupt reading, we also record current streampos and
    /// re-seek to there before returning.
    ///

    if (!is_open()) {
        throwErrorException("Cannot get size of stream which is not open");
    }

    clear();                  // clear any error state
    streampos curr = tellg(); // save current streampos
    seekg(0, std::ios::end);  // seek to end
    streampos end = tellg();  // get offset
    seekg(curr);              // seek back to original pos

    return long(end);
}

bool ibitstream::is_open(void) 
{
    ///
    /// Default implementation of is_open has the stream always open. Subclasses
    /// can customize this if they'd like.
    ///

    return true;
}

obitstream::obitstream()
    : std::ostream(NULL)
    , lastTell(0)
    , currByte(0)
    , pos(-1)
{
    ///
    /// Each obitstream tracks 3 integers as private data.
    /// - "lastTell" is streampos of the last byte that was written (this is 
    ///   used to detect when other non-writeBit activity has changed the tell).
    /// - "currByte" contains contents of byte currently being written.
    /// - "pos" is the bit position within currByte that is next to write.
    /// We set initial state for lastTell and currByte to 0, then pos is set at
    /// 8 so that next writeBit will start a new byte.
    ///
}

void obitstream::writeBit(const int bit)
{
    ///
    /// If bits remain to be written in currByte, add bit into byte and increment
    /// pos.
    /// Else if end of currByte (or some other write happened), then start a fresh
    /// byte at position 0.
    /// We write the byte out for each bit (backing up to overwrite as needed),
    /// rather than waiting for 8 bits. This is because the client might make
    /// 3 writeBit calls and then start using << so we can't wait til full-byte
    /// boundary to flush any partial-byte bits.
    ///
    
    if (bit != 0 && bit != 1) {
        throwErrorException("Must pass an integer argument of 0 or 1");
    }
    
    if (!is_open()) {
        throwErrorException("Stream is not open");
    }
    
    // if just filled currByte or if data written to stream after last
    // writeBit
    if (lastTell != tellp() || pos == -1) {
        currByte = 0; // zero out byte for next writes
        pos = 7;      // start writing to first bit of new byte
    }
    
    if (bit) {
        // only need to change if bit needs to be 1 (byte starts already zeroed)
        setNthBit(pos, currByte);
    }
    
    // only write if first bit in byte or changing 0 to 1
    if (pos == 7 || bit) {
        if (pos != 7) {
            seekp(-1, std::ios::cur); // back up to overwrite if pos > 0
        }
        put(currByte);
    }
    
    
    pos--; // advance to next bit position for next write
    lastTell = tellp();
    
}

int obitstream::reset(){
    if(pos != -1){
        pos = -1;
        return 1;
    }
    currByte = 0;
    return 0;
    }

void obitstream::writeByte(const BYTE byte)
{
    if (!is_open()) {
        throwErrorException("Stream is not open");
    }
    put((char)byte);
    fprintf(stderr, "write 0x%x\n", byte);
}

void obitstream::writeUint64(const uint64_t x)
{
    if (!is_open()) {
        throwErrorException("Stream is not open");
    }

    writeByte((x >> 56) & 0xFF);
    writeByte((x >> 48) & 0xFF);
    writeByte((x >> 40) & 0xFF);
    writeByte((x >> 32) & 0xFF);
    writeByte((x >> 24) & 0xFF);
    writeByte((x >> 16) & 0xFF);
    writeByte((x >>  8) & 0xFF);
    writeByte((x      ) & 0xFF);
}

void obitstream::writeBuffer(const BYTE &buf, const size_t n)
{
    if (!is_open()) {
        throwErrorException("Stream is not open");
    }

    write((char*)&buf, n);

    if (fail()) {
        throwErrorException("Stream error while writing buffer");
    }
}

long obitstream::size(void) 
{
    ///
    /// Seek to file end and use tell to retrieve position.
    /// In order to not disrupt writing, we also record cur streampos and
    /// re-seek to there before returning.
    ///

    if (!is_open()) {
        throwErrorException("Stream is not open");
    }

    clear();                  // clear any error state
    streampos curr = tellp(); // save current streampos
    seekp(0, std::ios::end);  // seek to end
    streampos end = tellp();  // get offset
    seekp(curr);              // seek back to original pos

    return long(end);
}

bool obitstream::is_open() 
{
    ///
    /// Default implementation of is_open has the stream always open. Subclasses
    /// can customize this if they'd like.
    ///

    return true;
}

ifbitstream::ifbitstream(void) 
{
    ///
    /// Wire up the stream class so that it knows to read data from disk.
    ///

    init(&fileBuffer);
}

ifbitstream::ifbitstream(const char *filename)
{
    ///
    /// Wire up the stream class so that it knows to read data from disk, then 
    /// open the given file.
    ///

    init(&fileBuffer);
    open(filename);
}

ifbitstream::ifbitstream(const std::string &filename)
{
    ///
    /// Wire up the stream class so that it knows to read data from disk, then 
    /// open the given file.
    ///

    init(&fileBuffer);
    open(filename);
}

void ifbitstream::open(const char *filename)
{
    ///
    /// Attempt to open the specified file, failing if unable to do so.
    ///

    if (!fileBuffer.open(filename, std::ios::in | std::ios::binary)) {
        setstate(std::ios::failbit);
    }
}

void ifbitstream::open(const std::string &filename)
{
    open(filename.c_str());
}

bool ifbitstream::is_open(void)
{
    return fileBuffer.is_open();
}

void ifbitstream::close(void)
{
    if (!fileBuffer.close()) {
        setstate(std::ios::failbit);
    }
}

ofbitstream::ofbitstream(void)
{
    ///
    /// Wire up the stream class so that it knows to write data to disk.
    ///

    init(&fileBuffer);
}

ofbitstream::ofbitstream(const char *filename)
{
    ///
    /// Wire up the stream class so that it knows to write data to disk, then
    /// open the given file.
    ///

    init(&fileBuffer);
    open(filename);
}

ofbitstream::ofbitstream(const std::string &filename)
{
    init(&fileBuffer);
    open(filename);
}

void ofbitstream::open(const char *filename)
{
    ///
    /// Attempt to open the specified file.
    ///

    if (!fileBuffer.open(filename, std::ios::out | std::ios::binary)) {
        setstate(std::ios::failbit);
    }
}
void ofbitstream::open(const std::string &filename)
{
    open(filename.c_str());
}

bool ofbitstream::is_open(void) 
{
    return fileBuffer.is_open();
}

void ofbitstream::close(void) 
{
    if (!fileBuffer.close()) {
        setstate(std::ios::failbit);
    }
}

istringbitstream::istringbitstream(const std::string &s)
{
    ///
    /// Sets the stream to use the string buffer, then sets the initial string
    /// to the specified value.
    ///

    init(&stringBuffer);
    stringBuffer.str(s);
}

void istringbitstream::str(const std::string &s)
{
    ///
    /// Sets the underlying string in the buffer to the specified string.
    ///

    stringBuffer.str(s);
}

ostringbitstream::ostringbitstream(void) {
    init(&stringBuffer); /// Set the stream to use the string buffer
}

std::string ostringbitstream::str(void)
{
    return stringBuffer.str();
}

