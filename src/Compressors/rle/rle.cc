/** @file rle.cc
 *  @brief This files contains the implementation of the RLE codec.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "rle.h"
#include "Common/debug.h"
#include "Common/Exceptions.h"
#include <iostream>

//
// Decoder specification:
// ----------------------
// A byte can have 256 values: from 0x00 = 0 up to 0xFF = 255. By definition
// the smallest run-length is 3. Therefore, with one byte, we can express
// run-lengths from 3 (= 0x00) up to 258 (= 0xFF).
// Assume, we only want to encode N symbols 0...(N-1). We assign the N byte
// values 0x00 - 0x** = (N-1) to our symbols. Then, we have (256-N) byte values
// left to encode run-lengths from 3 up to 258-N.
//
// Example for N = 4 symbols 0,1,2,3. We can encode 256 - 4 = 252 run-lengths
// ranging from 3 up to 254:
//   Symbol 0 is encoded as 0x00
//   Symbol 1 is encoded as 0x01
//   Symbol 2 is encoded as 0x02
//   Symbol 3 is encoded as 0x03
//   Run-length 3 is encoded as 0x04
//   Run-length 4 is encoded as 0x05
//   ...
//   Run-length 254 is encoded as 0xFF
//

static const unsigned int RLE_MIN_RUN_LENGTH = 3;
static const unsigned int RLE_MAX_RUN_LENGTH = 258;

void rle_encode(const std::string &in,
                std::string &out,
                const unsigned int &N,
                const unsigned int &offset)
{
    // Less than 2 symbols do not make sense, the same holds for more than
    // 255 symbols, because we need at least one byte to represent a run-length
    // of 3.
    if (N < 2 || N > 255) {
        throwErrorException("N out of range");
    }

    // We have at least 2 symbols. The highest possible values for these are
    // 0xFE = 254 and 0xFF = 255, respectively. Thus, an offset greater than
    // 254 cannot make sense.
    if (offset > 254) {
        throwErrorException("offset out of range");
    }

    const unsigned int maxRunLength = RLE_MAX_RUN_LENGTH - N;
    out = "";
    char prevChar = in[0] - offset;
    if ((unsigned int)prevChar > N-1) {
        throwErrorException("Symbol in stream out of range");
    }
    unsigned int n = 1;

    for (size_t i = 1; i < in.length(); i++) {
        char currChar = in[i] - offset;
        if ((unsigned int)currChar > N-1) {
            throwErrorException("Symbol in stream out of range");
        }
        if ((currChar != prevChar) || (n == maxRunLength)) {
            if (n < RLE_MIN_RUN_LENGTH) {
                for (size_t j = 0; j < n; j++) {
                    out += prevChar;
                }
            } else {
                out += (char)(n + N);
                out += prevChar;
            }
            n = 0;
        }
        prevChar = currChar;
        n++;
    }

    // Write out last run
    if (n < RLE_MIN_RUN_LENGTH) {
        for (size_t j = 0; j < n; j++) {
            out += (char)prevChar;
        }
    } else {
        out += (char)(n + N);
        out += prevChar;
    }

}

void rle_decode(const std::string &in,
                std::string &out,
                const unsigned int &N,
                const unsigned int &offset)
{
    // Less than 2 symbols do not make sense, the same holds for more than
    // 255 symbols, because we need at least one byte to represent a run-length
    // of 3.
    if (N < 2 || N > 255) {
        throwErrorException("N out of range");
    }

    // We have at least 2 symbols. The highest possible values for these are
    // 0xFE = 254 and 0xFF = 255, respectively. Thus, an offset greater than
    // 254 cannot make sense.
    if (offset > 254) {
        throwErrorException("offset out of range");
    }

    out = "";
    size_t i = 0;

    while (i < in.length()) {
        char currChar = in[i++];

        // Check, if currChar contains a run-length
        if ((unsigned int)currChar > N-1) {
            char runSymbol = in[i++] + (char)offset;
            auto runLength = currChar - N;
            for (size_t j = 0; j < runLength; j++) {
                out += runSymbol;
            }
        } else {
            out += currChar + (char)offset;
        }
    }
}

