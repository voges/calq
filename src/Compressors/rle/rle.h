/** @file rle.h
 *  @brief This files contains the interface to the RLE codec.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef RLE_H
#define RLE_H

#include <string>

/** @brief Function: rle_encode
 *
 *  Encodes a string message with run-length encoding.
 *  Example 1:
 *    in = "0112223334443"
 *    => N = 5
 *    => offset = '0' = 48
 *    ===> out = 0x 00 01 01 08 02 08 03 08 04 03
 *  Example 2:
 *    in = 0x 00 01 02 00 00 00 02
 *    => N = 3
 *    => offset = 0
 *    ===> out = 0x 00 01 02 06 00 02
 *
 *  @param in Input string (i.e. byte array) to encode
 *  @param out Output string (i.e. byte array) containing run-length encoding
 *  @param N Number of symbols in input string
 *  @param offset Offset of input string
 *  @return Void
 */
void rle_encode(const std::string &in,
                std::string &out,
                const unsigned int &N,
                const unsigned int &offset = (unsigned int)'0');

/** @brief Function: rle_decode
 *
 *  Decodes byte stream encoded with rle_encode.
 *
 *  @param in Input string (i.e. byte array) to decode
 *  @param out Output string (i.e. byte array) containing decoded data
 *  @param N Number of symbols in data (must equal N used for encoding)
 *  @param offset Offset of data (must equal offset used for encoding)
 *  @return Void
 */
void rle_decode(const std::string &in,
                std::string &out,
                const unsigned int &N,
                const unsigned int &offset = (unsigned int)'0');

#endif // RLE_H

