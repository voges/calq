/** @file rle.h
 *  @brief This files contains the interface to the RLE codec.
 * 
 *  Note it is up to the calling code to ensure that no overruns on input and
 *  output buffers occur!
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef RLE_H
#define RLE_H

#include <string>

/** @brief Function: rle_encode
 *
 *  Encodes a byte array with run-length encoding.
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
 *  @param in Byte array to encode
 *  @param in_sz Size in bytes of in
 *  @param out_sz Pointer to variable the output size is written to
 *  @param N Number of symbols in input string
 *  @param offset Offset of input string
 *  @return Returns a pointer to the byte buffer containing the encoded data,
 *          has to be freed by the caller (!)
 */
unsigned char * rle_encode(unsigned char       *in,
                           const size_t        in_sz,
                           size_t              *out_sz,
                           const unsigned char N,
                           const unsigned char offset = (unsigned char)'0');

/** @brief Function: rle_decode
 *
 *  Decodes byte stream encoded with rle_encode.
 *
 *  @param in Byte array to decode
 *  @param in_sz Size in bytes of in
 *  @param out_sz Pointer to variable the output size is written to
 *  @param N Number of symbols in data (must equal N used for encoding)
 *  @param offset Offset of data (must equal offset used for encoding)
 *  @return Returns a pointer to the byte buffer containing the decoded data,
 *          has to be freed by the caller (!)
 */
unsigned char * rle_decode(unsigned char       *in,
                           const size_t        in_sz,
                           size_t              *out_sz,
                           const unsigned char N,
                           const unsigned char offset = (unsigned char)'0');

#endif // RLE_H

