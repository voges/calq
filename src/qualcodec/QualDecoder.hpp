/** @file QualDecoder.hpp
 *  @brief QualDecoder class
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CQ_QUALDECODER_HPP
#define CQ_QUALDECODER_HPP

#include "sam/samrec.h"
#include "misc/str.h"
#include <stdio.h>
#include <vector>

namespace cq {

class QualDecoder {
  private:

  public:
    QualDecoder(void);
    ~QualDecoder(void);

    /** @brief Decodes a block of encoded quality score from the given stream and
     *         writes the decoded quality scores to qual.
     *  @param fp The file pointer to read from; has to be positioned at the
     *         beginning of an encoded quality score block
     *  @param qual An preallocated array of empty strings
     *  @param n Number of quality score vectors to decode
     */
    void decodeBlock(FILE *fp, str_t* qual[], size_t n);
};

}

#endif // CQ_QUALDECODER_HPP

