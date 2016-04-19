#ifndef CQ_QUALENCODER_HPP
#define CQ_QUALENCODER_HPP

#include "sam/samrec.h"
#include <stdio.h>
#include <vector>

namespace cq {

class QualEncoder {
  private:
    size_t recordCnt; // number of records processed in the current block
    size_t posMin;    // smallest mapping position encountered in block
    size_t posMax;    // largest mapping position encountered in block
    std::vector<size_t> depths;

  public:
    QualEncoder(void);
    ~QualEncoder(void);

    /** @brief Extracts the quality scores (plus some other information) from a
     *         SAM record and encodes them. The encoded quality score are kept 
     *         in a buffer which can be written to a file.
     *  @param samrec The SAM record to be encoded
     */
    void encodeRecord(const samrec_t * const samrec);

    /** @brief Writes encoded quality score to the given stream.
     *  @param fp Pointer to the output file
     *  @return Returns the number of written bytes
     */
    size_t finishBlock(FILE *fp);
};

}

#endif // CQ_QUALENCODER_HPP

