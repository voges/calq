/** @file QualEncoder.h
 *  @brief This file contains the definition of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_QUALENCODER_H_
#define CALQ_QUALCODEC_QUALENCODER_H_

#include <chrono>
#include <deque>
#include <map>
#include <string>

#include "config.h"
#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMPileupDeque.h"
#include "IO/SAM/SAMRecord.h"
#include "QualCodec/Genotyper.h"
#include "QualCodec/Quantizers/Quantizer.h"

namespace calq {

class QualEncoder {
public:
    static const int QUANTIZER_STEPS_MIN = 2;
    static const int QUANTIZER_STEPS_MAX = 8;
    static const int NR_QUANTIZERS = QUANTIZER_STEPS_MAX-QUANTIZER_STEPS_MIN+1;
    static const int QUANTIZER_IDX_MIN = 0;
    static const int QUANTIZER_IDX_MAX = NR_QUANTIZERS-1;

    explicit QualEncoder(const int &polyploidy,
                         const int &qualityValueMax,
                         const int &qualityValueMin,
                         const int &qualityValueOffset,
                         const std::map<int, Quantizer> &quantizers);
    ~QualEncoder(void);

    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    void finishBlock(void);
    size_t writeBlock(CQFile *cqFile);
#if MPEG_CE5_DESCRIPTOR_STREAMS_OUTPUT
    void writeMPEGBlock(File *mpegQVCIFile, File *mpegQVIFile);
#endif
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
    void writeContextAdaptiveMPEGBlock(CQFile *qvciFile, CQFile *qviFile);
#endif

    size_t compressedMappedQualSize(void) const;
    size_t compressedUnmappedQualSize(void) const;
    size_t compressedQualSize(void) const;
    size_t nrMappedRecords(void) const;
    size_t nrUnmappedRecords(void) const;
    size_t nrRecords(void) const;
    size_t uncompressedMappedQualSize(void) const;
    size_t uncompressedUnmappedQualSize(void) const;
    size_t uncompressedQualSize(void) const;
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
    size_t compressedQviSize(void) const;
    size_t compressedQvciSize(void) const;
    size_t uncompressedQvciSize(void) const;
    size_t uncompressedQviSize(void) const;
#endif

private:
    void encodeMappedQual(const SAMRecord &samRecord);
    void encodeUnmappedQual(const std::string &qual);

    // Timekeeping
    std::chrono::time_point<std::chrono::steady_clock> startTime_;
    std::chrono::time_point<std::chrono::steady_clock> stopTime_;

    // Sizes & counters
    size_t compressedMappedQualSize_;
    size_t compressedUnmappedQualSize_;
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
    size_t uncompressedMappedQualSize_;
    size_t uncompressedUnmappedQualSize_;
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
    size_t compressedQviSize_;
    size_t compressedQvciSize_;
    size_t uncompressedQvciSize_;
    size_t uncompressedQviSize_;
#endif

    int qualityValueOffset_;

    // 0-based position offset of this block
    uint32_t posOffset_;

    // Buffers
    std::string unmappedQualityValues_;
    std::deque<int> mappedQuantizerIndices_;
    std::deque<int> mappedQualityValueIndices_;
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
    std::vector< std::deque<int> > qvi_;
#endif

    // Pileup
    SAMPileupDeque samPileupDeque_;

    // Genotyper
    Genotyper genotyper_;

    // Quantizers
    std::map<int, Quantizer> quantizers_;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<SAMRecord> samRecordDeque_;
};

} // namespace calq

#endif // CALQ_QUALCODEC_QUALENCODER_H_

