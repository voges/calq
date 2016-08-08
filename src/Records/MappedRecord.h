/** @file MappedRecord.h
 *  @brief This file contains the definition of the MappedRecord class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef MAPPEDRECORD_H
#define MAPPEDRECORD_H

#include "Records/SAMRecord.h"
#include <string>
#include <vector>

/** @brief Class: MappedRecord
 *
 *  This class holds mapped records. Mapped records are a shortened version
 *  am SAM records, only containing alignment/mapping information.
 */
class MappedRecord {
public:
    MappedRecord(const SAMRecord &samRecord);
    ~MappedRecord(void);

    uint32_t posMin;
    uint32_t posMax;
    std::string nucleotides;
    std::string qualityValues;
    std::string cigar;

    void extractObservations(const uint32_t &observedPosMin,
                             std::vector<std::string> &observedNucleotides,
                             std::vector<std::string> &observedQualityValues);
};

#endif // MAPPEDRECORD_H

