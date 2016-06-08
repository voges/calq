/** @file MappedRecord.h
 *  @brief This file contains the definition of the MappedRecord class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef MAPPEDRECORD_H
#define MAPPEDRECORD_H

#include "Parsers/SAMRecord.h"
#include <string>
#include <vector>

/** @brief Class: MappedRecord
 *
 *  This class holds mapped records. Mapped records are a shortened version
 *  am SAM records, only containing alignment/mapping information.
 */
class MappedRecord {
public:
    MappedRecord(const SAMRecord &samRecord, const uint32_t &positionOffset);
    ~MappedRecord(void);

    uint32_t firstPosition;
    uint32_t lastPosition;
    uint32_t positionOffset;
    std::string nucleotides;
    std::string qualityValues;
    std::string cigar;

    void extractObservations(std::vector<std::string> &observedNucleotides,
                             std::vector<std::string> &observedQualityValues);
};

#endif // MAPPEDRECORD_H

