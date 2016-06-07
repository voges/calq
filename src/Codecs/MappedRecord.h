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

#include <string>

/** @brief Class: MappedRecord
 *
 *  This class holds mapped records. Mapped records are a shortened version
 *  am SAM records, only containing alignment/mapping information.
 */
class MappedRecord {
public:
    MappedRecord(const uint32_t &pos,
                 const std::string &cigar,
                 const std::string &seq,
                 const std::string &qual);
    ~MappedRecord(void);

    uint32_t firstPos;
    uint32_t lastPos;
    std::string alignedNucleotides;
    std::string insertedNucleotides;
    std::string alignedQualityValues;
    std::string insertedQualityValues;
    std::string cigar;
};

#endif // MAPPEDRECORD_H

