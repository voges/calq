/** @file SAMFile.h
 *  @brief This file contains the definition of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef SAMFILE_H
#define SAMFILE_H

#include "IO/File.h"
#include <deque>
#include <map>

class SAMRecord {
public:
    std::string qname; // Query template NAME
    uint16_t flag;   // bitwise FLAG (uint16_t)
    std::string rname; // Reference sequence NAME
    uint32_t pos;    // 1-based leftmost mapping POSition (uint32_t)
    uint8_t  mapq;   // MAPping Quality (uint8_t)
    std::string cigar; // CIGAR string
    std::string rnext; // Ref. name of the mate/NEXT read
    uint32_t pnext;  // Position of the mate/NEXT read (uint32_t)
    int64_t  tlen;   // observed Template LENgth (int64_t)
    std::string seq;   // segment SEQuence
    std::string qual;  // QUALity scores
    std::string opt;   // OPTional information

    void print(void) const;
};

class SAMFile : public File {
public:
    SAMFile(const std::string &path, const char *mode, const size_t &blockSize);
    ~SAMFile(void);

    void readBlock(void);

    std::string header;
    std::deque<SAMRecord> block;

private:
    size_t blockSize;
    char *line;
    size_t lineSize;
};

#endif // SAMFILE_H

