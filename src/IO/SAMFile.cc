/** @file SAMFile.cc
 *  @brief This file contains the implementation of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "SAMFile.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include <stdio.h>

#include <string.h>

void SAMRecord::print(void) const
{
    LOG("qname: %s", qname.c_str());
    LOG("flag: %d", flag);
    LOG("rname: %s", rname.c_str());
    LOG("pos: %d", pos);
    LOG("mapq: %d", mapq);
    LOG("cigar: %s", cigar.c_str());
    LOG("rnext: %s", rnext.c_str());
    LOG("pnext: %d", pnext);
    LOG("tlen: %d", tlen);
    LOG("seq: %s", seq.c_str());
    LOG("qual: %s", qual.c_str());
    LOG("opt: %s", opt.c_str());

}


SAMFile::SAMFile(const std::string &path, const char *mode, const size_t &blockSize)
    : File(path, mode)
    , blockSize(blockSize)
    , line(NULL)
    , lineSize(sizeof(char) * (1 * MB))
{
    // For fgets and (f)seek to work together properly, a SAMFile instance
    // has to be opened in binary mode
    // TODO: fgetpos/fsetpos
//     if (strcmp(mode, "rb") != 0) {
//         throwErrorException("rb mode required");
//     }

    // 1 million chars should be enough
    line = (char *)malloc(lineSize);

    bool foundHeader = false;

    // Read SAM header
    size_t alignmentSectionBegin = tell();
    for (;;) {
        alignmentSectionBegin = tell();
        if (fgets(line, lineSize, fp) != NULL) {
            if (line[0] == '@') {
                header += line;
            } else {
                
                break;
            }
        } else {
            throwErrorException("Could not read SAM header");
        }
    }
    seek(alignmentSectionBegin);
    if (header.empty() == true) {
        throwErrorException("SAM header is missing");
    }
    LOG("SAM header: %s", header.c_str());

    readBlock();
}

SAMFile::~SAMFile(void)
{
    free(line);
}

void SAMFile::readBlock(void) {
    block.clear();
    
    for (size_t i = 0; i < blockSize; i++) {
        if (fgets(line, lineSize, fp) != NULL) {
            SAMRecord samRecord;
            
            size_t l = strlen(line) - 1;

            while (l && (line[l] == '\r' || line[l] == '\n')) {
                line[l--] = '\0';
            }
            
            
            char *c =  line;
            char *pc = c;
            int f = 0;

            while (*c) {
                if (*c == '\t') {
                    *c = '\0';
                    
                    if (f == 0) samRecord.qname = pc;
                    if (f == 1) samRecord.flag = (uint16_t)atoi(pc);
                    if (f == 2) samRecord.rname = pc;
                    if (f == 3) samRecord.pos = (uint32_t)atoi(pc);
                    if (f == 4) samRecord.mapq = (uint8_t)atoi(pc);
                    if (f == 5) samRecord.cigar = pc;
                    if (f == 6) samRecord.rnext = pc;
                    if (f == 7) samRecord.pnext = (uint32_t)atoi(pc);
                    if (f == 8) samRecord.tlen = (int64_t)atoi(pc);
                    if (f == 9) samRecord.seq = pc;
                    if (f == 10) samRecord.qual = pc;
                    if (f == 11) samRecord.opt = pc;
                    f++;
                    //*c = '\0';
                    if (f == 12) break;
                    
                    pc = c+1;
                }
                
                c++;
            }

            if (f == 11) samRecord.opt = pc;
            
            
            
            
            block.push_back(samRecord);
        } else {
            LOG("Truncated block");
            break;
        }
    }
}

