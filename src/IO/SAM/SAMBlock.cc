/** @file SAMBlock.cc
 *  @brief This file contains the implementation of the SAMBlock class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAM/SAMBlock.h"

cq::SAMBlock::SAMBlock(void)
    : records()
    , nrMappedRecords_(0)
    , nrUnmappedRecords_(0)
{
    // empty
}

cq::SAMBlock::~SAMBlock(void)
{
    // empty
}

size_t cq::SAMBlock::nrMappedRecords(void) const
{
    return nrMappedRecords_;
}

size_t cq::SAMBlock::nrUnmappedRecords(void) const
{
    return nrUnmappedRecords_;
}

size_t cq::SAMBlock::nrRecords(void) const
{
    return (nrMappedRecords_ + nrUnmappedRecords_);
}

void cq::SAMBlock::reset(void)
{
    records.clear();
    nrMappedRecords_ = 0;
    nrUnmappedRecords_ = 0;
}

