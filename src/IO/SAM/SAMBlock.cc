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
    , m_numMappedRecords(0)
    , m_numUnmappedRecords(0)
{
    // empty
}

cq::SAMBlock::~SAMBlock(void)
{
    // empty
}

size_t cq::SAMBlock::numMappedRecords(void) const
{
    return m_numMappedRecords;
}

size_t cq::SAMBlock::numUnmappedRecords(void) const
{
    return m_numUnmappedRecords;
}

size_t cq::SAMBlock::numRecords(void) const
{
    return (m_numMappedRecords + m_numUnmappedRecords);
}

void cq::SAMBlock::reset(void)
{
    records.clear();
    m_numMappedRecords = 0;
    m_numUnmappedRecords = 0;
}

