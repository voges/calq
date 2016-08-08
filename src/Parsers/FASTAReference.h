/** @file FASTAReference.h
 *  @brief This file contains the definition of the FASTAReference class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef FASTAREFERENCE_H
#define FASTAREFERENCE_H

#include <string>

/** @brief Class: FASTAReference
 *
 *  This class holds a FASTA reference sequence consisting of a header and the
 *  actual sequence.
 */
class FASTAReference {
public:
    FASTAReference(void) 
        : header("")
        , sequence("")
    {}
    ~FASTAReference(void) {}

    std::string header;
    std::string sequence;
};

#endif // FASTAREFERENCE_H

