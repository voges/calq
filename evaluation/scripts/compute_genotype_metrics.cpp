/*
 * Copyright (c) 2015, Fonleap Ltd
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this
 *    list of conditions and the following disclaimer in the documentation and/or
 *    other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <errno.h>

#define COMPARE_SNP   0
#define COMPARE_INDEL 1
#define COMPARE_ALL   2

// Parse reference FA file
class ParseRefFA {
protected:
    uint8_t *refSeq;
    unsigned int refSeqSize;
    std::map<std::string, uint64_t> chromStartPos;

public:
    ParseRefFA() {
        refSeq = NULL;
        refSeqSize = 0;
    }

    ~ParseRefFA() {
        if (refSeq)
            free(refSeq);
    }

    bool readFile(const char *filename) {
        char buf[65536];

        FILE *file = fopen(filename, "rb");
        if (file == NULL)
            return false;

        refSeqSize = 32*1024*1024;
        refSeq = (uint8_t *)malloc(32*1024*1024);
        uint64_t seqPos = 0;

        while (fgets(buf, 65536, file) != NULL) {
            // read the line
            if (buf[0] == '>') {
                // marks the chromosome
                for (int i=1; i<65536; ++i) {
                    if (buf[i]==' ' || buf[i]=='\t') {
                        buf[i]=0;
                        break;
                    }
                }
                std::string chrName;
                if (buf[1]=='c' && buf[2]=='h' && buf[3]=='r')
                    chrName = buf+4;
                else
                    chrName = buf+1;
                chromStartPos[chrName] = seqPos;
                chromStartPos["chr" + chrName] = seqPos;
            }
            else {
                if (seqPos+65536 > refSeqSize) {
                    refSeqSize += 32*1024*1024;
                    uint8_t *newBuf = (uint8_t *)realloc(refSeq, refSeqSize);
                    if (newBuf == NULL)
                        return false;
                    refSeq = newBuf;
                }
                for (int i=0; i<65536; ++i) {
                    if (buf[i] == 0 || buf[i] == '\n' || buf[i] == '\r')
                        break;
                    if (buf[i] > 'Z') // convert lower case to upper case
                        buf[i] = buf[i] - 'a' + 'A';
                    refSeq[seqPos++] = buf[i];
                }
            }
        }
        return true;
    }

    // get sequence string immediately preceding chromosome:pos of length len
    std::string getPrefix(std::string &chromosome, uint64_t pos, uint64_t len) {
        std::map<std::string, uint64_t>::iterator it = chromStartPos.find(chromosome);
        if (it == chromStartPos.end())
            return "";
        uint64_t absPos = it->second + pos - 1;
        if (len > absPos)
            len = absPos;
        std::string prefix((const char *)(&refSeq[absPos-len]), len);
        return prefix;
    }
};

/* For example:
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3
*/

struct Variant {
    int chrom;
    uint64_t pos;
    unsigned int srcLen;
    std::string dst;
    double qual;
    bool isTruePositive;

    Variant() {
        chrom = 0;
        pos = 0;
        srcLen = 0;
        dst.clear();
        qual = 0.0;
        isTruePositive = false;
    }

    unsigned int commonHead(const std::string &ref) {
        unsigned int minLen = dst.size();
        if (minLen > ref.size())
            minLen = ref.size();
        for (unsigned int i=0; i<minLen; ++i) {
            if (ref[i] != dst[i])
                return i;
        }
        return minLen;
    }

    unsigned int commonTail(const std::string &ref) {
        unsigned int minLen = dst.size();
        if (minLen > ref.size())
            minLen = ref.size();
        for (unsigned int i=0; i<minLen; ++i) {
            if (ref[ref.size()-i-1] != dst[dst.size()-i-1])
                return i;
        }
        return minLen;
    }

    std::string canonicalizeR(const std::string &ref) {
        unsigned int commonTailLen = commonTail(ref);
        if (commonTailLen==0)
            return ref;
        dst.resize(dst.size()-commonTailLen);
        srcLen -= commonTailLen;
        return ref.substr(0, ref.size()-commonTailLen);
    }

    void canonicalizeL(const std::string &ref) {
        // remove any commonality at head and tail
        unsigned int commonHeadLen = commonHead(ref);
        if (commonHeadLen==0)
            return;
        dst = dst.substr(commonHeadLen);
        pos += commonHeadLen;
        srcLen -= commonHeadLen;
    }

    bool hasCommonSuffix(const std::string &ref) {
        return (ref[ref.size()-1] == dst[dst.size()-1]);
    }

    /* Converts variant into a canonical form. Here are some examples:
     1 100727133 . CCTCT CCT,C 152 => (+3,CT,{}),(+1,CTCT,{})
     1 10087983 . CCTGA TCTGA,C 43 => (0,C,T),(+1,CTGA,{})
     1 100926991 . CT TT,C 1066 => (0,C,T),(+1,T,{})
     1 100976646 . T TGTTCA,TGTTCG 124 => (+1,{},GTTCA),(+1,{},GTTCG}
    */
    void canonicalize(const std::string &ref) {
        canonicalizeR(ref);
        canonicalizeL(ref);
    }

    void prepend(const std::string &prefix) {
        unsigned int len = prefix.size();
        dst = prefix + dst;
        pos -= len;
        srcLen += len;
    }

    inline bool operator < (const Variant &X) const {
        if (pos != X.pos)
            return pos < X.pos;
        if (chrom != X.chrom)
            return chrom < X.chrom;
        if (srcLen != X.srcLen)
            return srcLen < X.srcLen;
        return dst < X.dst;
    }

    inline bool operator == (const Variant &X) const {
        if (pos != X.pos)
            return false;
        if (chrom != X.chrom)
            return false;
        if (srcLen != X.srcLen)
            return false;
        return dst == X.dst;
    }
};

struct ROCItem {
    unsigned int truePositives;
    unsigned int falsePositives;

    ROCItem() {
        truePositives = 0;
        falsePositives = 0;
    }

    void addTruePositive() {
        truePositives++;
    }

    void addFalsePositive() {
        falsePositives++;
    }

    void add(const ROCItem &X) {
        truePositives += X.truePositives;
        falsePositives += X.falsePositives;
    }
};

// Chromosome, position, ID, Ref, [Alt,...], Quality, Filter
class VCFParser {
protected:
    std::map<std::string, int> chromMap;
    std::set<std::string> ignoredChrom;
    uint64_t countTruePositives;

    // true if no error, false if error
    bool parseLine(const int compareType, char *cline, ParseRefFA *parseRefFA = NULL, VCFParser *positiveSet = NULL, std::set<Variant> *positiveVariants = NULL) {
        if (cline[0] == 0)
            return true;
        if (cline[0] == '#')
            return true;

        // tokenize
        char *saveptr, *token1, *token2, /* *token3,*/ *token4, *token5, *token6, *altToken;
        token1 = strtok_r(cline, " \t", &saveptr);
        token2 = strtok_r(NULL, " \t", &saveptr);
        strtok_r(NULL, " \t", &saveptr);
        token4 = strtok_r(NULL, " \t", &saveptr);
        token5 = strtok_r(NULL, " \t", &saveptr);
        token6 = strtok_r(NULL, " \t", &saveptr);
        if (token6 != NULL) {
            std::string chrom = token1;
            std::map<std::string, int>::iterator iter = chromMap.find(chrom);
            std::string refStr = token4;
            if (iter == chromMap.end()) {
                if (ignoredChrom.insert(chrom).second) {
                    //printf("Warning - unable to find chromosome [%s]\n", chrom.c_str());
                }
                return true;
            }
            else {
                Variant v, v2;
                v.chrom = iter->second;
                v.pos = atol(token2);
                v.srcLen = refStr.size();
                v.qual = atof(token6);
                altToken = strtok_r(token5, ",", &saveptr);
                while (altToken != NULL) {
                    v2 = v;
                    v2.dst = altToken;
                    std::string refStr2 = refStr;
                    if (parseRefFA != NULL) {
                        // use reference sequence to normalise to a canonical form with minimum possible 'pos'
                        while (v2.hasCommonSuffix(refStr2)) {
                            std::string prefix = parseRefFA->getPrefix(chrom, v2.pos, 100);
                            // prepend both the variant and reference with the real prefix sequence data
                            v2.prepend(prefix);
                            refStr2 = v2.canonicalizeR(prefix + refStr2); // discard RHS common portion
                        }
                    }
                    // final RHS and LHS canonicalization
                    v2.canonicalize(refStr2);
                    if (positiveSet != NULL) {
                        v2.isTruePositive = (positiveSet->variantDB.find(v2) != positiveSet->variantDB.end());
                    }
                    altToken = strtok_r(NULL, ",", &saveptr);
                    bool isSNP = (v2.srcLen == 1 && v2.dst.size() == 1);
                    if (compareType == COMPARE_SNP && !isSNP)
                        continue;
                    if (compareType == COMPARE_INDEL && isSNP)
                        continue;
                    bool isNewVariant = variantDB.insert(v2).second;
                    if (v2.isTruePositive && isNewVariant) {
                        countTruePositives++;
                        if (positiveVariants != NULL)
                            positiveVariants->insert(v2);
                    }
                }
            }
        }
        return true;
    }

    static bool _sortFn(const Variant &X, const Variant &Y) {
        if (X.qual != Y.qual)
            return X.qual > Y.qual;
        return X.pos < Y.pos;
    }

public:
    // need to convert chromosome and position to an actual full position
    std::set<Variant> variantDB;

    VCFParser() {
        // ensure consistent mapping of chromosome names
        chromMap["X"] = -1;
        chromMap["chrX"] = -1;
        chromMap["Y"] = -2;
        chromMap["chrY"] = -2;
        chromMap["M"] = -3;
        chromMap["chrM"] = -3;
        chromMap["MT"] = -3;
        chromMap["chrMT"] = -3;
        for (unsigned int i=1; i<1000; ++i) {
            char namebuf[256];
            snprintf(namebuf, 256, "%u", i);
            std::string name = namebuf;
            chromMap[name] = i;
            chromMap["chr"+name] = i;
        }
        countTruePositives = 0;
    }

    size_t size() const {
        return variantDB.size();
    }

    uint64_t numTruePositives() const {
        return countTruePositives;
    }

    bool parseFile(const int compareType, const std::string &filename, ParseRefFA *parseRefFA = NULL, VCFParser *positiveSet = NULL, std::set<Variant> *positiveVariants = NULL) {
        FILE *in = fopen(filename.c_str(), "r");
        if (in == NULL) {
            printf("Unable to open input file: %s [%s]\n", filename.c_str(), strerror(errno));
            return false;
        }
        char line[4*1024*1024];
        while (fgets(line, sizeof(line), in)) {
            if (!parseLine(compareType, line, parseRefFA, positiveSet, positiveVariants)) {
                fclose(in);
                return false;
            }
        }
        fclose(in);
        return true;
    }

    // returns sorted variants from largest quality to smallest quality
    void getQualSortedVariants(std::map<double, ROCItem> *sorted) const {
        if (sorted == NULL)
            return;
        sorted->clear();
        for (std::set<Variant>::iterator it = variantDB.begin(); it != variantDB.end(); ++it) {
            if (it->isTruePositive)
                (*sorted)[it->qual].addTruePositive();
            else
                (*sorted)[it->qual].addFalsePositive();
        }

        std::map<double, std::set<Variant> > holepos, holeneg;
        for (std::set<Variant>::iterator it = variantDB.begin(); it != variantDB.end(); ++it) {
            if (it->isTruePositive)
                holepos[it->qual].insert(*it);
            else
                holeneg[it->qual].insert(*it);
        }

/*
        for (std::set<Variant>::iterator iter = holeneg.rbegin()->second.begin(); iter != holeneg.rbegin()->second.end(); ++iter) {
            printf("NEG: qual %lf chrom %d pos %"PRId64" srcLen %d dst %s\n", iter->qual, iter->chrom, iter->pos, iter->srcLen, iter->dst.c_str());
        }
        for (std::set<Variant>::iterator iter = holepos.rbegin()->second.begin(); iter != holepos.rbegin()->second.end(); ++iter) {
            printf("POS: qual %lf chrom %d pos %"PRId64" srcLen %d dst %s\n", iter->qual, iter->chrom, iter->pos, iter->srcLen, iter->dst.c_str());
        }
*/
    }
};

int main(int argc, char *argv[]) {
    ParseRefFA *refFile = NULL;
    if (argc < 5) {
        goto failed;
    }
    int compareType;
    if (strcmp(argv[1], "snp") == 0) {
        compareType = COMPARE_SNP;
    }
    else if (strcmp(argv[1], "indel") == 0) {
        compareType = COMPARE_INDEL;
    }
    else if (strcmp(argv[1], "all") == 0) {
        compareType = COMPARE_ALL;
    }
    else {
        goto failed;
    }

    // read fa file
    if (strcmp(argv[2], "none") != 0) {
        refFile = new ParseRefFA();
        if (!refFile->readFile(argv[2])) {
            printf("Error reading file %s: %s\n", argv[2], strerror(errno));
            goto errfailed;
        }
    }

    if (argc == 5) {
        VCFParser goldSet, sampleSet;
        if (!goldSet.parseFile(compareType, argv[3], refFile)) {
            printf("Error - unable to parse file %s\n", argv[3]);
            goto errfailed;
        }
        if (!sampleSet.parseFile(compareType, argv[4], refFile, &goldSet)) {
            printf("Error - unable to parse file %s\n", argv[4]);
            goto errfailed;
        }

        uint64_t trueSetSize = goldSet.size();
        uint64_t sampleSetSize = sampleSet.size();

        // items in sample that are in the gold set
        uint64_t numTruePositives = sampleSet.numTruePositives();

        // items in sample that are not in the gold set
        uint64_t numFalsePositives = sampleSetSize - numTruePositives;

        // items in gold set that are not in the sample
        uint64_t numFalseNegatives = trueSetSize - numTruePositives;

        printf("Precision = %lf, Recall = %lf, F-score = %lf\n", (double)numTruePositives / (double)sampleSetSize, (double)numTruePositives / (double)trueSetSize, (double)(2*numTruePositives)/(double)(2*numTruePositives+numFalsePositives+numFalseNegatives));
        printf("numTruePositives = %lf, sampleSetSize = %lf, trueSetSize = %lf\n", (double)numTruePositives, (double)sampleSetSize, (double)trueSetSize);
        return 0;
    }
    if (argc > 5) {
        VCFParser goldSet;
        std::vector<VCFParser *> sampleSets;

        if (!goldSet.parseFile(compareType, argv[3], refFile)) {
            printf("Error - unable to parse file %s\n", argv[3]);
            goto errfailed;
        }
        uint64_t trueSetSize = goldSet.size();

        unsigned int numSamples = argc-4;
        std::set<Variant> positiveVariants; // union of all true positives found across samples
        uint64_t sampleSetSize[numSamples];
        uint64_t numTruePositives[numSamples];
        uint64_t numFalsePositives[numSamples];
        uint64_t numFalseNegatives[numSamples];
        for (int i=0; i<numSamples; ++i) {
            sampleSets.push_back(new VCFParser());
            if (!sampleSets[i]->parseFile(compareType, argv[i+4], refFile, &goldSet, &positiveVariants)) {
                printf("Error - unable to parse file %s\n", argv[i+4]);
                goto errfailed;
            }

            sampleSetSize[i] = sampleSets[i]->size();
            // items in sample that are in the gold set
            numTruePositives[i] = sampleSets[i]->numTruePositives();
            // items in sample that are not in the gold set
            numFalsePositives[i] = sampleSetSize[i] - numTruePositives[i];
            // items in gold set that are not in the sample
            numFalseNegatives[i] = trueSetSize - numTruePositives[i];

        }
/*
        for (int i=0; i<numSamples; ++i)
            printf("%s ", argv[i+4]);
        printf("\n");
*/
        unsigned int numUnionPositiveVariants = positiveVariants.size();
        std::set<unsigned int> allSortedPos;
        std::map<unsigned int, unsigned int> sortedSamples[numSamples];
        double AUC[numSamples];

        allSortedPos.insert(0); // initial datapoint has no false positives
        for (int i=0; i<numSamples; ++i) {
            AUC[i] = 0.0;
            std::map<double, ROCItem> sorted;
            sampleSets[i]->getQualSortedVariants(&sorted);
            ROCItem cumm;
            for (std::map<double, ROCItem>::reverse_iterator iter = sorted.rbegin(); iter != sorted.rend(); ++iter) {
                cumm.add(iter->second);
                allSortedPos.insert(cumm.falsePositives);
                sortedSamples[i][cumm.falsePositives] = cumm.truePositives;
            }
            unsigned int prevTruePos = 0;
            unsigned int prevFalsePos = 0;
            for (std::map<unsigned int, unsigned int>::iterator it=sortedSamples[i].begin(); it != sortedSamples[i].end(); ++it) {
                // trapezoidal AUC calculation
                AUC[i] += ((double)(it->second + prevTruePos) * (double)(it->first - prevFalsePos)) / 2.0;
                prevTruePos = it->second;
                prevFalsePos = it->first;
            }
            sortedSamples[i][0] = 0;
        }
        unsigned int maxFalsePositives = *allSortedPos.rbegin();
        // interpolate from last given data point to (1.0,1.0) i.e. 100% true positive found at 100% false positive rate
        for (int i=0; i<numSamples; ++i) {
            // trapezoidal AUC calculation
            AUC[i] += ((double)(numUnionPositiveVariants + sortedSamples[i].rbegin()->second) * (double)(maxFalsePositives - sortedSamples[i].rbegin()->first)) / 2.0;
            sortedSamples[i][maxFalsePositives] = numUnionPositiveVariants;
        }
        // normalise x and y scale of AUC
        for (int i=0; i<numSamples; ++i) {
            AUC[i] /= numUnionPositiveVariants;
            AUC[i] /= maxFalsePositives;
        }
        // print AUC along with filename
        printf("NAMED_AUC:");
        for (int i=0; i<numSamples; ++i) {
            printf(" %s=%lf", argv[i+4],AUC[i]);
        }
        printf("\n");
        printf("AUC:");
        for (int i=0; i<numSamples; ++i) {
            printf(" %lf", AUC[i]);
        }
        printf("\n");
        printf("PREC:");
        for (int i=0; i<numSamples; ++i) {
            printf(" %lf", (double)numTruePositives[i] / (double)sampleSetSize[i]);
        }
        printf("\n");
        printf("RECALL:");
        for (int i=0; i<numSamples; ++i) {
            printf(" %lf", (double)numTruePositives[i] / (double)trueSetSize);
        }
        printf("\n");
        printf("FSCORE:");
        for (int i=0; i<numSamples; ++i) {
            printf(" %lf", (double)(2*numTruePositives[i])/(double)(2*numTruePositives[i]+numFalsePositives[i]+numFalseNegatives[i]));
        }
        printf("\n");
        printf("CSV[AUC/PREC/RECALL/FSCORE]:");
        for (int i=0; i<numSamples; ++i) {
            printf(" %s=%lf,%lf,%lf,%lf", argv[i+4] ,AUC[i],(double)numTruePositives[i] / (double)sampleSetSize[i], (double)numTruePositives[i] / (double)trueSetSize, (double)(2*numTruePositives[i])/(double)(2*numTruePositives[i]+numFalsePositives[i]+numFalseNegatives[i]));
        }
        printf("\n");
        for (std::set<unsigned int>::iterator iter = allSortedPos.begin(); iter != allSortedPos.end(); ++iter) {
            // for each sample find the linearly interpolated position
            unsigned int falsePositivePos = *iter;
            printf("%lf", (double)falsePositivePos / (double)maxFalsePositives);
            for (int i=0; i<numSamples; ++i) {
                std::map<unsigned int, unsigned int>::iterator it = sortedSamples[i].lower_bound(falsePositivePos);
                unsigned int foundFalsePos = it->first;
                unsigned int foundTruePos = it->second;
                if (foundFalsePos == falsePositivePos) {
                    printf(" %lf", (double)foundTruePos/(double)numUnionPositiveVariants);
                }
                else { // interpolate with previous data entry
                    --it;
                    unsigned int prevFoundFalsePos = it->first;
                    unsigned int prevFoundTruePos = it->second;
                    double interp = ((double)(falsePositivePos - prevFoundFalsePos)) / ((double)(foundFalsePos - prevFoundFalsePos));
                    double interpTruePos = interp * (double)foundTruePos + (1.0-interp) * (double)prevFoundTruePos;
                    printf(" %lf", interpTruePos/(double)numUnionPositiveVariants);
                }
            }
            printf("\n");
        }
        for (int i=0; i<numSamples; ++i) {
            delete sampleSets[i];
        }
    }
    return 0;
failed:;
    printf("Calculate Precision, Recall, F-score:\n\%s [all | snp | indel] [reference.fa | none] [goldset.vcf] [sample.vcf]\n\n", argv[0]);
    printf("Calculate ROC AUC plot:\n\%s [all | snp | indel] [reference.fa | none] [goldset.vcf] [[sample1.vcf] [sample2.vcf] ...]\n\n", argv[0]);
    printf("Note: To correctly match indels, the reference fasta that the variant files are based on should be specified.\n");
    printf("      First column of ROC AUC plot is false-positive rate, other columns are true-positive rate in order of sample\n");
    printf("      First lines of plot have text describing Precision, Recall, F-score and ROC AUC numbers for each sample\n");
    return 0;
errfailed:;
   if (refFile)
       delete refFile;
   return -1;
}
