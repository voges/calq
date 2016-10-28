/** @file QuantCodec.cc
 *  @brief This file contains the implementation of the QuantEncoder,
 *         QuantEncoderScore, QuantEncoderUniform, QuantDecoder classes.
 *  @author Franck Awounang (awounang)
 *  @bug No known bugs
 */


#include "QuantCodec.h"
#include <cmath>
#include <string>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <iostream>



QuantEncoder::QuantEncoder():intervals(), levels(){
    //Empty
};

int QuantEncoder::getIndex(int v) {
    for (int i=0; i<levels.size(); i++) {
        if (v==levels[i]) {
            return i;
        }
    }
    for (int i=intervals.size()-1; i>=0; i--) {
        if (v >= intervals[i]) {
            return i+1;
        }
    }
    return 0;
}

int QuantEncoder::getValue(int v) {
    return levels[getIndex(v)];
}

const std::vector<int>& QuantEncoder::getLevels() {
    return levels;
}


QuantEncoder::~QuantEncoder() {
    // empty
}

QuantEncoderUniform::QuantEncoderUniform(const std::string& qual, int nbLevels){
    qMin = qual[0];
    qMax = qual[0];
    for (size_t i=0; i<qual.size(); i++) {
        if (qMin>qual[i]) {
            qMin = qual[i];
        } else if (qMax<qual[i]) {
            qMax = qual[i];
        }
    }
    update(qMin, qMax, nbLevels);
}
    
QuantEncoderUniform::QuantEncoderUniform(const std::vector<int>& data, int nbLevels){
    qMin = data[0];
    qMax = data[0];
    for (size_t i=0; i<data.size(); i++) {
        if (qMin>data[i]) {
            qMin = data[i];
        } else if (qMax<data[i]) {
            qMax = data[i];
        }
    }
    update(qMin, qMax, nbLevels);
}

QuantEncoderUniform::QuantEncoderUniform(int qMin, int qMax, int nbLevels) {
    update(qMin, qMax, nbLevels);
}

void QuantEncoderUniform::update(int qMin, int qMax, int nbLevels) {
    this->qMin = qMin;
    this->qMax = qMax;
    for (int i=0; i<nbLevels; i++) {
        double value = qMin + (qMax - qMin) * (i+0.5)/nbLevels;
        levels.push_back(value);
    }
    for (size_t i=1; i<levels.size(); i++) {
        intervals.push_back((levels[i] + levels[i-1])/2.0);
    }
}


QuantEncoderSemiUniform::QuantEncoderSemiUniform(const std::string& qual, int nbLevels) {
    int dMin = qual[0];
    int dMax = qual[0];
    for (size_t i=1; i<qual.size(); i++) {
        if (dMin > qual[i]) {
            dMin = qual[i];
        } else if (dMax < qual[i]) {
            dMax = qual[i];
        }
    }
    qMin = dMin;
    qMax = dMax;
    std::vector<int> occ;
    std::vector<float> scores;
    for (size_t i=0; i<nbLevels; i++) {
        scores.push_back(0);
        occ.push_back(0);
    }
    float dist = dMax-dMin+1.0;
    for (int d: qual) {
        int idx = (int)((d-dMin)*nbLevels/dist);
        scores[idx] += d;
        occ[idx] += 1;
    }
    update(scores, occ, dMin, dMax, nbLevels);
}
QuantEncoderSemiUniform::QuantEncoderSemiUniform(const std::vector<int>& data, int nbLevels) {
    int dMin = data[0];
    int dMax = data[0];
    for (size_t i=1; i<data.size(); i++) {
        if (dMin > data[i]) {
            dMin = data[i];
        } else if (dMax < data[i]) {
            dMax = data[i];
        }
    }
    qMin = dMin;
    qMax = dMax;
    std::vector<int> occ;
    std::vector<float> scores;
    for (size_t i=0; i<nbLevels; i++) {
        scores.push_back(0);
        occ.push_back(0);
    }
    float dist = dMax-dMin+1.0;
    for (int d: data) {
        int idx = (int)((d-dMin)*nbLevels/dist);
        scores[idx] += d;
        occ[idx] += 1;
    }
    update(scores, occ, dMin, dMax, nbLevels);
}

void QuantEncoderSemiUniform::update(std::vector<float>& values, std::vector<int>& occurencies, int& dMin, int& dMax, int nbLevels) {
    float dist = dMax - dMin + 1.0;
    levels.clear();
    for (size_t i=0; i<nbLevels; i++) {
        if (occurencies[i] != 0) {
            levels.push_back((int) (values[i]/occurencies[i]+0.5));
        } else {
            levels.push_back((int)(dMin + (0.5+i)*(dMax-dMin)/nbLevels+0.5));
        }
    }
    for (size_t i=1; i<levels.size(); i++) {
        intervals.push_back((levels[i] + levels[i-1])/2.0);
    }
}

QuantEncoderScore::QuantEncoderScore(const std::vector<int>& data, int nbLevels) {
    int dMin = data[0];
    int dMax = data[0];
    for (size_t i=1; i<data.size(); i++) {
        if (dMin > data[i]) {
            dMin = data[i];
        } else if (dMax < data[i]) {
            dMax = data[i];
        }
    }
    std::vector<float> scores(dMax-dMin+1, 0.0);
    std::vector<float> pmf;
    const size_t qualSize(data.size());
    for (size_t i=0; i<data.size(); i++) {
        scores[(int)data[i]-dMin] += 1;
    }
    for (size_t i=0; i<scores.size(); i++) {
        if (scores[i]>0) {
            float p = scores[i]/qualSize;
            pmf.push_back(p);
            levels.push_back(i);           
        }
    }
    update(pmf, nbLevels, dMin);          
}

QuantEncoderScore::QuantEncoderScore(const std::string& qual, int nbLevels) {
    int dMin = qual[0];
    int dMax = qual[0];
    for (size_t i=1; i<qual.size(); i++) {
        if (dMin > qual[i]) {
            dMin = qual[i];
        } else if (dMax < qual[i]) {
            dMax = qual[i];
        }
    }
    std::vector<float> scores(dMax-dMin+1, 0.0);
    std::vector<float> pmf;
    const size_t qualSize(qual.size());
    for (size_t i=0; i<qual.size(); i++) {
        scores[(int)qual[i]-dMin] += 1;
    }
    for (size_t i=0; i<scores.size(); i++) {
        if (scores[i]>0) {
            float p = scores[i]/qualSize;
            pmf.push_back(p);
            levels.push_back(i);           
        }
    }
    update(pmf, nbLevels, dMin);
}
    
void QuantEncoderScore::update(std::vector<float> data, int nbLevels, int offset) {    
    std::vector<float> pmf(data);
    std::vector<float> scores(pmf.size(), 0.0);
    //float pBonus=0;
    std::vector<int> topLevels;
    //TODO find a way of selecting the best N instead of eleminating the worse M
    if ((size_t)nbLevels > pmf.size()/2 || true) {
        int nbRounds = pmf.size() - nbLevels;
        for (int round=0; round<nbRounds; round++) {
            scores.clear();
            for (size_t idx=0; idx<pmf.size(); idx++) {
                float score = 0;
                float dist = 0;
                for (size_t idx2=0; idx2<pmf.size(); idx2++) {
                    if (idx!=idx2) {
                        dist += std::pow(std::abs(levels[idx2]-levels[idx]), -3);
                    }
                }
                //TODO use neighbour's probabilies in score
                /*
                if (idx==0) {
                    pBonus = 0.1 * pmf[1];
                } else if (idx==pmf.size()-1) {
                    pBonus = 0.1 * pmf[idx];
                } else {
                    pBonus = 0.05 * (pmf[idx-1] + pmf[idx+1]);
                }
                */
                dist = 1/dist;
                score = (pmf[idx]) * dist;
                scores.push_back(score);
            }
            unsigned int idMin = 0;
            for (size_t idx2=0; idx2<scores.size(); idx2++) {
                if (scores[idMin] > scores[idx2]) {
                    idMin = idx2;
                }
            }
            //TODO redistribute pmf[idMin] to neighbours, according to neighbour's probabilies.
            /*
            if (idMin==0) {
                pmf[1] += pmf[0];
            } else if (idMin==pmf.size()-1) {
                pmf[idMin-1] = pmf[idMin];
            } else {
                if (pmf[idMin-1] == pmf[idMin+1]) {
                    pmf[idMin-1] += 0.5 * pmf[idMin];
                    pmf[idMin+1] += 0.5 * pmf[idMin];
                } else {
                    float coef = pmf[idMin-1]/(pmf[idMin-1] + pmf[idMin+1]);
                    pmf[idMin-1] += coef * pmf[idMin];
                    pmf[idMin+1] += (1-coef) * pmf[idMin];
                }
            }
            */
            pmf.erase(pmf.begin()+idMin);
            levels.erase(levels.begin()+idMin);
            scores.erase(scores.begin()+idMin);      
        }
    } else {
        for (int round=0; round<nbLevels; round++) {
            scores.clear();
            for (size_t idx=0; idx<pmf.size(); idx++) {
                float score = 0;
                float dist = 0;
                for (size_t idx2=0; idx2<pmf.size(); idx2++) {
                    if (idx!=idx2) {
                        dist += std::pow(std::abs(levels[idx2]-levels[idx]), -3);
                    }
                }
                dist = 1/dist;
                score = (pmf[idx]) * dist;
                scores.push_back(score);
            }
            unsigned int idMax = 0;
            for (size_t idx2=1; idx2<scores.size(); idx2++) {
                if (scores[idMax] < scores[idx2]) {
                    idMax = idx2;
                }
            }
            int val = levels[idMax];
            topLevels.push_back(val);
            pmf.erase(pmf.begin()+idMax);
            levels.erase(levels.begin()+idMax);
            scores.erase(scores.begin()+idMax);
        }
        std::sort(topLevels.begin(), topLevels.end());
        levels = topLevels;
    }
    for (size_t i=0; i<levels.size(); i++) {
        levels[i] += offset;
    }
    for (size_t i=1; i<levels.size(); i++) {
        intervals.push_back((levels[i] + levels[i-1])/2.0);
    }
}

QuantEncoderLloydMax::QuantEncoderLloydMax(const std::vector<int>& data, int nbLevels) {
    std::vector<Item> items;
    int dMin = data[0];
    int dMax = data[0];
    for (size_t i=0; i<data.size(); i++) {
        Item item(data[i], 0);
        items.push_back(item);
        if (data[i] < dMin) {
            dMin = data[i];
        } else if (data[i] > dMax) {
            dMax = data[i];
        }
    }
    update(items, nbLevels, dMin, dMax-dMin+1);
}


QuantEncoderLloydMax::QuantEncoderLloydMax(const std::string& qual, int nbLevels) {
    std::vector<Item> items;
    int dMin = qual[0];
    int dMax = qual[0];
    for (size_t i=0; i<qual.size(); i++) {
        Item item(qual[i], 0);
        items.push_back(item);
        if (qual[i] < dMin) {
            dMin = qual[i];
        } else if (qual[i] > dMax) {
            dMax = qual[i];
        }
    }
    update(items, nbLevels, dMin, dMax-dMin+1);
}


void QuantEncoderLloydMax::update(std::vector<Item>& items, int nbLevels, int offset,
                             int alphabetSize) {
    const int NB_ROUNDS_MAX = 10;
    const int NB_ROUNDS_BASE = 6;
    int nbRounds = NB_ROUNDS_BASE;
    std::vector<int> toDel;
    std::vector<double> clusters;
    int nbBins = (nbLevels * 1);    // replacing 2*nbLevels instead of 1 improves performance for poorer MSE
    for (int i=0; i<nbBins; i++) {
        int lvl = offset + 0.5 + ((i+0.5)*(alphabetSize))/nbBins;
        clusters.push_back(lvl);
    }
    for (size_t round=0; round<nbRounds; round++) {
        // update items;
        for (size_t i=0; i<items.size(); i++) {
            int cluster = 0;
            for(size_t j=1; j<clusters.size(); j++) {
                if (std::abs(items[i].value - clusters[cluster]) > std::abs(items[i].value - clusters[j])) {
                    cluster = j;
                }
            }
            items[i].cluster = cluster;
        }
        // update clusters;
        for (size_t k=0; k<clusters.size(); k++) {
            double score = 0;
            int cpt = 0;
            for (size_t i=0; i<items.size(); i++) {
                if (items[i].cluster == k) {
                    score += items[i].value;
                    cpt++;
                }
            }
            if (cpt==0) {
                toDel.push_back(k);
            } else {
                double lvlF = score/cpt;
                clusters[k] = lvlF;
            }
        }
        for (int i=toDel.size()-1; i>=0; i--) {
            clusters.erase(clusters.begin()+toDel[i]);
        }
        
        if (clusters.size() < nbLevels) {
            int steps = nbLevels-clusters.size();
            for (int i=0; i<steps; i++) {
                int lvl = offset + rand() % alphabetSize;
                clusters.push_back(lvl);
            }
            if (round == NB_ROUNDS_BASE-1) {
                nbRounds = NB_ROUNDS_MAX;
            }
        }
        toDel.clear();
    }
    
    std::sort(clusters.begin(), clusters.end()); 
    for (size_t i=0; i<clusters.size(); i++) {
        levels.push_back((int)(0.5 + clusters[i]));
    }
    for (size_t i=1; i<levels.size(); i++) {
        intervals.push_back((levels[i] + levels[i-1] )/2.0);
    }
}




QuantDecoder::QuantDecoder(const std::vector<int>& pLevels) {
    for (int value: pLevels) {
        levels.push_back(value);
    }
}

int QuantDecoder::decode(int value) {
    return levels[value];
}

QuantDecoder::~QuantDecoder() {
}
