#include "BaseSpreader.h"
#include <iostream>


double BaseSpreader::push(double score, size_t softclips) {
    if(softclips >= MIN_HQ_SOFTCLIPS){
        int clipped = std::min(softclips, MAX_PROPAGATION);
        forwardSpread.push_back(std::pair<size_t, double>(clipped+1,score));
        for(int i=buffer.size()-clipped;i<buffer.size();++i){
            buffer[i] += score;

        }
    }

    double ownscore = score;

    for(auto it = forwardSpread.begin();it != forwardSpread.end();){
        ownscore += (*it).second;
        --(*it).first;
        if((*it).first == 0)
            it = forwardSpread.erase(it);
        else
            ++it;
    }

    return buffer.push(ownscore);
}

size_t BaseSpreader::getOffset() const {
    return MAX_PROPAGATION;
}

BaseSpreader::BaseSpreader(size_t max_prop, size_t min_hq_clips) : buffer(max_prop,0.0), MAX_PROPAGATION(max_prop), MIN_HQ_SOFTCLIPS(min_hq_clips) {

}
