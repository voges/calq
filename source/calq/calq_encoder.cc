#include "calq/calq_encoder.h"

#include <chrono>

#include "calq/error_exception_reporter.h"
#include "calq/log.h"
#include "calq/helpers.h"
#include "calq/fasta_file.h"
#include "calq/qual_encoder.h"
#include "calq/probability_distribution.h"
#include "calq/uniform_min_max_quantizer.h"
#include "calq/lloyd_max_quantizer.h"
#include "calq/sam_file.h"


namespace calq {

uint32_t computeLength(const std::string& cigar){
    // Compute 0-based first position and 0-based last position this record
    // is mapped to on the reference used for alignment
    uint32_t posMax = 0;

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    uint32_t opLen = 0;  // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (uint32_t) cigar[cigarIdx] - (uint32_t) '0';
            continue;
        }
        switch (cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                posMax += opLen;
                break;
            case 'I':
            case 'S':
                break;
            case 'D':
            case 'N':
                posMax += opLen;
                break;
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }
    posMax -= 1;
    return posMax;
}

void encode(const EncodingOptions& opt,
            const EncodingSideInformation& sideInformation,
            const EncodingBlock& input,
            DecodingBlock *output
){
    ProbabilityDistribution pdf(opt.qualityValueMin, opt.qualityValueMax);

    // Check quality value range
    for (auto const& samRecord : input.qvalues)
    {
        for (auto const& q : samRecord)
        {
            if ((static_cast<int>(q) - opt.qualityValueOffset) < opt.qualityValueMin)
            {
                throwErrorException("Quality value too small");
            }
            if ((static_cast<int>(q) - opt.qualityValueOffset) > opt.qualityValueMax)
            {
                throwErrorException("Quality value too large");
            }
            pdf.addToPdf((static_cast<size_t>(q) - opt.qualityValueOffset));
        }
    }

    std::map<int, Quantizer> quantizers;

    for (auto i = static_cast<int>(opt.quantizationMin); i <= static_cast<int>(opt.quantizationMax); ++i)
    {
        if (opt.quantizerType == EncodingOptions::QuantizerType::UNIFORM)
        {
            UniformMinMaxQuantizer quantizer(
                    static_cast<const int&>(opt.qualityValueMin),
                    static_cast<const int&>(opt.qualityValueMax), i
            );
            quantizers.insert(
                    std::pair<int, Quantizer>(
                            static_cast<const int&>(i - opt.quantizationMin),
                            quantizer
                    ));
        }
        else if (opt.quantizerType == EncodingOptions::QuantizerType::LLOYD_MAX)
        {
            LloydMaxQuantizer quantizer(static_cast<size_t>(i));
            quantizer.build(pdf);
            quantizers.insert(
                    std::pair<int, Quantizer>(
                            static_cast<const int&>(i - opt.quantizationMin),
                            quantizer
                    ));
        }
        else
        {
            throwErrorException("Quantization Type not supported");
        }
    }

    // Encode the quality values
    QualEncoder qualEncoder(opt, quantizers);
    for (size_t i = 0; i < sideInformation.positions.size(); ++i)
    {
        EncodingRead r = {sideInformation.positions[i], sideInformation.positions[i] +  computeLength(sideInformation.cigars[i]), input.qvalues[i], sideInformation.cigars[i], sideInformation.sequences[i]};
        qualEncoder.addMappedRecordToBlock(r, sideInformation.reference.substr(r.posMin - sideInformation.positionStart, r.posMax - r.posMin));
    }
}

}  // namespace calq
