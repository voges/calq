Parameters Set
==============

The files *.mpeg_qv_parameters_set contain the parameters as defined in
Section 5 "Parameters Set" of document N16758.

Syntax:
-------
// The test data contains 1 QV per base.
uint32 QVIndexDimension = 0x00000001;

// QVIndex is in the set {0,1,2,3,4,5,6,7}.
uint32 QVIndexAlphabetSize = 0x00000008;

// The data was coded using 7 codebook identifiers in the set {0,1,2,3,4,5,6}.
// Note that the codebook identifier '-1' may also appear (this codebook
// identifier is sent at loci with zero sequencing depth). Thus, we incremented
// the descriptor QVCodebookIdentifier by 1. The descriptor stream therefore
// contains data in the set {0,1,2,3,4,5,6,7}.
uint32 QVCodebookIdentifierAlphabetSize = 0x00000007;

for (QVCodebookIdentifier = 0; QVCodebookIdentifier < QVCodebookIdentifierAlphabetSize; QVCodebookIdentifier++) {
    uint64 QVCodebookIdentifier;
    uint64 QVCodebookSize;
    for (QVIndex = 0; QVIndex < QVCodebookSize; QVIndex++) {
        uint8 QVIndex;
        uint8 QVReconstructed;
    }
}


Descriptor streams
==================

All (aligned) data was coded in the matrix mode defined in Section 9.1.2 of
document N16758.

The files *.mpeg_qvci_0 contain the descriptor QVCodebookIdentifier stored in
multiple blocks.

Syntax:
-------
while (EndOfFile == false) {
    // 0-based position offset of the current block.
    uint32 PositionOffset;

    // Size of the QVCodebookIdentifier descriptor stream in bytes.
    uint64 QVCodebookIdentifierBufferSize;

    // Buffer containing the QVCodebookIdentifier descriptor stream.
    // The data is stored in the ASCII format; thus the QVCodebookIdentifier
    // '0' is actually stored as 0x30 and similarly for the other values.
    uint8[QVCodebookIdentifierBufferSize] QVCodebookIdentifier;
}

The files *.mpeg_qvi_0 contain the descriptor QVIndex stored in multiple blocks.

Syntax:
-------
while (EndOfFile == false) {
    // Size of the QVIndex descriptor stream in bytes.
    uint64 QVIndexBufferSize;

    // Buffer containing the QVIndex descriptor stream.
    // The data is stored in the ASCII format; thus the QVIndex
    // '0' is actually stored as 0x30 and similarly for the other values.
    uint8[QVIndexBufferSize] QVIndex;
}


Please address all questions to Jan Voges (voges@tnt.uni-hannover.de).

