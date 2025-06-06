#ifndef NW_PARALLEL_H
#define NW_PARALLEL_H

#include "../utils/types.h"
#include <string>
#include <vector>
#include <utility>


AlignmentResult NeedlemanWunschWavefront(
    const std::string& A,
    const std::string& B,
    int mismatch_penalty,
    int match_score,
    int gap_penalty,
    size_t num_threads
);

#endif // PARALLEL_NW_H