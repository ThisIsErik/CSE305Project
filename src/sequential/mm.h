#ifndef MM_H
#define MM_H

#include <vector>
#include <string>
#include <utility>
#include "nw.h"

struct AlignmentData {
    std::vector<int> M;   // Match/mismatch scores
    std::vector<int> Ix;  // Gap in A
    std::vector<int> Iy;  // Gap in B
};


struct AlignmentResult {
    std::string seq_a;
    std::string seq_b;
    int score;
};

std::tuple<std::string, std::string, int> myers_miller_align(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend);

#endif // MM_H
