#ifndef MM_DIAGONAL_SCORE_ONLY_H
#define MM_DIAGONAL_SCORE_ONLY_H

#include <string>
#include <tuple>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <iostream>
#include <climits>  // Added for INT_MIN
#include "../../utils/types.h"

void MMAntidiagonalAux_ScoreOnly(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    int diagonal,
    int start_i,
    int end_i,
    const std::vector<int>& prev_diag,
    const std::vector<int>& prev_prev_diag,
    int wall_case,
    std::vector<int>& curr_diag,
    std::pair<int, int>& midpoint,
    int mid_a
);

// Returns {score, mid_a, mid_j}
std::tuple<int, int, int> MyersMillerWavefront_ScoreOnly(
    const std::string& A, 
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
);

#endif // MM_DIAGONAL_SCORE_ONLY_H