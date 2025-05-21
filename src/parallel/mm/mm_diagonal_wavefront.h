#ifndef MM_DIAGONAL_WAVEFRONT_H
#define MM_DIAGONAL_WAVEFRONT_H

#include <vector>
#include <thread>
#include <algorithm>
#include <string>
#include <iostream>
#include <climits>  // Added for INT_MIN/MAX if needed
#include "../../utils/types.h"

void MMAntidiagonalAux(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& forward_dp,
    std::vector<std::vector<int>>& reverse_dp,
    std::pair<int, int>& midpoint
);

MMResult MyersMillerWavefront(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
);

#endif