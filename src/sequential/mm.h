#ifndef MM_H
#define MM_H

#include <string>
#include <vector>
#include <utility>
#include "../utils/types.h"

// Main API-compatible with SW
std::pair<std::vector<std::vector<int>>, std::pair<int, int>> myers_miller_dp(
    const std::string& A,
    const std::string& B,
    int mismatch,
    int match,
    int gap
);

std::pair<std::string, std::string> myers_miller_traceback(
    const std::string& A,
    const std::string& B,
    int mismatch,
    int match,
    int gap
);

// Internals
std::pair<std::string, std::string> myers_miller_recursive(
    const std::string& A,
    const std::string& B,
    int start_a,
    int end_a,
    int start_b,
    int end_b,
    int match,
    int mismatch,
    int gap
);

int find_midpoint(
    const std::string& A,
    const std::string& B,
    int start_a,
    int mid_a,
    int end_a,
    int start_b,
    int end_b,
    int match,
    int mismatch,
    int gap
);

std::vector<int> forward_score(
    const std::string& A,
    const std::string& B,
    int start_a,
    int end_a,
    int start_b,
    int end_b,
    int match,
    int mismatch,
    int gap
);

std::vector<int> reverse_score(
    const std::string& A,
    const std::string& B,
    int start_a,
    int end_a,
    int start_b,
    int end_b,
    int match,
    int mismatch,
    int gap
);

#endif
