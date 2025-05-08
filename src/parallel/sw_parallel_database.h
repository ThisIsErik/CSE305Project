#ifndef SW_P_H
#define SW_P_H

#include <vector>
#include <string>
#include <utility>
#include "local_max.h"

std::vector<SWResult> smith_waterman_parallel(
    const std::string& query,
    const std::vector<std::string>& database,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = -1,
    size_t num_threads = 8
);

#endif 