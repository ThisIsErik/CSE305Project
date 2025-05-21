#ifndef MM_PARALLEL_DATABASE_H
#define MM_PARALLEL_DATABASE_H

#include <vector>
#include <string>
#include <thread>
#include <algorithm>

#include "../../utils/types.h"
#include "../../sequential/mm.h" // for myers_miller_dp

/**
 * @brief Perform parallelized Myers-Miller alignment for a database of sequences.
 *
 * @param query The query sequence
 * @param database Vector of sequences to compare
 * @param match_score Score for a match
 * @param mismatch_penalty Penalty for a mismatch
 * @param gap_penalty Penalty for a gap
 * @param num_threads Number of threads to use
 * @return Vector of MMResult objects
 */
std::vector<MMResult> myers_miller_parallel(
    const std::string& query,
    const std::vector<std::string>& database,
    int match_score = 2,
    int mismatch_penalty = -1,
    int gap_penalty = -1,
    size_t num_threads = 8
);

#endif // MM_PARALLEL_DATABASE_H
