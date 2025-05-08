#include <iostream>
#include <thread>
#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <algorithm>

#include "../sequential/sw.h" // smith_waterman_dp
#include "types.h"

void SWWorker(
    std::vector<std::string>::const_iterator begin,
    std::vector<std::string>::const_iterator end,
    const std::string& query,
    int mi, int ma, int g,
    std::vector<SWResult>& results,
    size_t result_offset
) {
    for (auto it = begin; it != end; ++it, ++result_offset) {
        results[result_offset] = smith_waterman_dp(query, *it, mi, ma, g);
    }
}

std::vector<SWResult> smith_waterman_parallel(
    const std::string& query,
    const std::vector<std::string>& database,
    int match_score, int mismatch_penalty, int gap_penalty,
    size_t num_threads
) {
    size_t length = database.size();
    if (length == 0) return {};
    num_threads = std::min(num_threads, length);

    size_t block_size = length / num_threads;

    std::vector<SWResult> results(length);
    std::vector<std::thread> workers(num_threads - 1);

    auto start_block = database.begin();
    for (size_t i = 0; i < num_threads - 1; ++i) {
        auto end_block = start_block + block_size;
        size_t result_offset = std::distance(database.begin(), start_block);
        
        workers[i] = std::thread(
            SWWorker,
            start_block, end_block,
            std::ref(query),
            match_score, mismatch_penalty, gap_penalty,
            std::ref(results),
            result_offset
        );
        
        start_block = end_block;
    }

    size_t result_offset = std::distance(database.begin(), start_block);
    SWWorker(start_block, database.end(), query, 
            match_score, mismatch_penalty, gap_penalty,
            results, result_offset);

    for (auto& worker : workers) {
        worker.join();
    }

    return results;
}