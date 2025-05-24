#include "mm_parallel_database.h"
#include "mm_diagonal_wavefront_tp.h" // Include your MM implementation

void MMWorker(
    std::vector<std::string>::const_iterator begin,
    std::vector<std::string>::const_iterator end,
    const std::string& query,
    int match_score, int mismatch_penalty, int gap_penalty,
    std::vector<MMResult>& results,
    size_t result_offset
) {
    for (auto it = begin; it != end; ++it, ++result_offset) {
        // Call your Myers-Miller wavefront implementation instead of myers_miller
        results[result_offset] = MyersMillerWavefrontTp(query, *it, match_score, mismatch_penalty, gap_penalty, 1);
    }
}

std::vector<MMResult> myers_miller_parallel(
    const std::string& query,
    const std::vector<std::string>& database,
    int match_score, int mismatch_penalty, int gap_penalty,
    size_t num_threads
) {
    size_t length = database.size();
    if (length == 0) return {};
    
    num_threads = std::min(num_threads, length);
    size_t block_size = length / num_threads;
    std::vector<MMResult> results(length);
    std::vector<std::thread> workers(num_threads - 1);
    
    auto start_block = database.begin();
    for (size_t i = 0; i < num_threads - 1; ++i) {
        auto end_block = start_block + block_size;
        size_t result_offset = std::distance(database.begin(), start_block);
        
        workers[i] = std::thread(
            MMWorker,
            start_block, end_block,
            std::ref(query),
            match_score, mismatch_penalty, gap_penalty,
            std::ref(results),
            result_offset
        );
        
        start_block = end_block;
    }
    
    size_t result_offset = std::distance(database.begin(), start_block);
    MMWorker(start_block, database.end(), query,
             match_score, mismatch_penalty, gap_penalty,
             results, result_offset);
    
    for (auto& worker : workers) {
        worker.join();
    }
    
    return results;
}