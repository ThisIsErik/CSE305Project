#ifndef MM_DIAGONAL_WAVEFRONT_TP_H
#define MM_DIAGONAL_WAVEFRONT_TP_H

#include <vector>
#include <thread>
#include <algorithm>
#include <string>
#include <iostream>
#include <climits>
#include <future>
#include "../../utils/types.h"
#include "../../utils/thread_pool.h"
//#include "thread_pool/thread_pool.h"

/**
 * @brief Performs Myers-Miller algorithm using wavefront parallelism with thread pool
 * 
 * @param A First sequence
 * @param B Second sequence
 * @param match Score for matching characters
 * @param mismatch Penalty for mismatched characters
 * @param gap Penalty for gaps
 * @param num_threads Number of threads to use
 * @return MMResult containing alignment score, midpoint and matrices
 */
MMResult MyersMillerWavefrontTp(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
);

#endif