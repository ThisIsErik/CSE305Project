#ifndef MM_DIAGONAL_WAVEFRONT_TP_H
#define MM_DIAGONAL_WAVEFRONT_TP_H

#include <string>
#include <vector>
#include "utils/types.h"
#include "thread_pool/thread_pool.h" 

MMResult MyersMillerWavefrontTp(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
);

#endif // MM_DIAGONAL_WAVEFRONT_TP_H
