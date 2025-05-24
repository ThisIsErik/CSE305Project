#ifndef MM_DIAGONAL_WAVEFRONT_TP_H
#define MM_DIAGONAL_WAVEFRONT_TP_H

#include <string>
#include "utils/types.h"

// mm_diagonal_wavefront_tp.h
MMResult MyersMillerWavefrontTp(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads);


#endif // MM_DIAGONAL_WAVEFRONT_TP_H

