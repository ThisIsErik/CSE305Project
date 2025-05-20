#ifndef SW_WF_H
#define SW_WF_H

#include <vector>
#include <string>
#include <utility>
#include <mutex>
#include "utils/types.h"

void AntiDiagonalAux(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    const std::vector<int>& prev_diag,
    const std::vector<int>& prev_prev_diag,
    int prev_i_start,
    int prev_prev_i_start,
    std::vector<int>& curr_diag,
    LocalMax& local_max
);

SWResult SmithWatermanWavefront(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
);

#endif 