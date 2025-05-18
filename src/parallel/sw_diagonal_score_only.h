#ifndef SW_DIAGONAL_SCORE_ONLY
#define SW_DIAGONAL_SCORE_ONLY

#include <string>
#include <tuple>
#include <vector>
#include "utils/types.h"

void AntiDiagonalAux_ScoreOnly(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    const std::vector<int>& prev_diag,
    const std::vector<int>& prev_prev_diag,
    int wall_case,
    std::vector<int>& curr_diag,
    LocalMax& local_max
);

std::tuple<int, int, int> SmithWatermanWavefront_ScoreOnly(
    const std::string& A,
    const std::string& B,
    int mismatch_penalty,
    int match_award,
    int gap_penalty,
    size_t num_threads = 1
);

#endif // SW_DIAGONAL_SCORE_ONLY
