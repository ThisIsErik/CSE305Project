#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <string>
#include <iostream>
#include "local_max.h"

void AntiDiagonalAux_ScoreOnly(
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
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n)
            continue;

        int diag_idx         = i - start_i;               // for curr_diag
        int diag_idx_prev    = i - prev_i_start;          // for dp[i-1][j] 
        int diag_idx_prev_l  = i - 1 - prev_i_start;      // for dp[i][j-1]
        int diag_idx_prev2   = i - 1 - prev_prev_i_start; // for dp[i-1][j-1]

        int match_score = (A[i - 1] == B[j - 1]) ? ma : mi;

        int val = std::max({
            (diag_idx_prev2 >= 0 && diag_idx_prev2 < static_cast<int>(prev_prev_diag.size()) 
                ? prev_prev_diag[diag_idx_prev2] + match_score 
                : match_score),  
            (diag_idx_prev_l >= 0 && diag_idx_prev_l < static_cast<int>(prev_diag.size()) 
                ? prev_diag[diag_idx_prev_l] + g 
                : g),  
            (diag_idx_prev >= 0 && diag_idx_prev < static_cast<int>(prev_diag.size()) 
                ? prev_diag[diag_idx_prev] + g 
                : g),  
            0
        });

        curr_diag[diag_idx] = val;

        if (val > local_max.val || (val == local_max.val && std::make_pair(i, j) > std::make_pair(local_max.i, local_max.j))) {
            local_max.val = val;
            local_max.i = i;
            local_max.j = j;
        }
    }
}

std::tuple<int, int, int> SmithWatermanWavefront_ScoreOnly(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
) {
    const int m = A.size();
    const int n = B.size();

    int max_val = 0;
    std::pair<int, int> max_pos = {0, 0};

    std::vector<LocalMax> local_maxes(num_threads); 
    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    std::vector<int> prev_diag, prev_prev_diag, curr_diag;

    int prev_i_start = 0;
    int prev_prev_i_start = 0;

    for (int d = 1; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int len = i_end - i_start + 1;

        prev_prev_diag = std::move(prev_diag);
        prev_diag = std::move(curr_diag);
        curr_diag.assign(len, 0);

        prev_prev_i_start = prev_i_start;
        prev_i_start = i_start;

        if (len < static_cast<int>(num_threads * 2)) {
            LocalMax local;
            AntiDiagonalAux_ScoreOnly(A, B, mi, ma, g, d, i_start, i_end,
                                      prev_diag, prev_prev_diag, prev_i_start, prev_prev_i_start,
                                      curr_diag, local);
            if (local.val > max_val) {
                max_val = local.val;
                max_pos = {local.i, local.j};
            }
            continue;
        }

        workers.clear();
        int block_size = len / num_threads;
        int current_i = i_start;

        for (size_t t = 0; t < num_threads - 1; ++t) {
            int block_end = current_i + block_size;
            workers.emplace_back(
                AntiDiagonalAux_ScoreOnly,
                std::ref(A), std::ref(B),
                mi, ma, g, d,
                current_i, block_end,
                std::cref(prev_diag),
                std::cref(prev_prev_diag),
                prev_i_start,
                prev_prev_i_start,
                std::ref(curr_diag),
                std::ref(local_maxes[t])
            );
            current_i = block_end + 1;
        }

        AntiDiagonalAux_ScoreOnly(A, B, mi, ma, g, d, current_i, i_end,
                                  prev_diag, prev_prev_diag,
                                  prev_i_start, prev_prev_i_start,
                                  curr_diag, local_maxes[num_threads - 1]);

        for (auto& t : workers)
            t.join();

        for (const auto& local : local_maxes) {
            if (local.val > max_val || (local.val == max_val && std::make_pair(local.i, local.j) > max_pos)) {
                max_val = local.val;
                max_pos = {local.i, local.j};
            }
        }
    }

    return {max_val, max_pos.first, max_pos.second};
}
