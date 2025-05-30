#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <string>
#include <iostream>
#include "utils/types.h"
#include "thread_pool/thread_pool.h"

void AntiDiagonalAux_ScoreOnly_Tp(
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
    LocalMax& local_max,
    int global_i_start
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n)
            continue;
        
        int diag_idx = i - global_i_start; 
        int diag_left, diag_above;
        int diag_above_left;

        if (wall_case == 0) { //left wall
            diag_left = diag_idx;
            diag_above = diag_idx - 1;
            diag_above_left = diag_idx - 1;

        } else if (wall_case == 1) { //transition
            diag_left = diag_idx + 1;
            diag_above = diag_idx;
            diag_above_left = diag_idx;
        }
        else { //right wall
            diag_left = diag_idx + 1;
            diag_above = diag_idx;
            diag_above_left = diag_idx + 1;
        }

        int match_score = (A[i - 1] == B[j - 1]) ? ma : mi;

        int val = std::max({
            (diag_above_left >= 0 && diag_above_left < prev_prev_diag.size() 
                ? prev_prev_diag[diag_above_left] + match_score 
                : match_score),
            (diag_left >= 0 && diag_left < prev_diag.size()
                ? prev_diag[diag_left] + g 
                : g),
            (diag_above >= 0 && diag_above < prev_diag.size()
                ? prev_diag[diag_above] + g 
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

SWResultScore SmithWatermanWavefront_ScoreOnly_Tp(
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
    dp::thread_pool pool(num_threads);

    std::vector<int> prev_diag, prev_prev_diag, curr_diag;

    for (int d = 1; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int len = i_end - i_start + 1;

        int wall_case;
        if (i_start == 1) {
            wall_case = 0;
        } else if (i_start == 2) {
            wall_case = 1;
        } else {
            wall_case = 2;
        }

        prev_prev_diag = std::move(prev_diag);
        prev_diag = std::move(curr_diag);
        curr_diag.assign(len, 0);

        if (len < static_cast<int>(num_threads * 100)) {
            LocalMax local;
            AntiDiagonalAux_ScoreOnly_Tp(A, B, mi, ma, g, d, i_start, i_end,
                                         prev_diag, prev_prev_diag,
                                         wall_case,
                                         curr_diag, local, i_start);
            if (local.val > max_val ||
                (local.val == max_val && std::make_pair(local.i, local.j) > max_pos)) {
                max_val = local.val;
                max_pos = {local.i, local.j};
            }
            continue;
        }

        int block_size = len / num_threads;
        int current_i = i_start;

        for (size_t t = 0; t < num_threads - 1; ++t) {
            int block_end = current_i + block_size;
            auto task = [=, &A, &B, &prev_diag, &prev_prev_diag, &curr_diag, &local_maxes] {
                AntiDiagonalAux_ScoreOnly_Tp(A, B, mi, ma, g, d, current_i, block_end,
                                             prev_diag, prev_prev_diag,
                                             wall_case,
                                             curr_diag, local_maxes[t], i_start);
            };
            pool.enqueue_detach(task);
            current_i = block_end + 1;
        }

        auto task = [=, &A, &B, &prev_diag, &prev_prev_diag, &curr_diag, &local_maxes] {
            AntiDiagonalAux_ScoreOnly_Tp(A, B, mi, ma, g, d, current_i, i_end,
                                         prev_diag, prev_prev_diag,
                                         wall_case,
                                         curr_diag, local_maxes[num_threads - 1], i_start);
        };
        pool.enqueue_detach(task);

        pool.wait_for_tasks();

        for (const auto& local : local_maxes) {
            if (local.val > max_val ||
                (local.val == max_val && std::make_pair(local.i, local.j) > max_pos)) {
                max_val = local.val;
                max_pos = {local.i, local.j};
            }
        }
    }

    return {max_val, max_pos.first, max_pos.second};
}

