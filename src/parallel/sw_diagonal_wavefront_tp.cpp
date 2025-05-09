#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <string>
#include <iostream>
#include "types.h"
#include "thread_pool/thread_pool.h"


void AntiDiagonalAuxTp(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& dp,
    LocalMax& local_max
) {
    int m = A.size();
    int n = B.size();
    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n)
            continue;
        int p = (A[i-1] == B[j-1]) ? ma : mi;
        int val = std::max({
            dp[i-1][j-1] + p,
            dp[i][j-1] + g,
            dp[i-1][j] + g,
            0
        });
        dp[i][j] = val;

        if (val > local_max.val || (val == local_max.val && std::make_pair(i, j) > std::make_pair(local_max.i, local_max.j))) {
            local_max.val = val;
            local_max.i = i;
            local_max.j = j;
        }
    }
}

SWResult SmithWatermanWavefrontTp(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
) {
    const int m = A.size();
    const int n = B.size();

    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    int max_val = 0;
    std::pair<int, int> max_pos = {0, 0};

    std::vector<LocalMax> local_maxes(num_threads); 
    dp::thread_pool pool(num_threads);
    
    //For each antidiagonal
    for (int d = 1; d <= m + n; ++d) {
        
        //Bounds for the curr antidiagonal
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        if (num_cells < THRESHOLD) {
            LocalMax local;
            auto task = [=, &A, &B, &dp, &local] {
                AntiDiagonalAuxTp(A, B, mi, ma, g, d, i_start, i_end, dp, local);
            };
            pool.enqueue_detach(task);
            pool.wait_for_tasks();
            
            if (local.val > max_val || (local.val == max_val && std::make_pair(local.i, local.j) > max_pos)) {
                max_val = local.val;
                max_pos = {local.i, local.j};
            }
            continue;
        }

        // Split work across threads
        int block_size = num_cells / num_threads;
        int current_i = i_start;
        for (size_t i = 0; i < num_threads - 1; ++i) {
            int block_end = current_i + block_size;
            auto task = [=, &A, &B, &dp, &local_maxes] {
                AntiDiagonalAuxTp(A, B, mi, ma, g, d, current_i, block_end, dp, local_maxes[i]);
            };
            pool.enqueue_detach(task);
            current_i = block_end + 1;
        }

        auto task = [=, &A, &B, &dp, &local_maxes] {
                AntiDiagonalAuxTp(A, B, mi, ma, g, d, current_i, i_end, dp, local_maxes[num_threads-1]);
            };
        pool.enqueue_detach(task);
        pool.wait_for_tasks();

        for (const auto& local : local_maxes) {
            if (local.val > max_val || (local.val == max_val && std::make_pair(local.i, local.j) > max_pos)) {
                max_val = local.val;
                max_pos = {local.i, local.j};
            }
        }
    }

    return {dp, max_pos};
}