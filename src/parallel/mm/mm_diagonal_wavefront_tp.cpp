#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <string>
#include <iostream>
#include "../../utils/types.h"
#include "thread_pool/thread_pool.h"

void MMAntidiagonalAuxTp(
    const std::string& A,
    const std::string& B,
    int match,
    int mismatch,
    int gap,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& forward_dp,
    std::vector<std::vector<int>>& reverse_dp,
    int mid_a,
    std::pair<int, int>& midpoint
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n)
            continue;

        int score = (A[i - 1] == B[j - 1]) ? match : mismatch;

        forward_dp[i][j] = std::max({
            forward_dp[i - 1][j - 1] + score,
            forward_dp[i][j - 1] + gap,
            forward_dp[i - 1][j] + gap
        });

        // Compute reverse DP matrix
        int rev_i = m - i + 1;
        int rev_j = n - j + 1;
        if (rev_i <= m && rev_j <= n && rev_i > 0 && rev_j > 0) {
            int rev_score = (A[m - rev_i] == B[n - rev_j]) ? match : mismatch;
            reverse_dp[rev_i][rev_j] = std::max({
                reverse_dp[rev_i - 1][rev_j - 1] + rev_score,
                reverse_dp[rev_i][rev_j - 1] + gap,
                reverse_dp[rev_i - 1][rev_j] + gap
            });
        }

        // Check midpoint condition
        if (i == mid_a) {
            midpoint = {i, j};
        }
    }
}

MMResult MyersMillerWavefrontTp(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
) {
    const int m = A.size();
    const int n = B.size();
    int mid_a = m / 2;

    std::vector<std::vector<int>> forward_dp(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> reverse_dp(m + 1, std::vector<int>(n + 1, 0));

    for (int i = 0; i <= m; ++i) {
        forward_dp[i][0] = i * gap;
        reverse_dp[i][0] = i * gap;
    }
    for (int j = 0; j <= n; ++j) {
        forward_dp[0][j] = j * gap;
        reverse_dp[0][j] = j * gap;
    }

    dp::thread_pool pool(num_threads);
    std::pair<int, int> global_midpoint = {mid_a, 0};

    for (int d = 1; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        if (num_cells < THRESHOLD) {
            MMAntidiagonalAuxTp(A, B, match, mismatch, gap, d, i_start, i_end,
                                forward_dp, reverse_dp, mid_a, global_midpoint);
            continue;
        }

        std::mutex midpoint_mutex;
        std::pair<int, int> best_midpoint = global_midpoint;

        int block_size = std::max(1, num_cells / static_cast<int>(num_threads));
        int current_i = i_start;

        for (size_t t = 0; t < num_threads && current_i <= i_end; ++t) {
            int block_end = (t == num_threads - 1 || current_i + block_size - 1 > i_end)
                            ? i_end : (current_i + block_size - 1);

            pool.enqueue_detach([=, &A, &B, &forward_dp, &reverse_dp, &midpoint_mutex, &best_midpoint]() mutable {
                std::pair<int, int> local_mid = {mid_a, 0};
                MMAntidiagonalAuxTp(A, B, match, mismatch, gap, d, current_i, block_end,
                                    forward_dp, reverse_dp, mid_a, local_mid);

                int total = forward_dp[mid_a][local_mid.second] +
                            reverse_dp[m - mid_a][n - local_mid.second + 1];

                std::lock_guard<std::mutex> lock(midpoint_mutex);
                int best = forward_dp[mid_a][best_midpoint.second] +
                           reverse_dp[m - mid_a][n - best_midpoint.second + 1];

                if (total > best) {
                    best_midpoint = local_mid;
                }
            });

            current_i = block_end + 1;
        }

        pool.wait_for_tasks();
        global_midpoint = best_midpoint;
    }

    return {std::move(forward_dp), global_midpoint};
}