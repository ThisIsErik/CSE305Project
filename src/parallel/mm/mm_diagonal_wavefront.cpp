#include "mm_diagonal_wavefront.h"
#include <thread>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

void MMAntidiagonalAux(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& forward_dp,
    std::vector<std::vector<int>>& reverse_dp,
    std::pair<int, int>& midpoint
) {
    int m = A.size();
    int n = B.size();
    int mid_a = m / 2;

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i < 0 || i > m || j < 0 || j > n) continue;

        // Forward pass
        if (i <= mid_a && i > 0 && j > 0) {
            int p = (A[i - 1] == B[j - 1]) ? match : mismatch;
            forward_dp[i][j] = std::max({
                forward_dp[i - 1][j - 1] + p,
                forward_dp[i][j - 1] + gap,
                forward_dp[i - 1][j] + gap
            });
        }

        // Reverse pass
        if (i >= mid_a && i < m && j < n && i > 0 && j > 0) {
            int rev_i = m - i;
            int rev_j = n - j;
            int p = (A[rev_i] == B[rev_j]) ? match : mismatch;
            reverse_dp[rev_i][rev_j] = std::max({
                reverse_dp[rev_i + 1][rev_j + 1] + p,
                reverse_dp[rev_i][rev_j + 1] + gap,
                reverse_dp[rev_i + 1][rev_j] + gap
            });

            if (i == mid_a) {
                int combined = forward_dp[mid_a][j] + reverse_dp[m - mid_a][n - j];
                int current_best = forward_dp[mid_a][midpoint.second] + reverse_dp[m - mid_a][n - midpoint.second];
                if (combined > current_best) {
                    midpoint = {mid_a, j};
                }
            }
        }
    }
}

MMResult MyersMillerWavefront(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
) {
    int m = A.size(), n = B.size();
    int mid_a = m / 2;

    std::vector<std::vector<int>> forward_dp(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> reverse_dp(m + 1, std::vector<int>(n + 1, 0));

    // Init forward matrix
    for (int i = 0; i <= m; ++i) forward_dp[i][0] = i * gap;
    for (int j = 0; j <= n; ++j) forward_dp[0][j] = j * gap;

    // Init reverse matrix (from m,n to mid)
    for (int i = m; i >= 0; --i) reverse_dp[i][n] = (m - i) * gap;
    for (int j = n; j >= 0; --j) reverse_dp[m][j] = (n - j) * gap;

    std::pair<int, int> midpoint = {mid_a, 0};
    std::vector<std::thread> workers(num_threads - 1);

    for (int d = 2; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        std::vector<std::pair<int, int>> thread_midpoints(num_threads, midpoint);

        if (num_cells < THRESHOLD) {
            MMAntidiagonalAux(A, B, match, mismatch, gap, d, i_start, i_end,
                              forward_dp, reverse_dp, midpoint);
            continue;
        }

        int block_size = num_cells / num_threads;
        int current_i = i_start;

        for (size_t i = 0; i < num_threads - 1; ++i) {
            int block_end = current_i + block_size;
            workers[i] = std::thread(
                MMAntidiagonalAux,
                std::ref(A), std::ref(B),
                match, mismatch, gap,
                d, current_i, block_end,
                std::ref(forward_dp), std::ref(reverse_dp),
                std::ref(thread_midpoints[i])
            );
            current_i = block_end + 1;
        }

        MMAntidiagonalAux(A, B, match, mismatch, gap, d, current_i, i_end,
                          forward_dp, reverse_dp, thread_midpoints[num_threads - 1]);

        for (auto& t : workers) t.join();

        for (const auto& mp : thread_midpoints) {
            int combined = forward_dp[mid_a][mp.second] + reverse_dp[m - mid_a][n - mp.second];
            int best = forward_dp[mid_a][midpoint.second] + reverse_dp[m - mid_a][n - midpoint.second];
            if (combined > best) {
                midpoint = mp;
            }
        }
    }

    return {std::move(forward_dp), midpoint};
}
