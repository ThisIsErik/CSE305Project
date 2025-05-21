#include "mm_diagonal_wavefront_tp.h"
#include <future>
#include <limits>

void MMAntidiagonalAuxTp(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& forward_dp,
    std::vector<std::vector<int>>& reverse_dp,
    int mid_a,
    std::pair<int, int>& local_midpoint
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n) continue;

        // Forward DP (first half)
        if (i <= mid_a) {
            int p = (A[i - 1] == B[j - 1]) ? match : mismatch;
            forward_dp[i][j] = std::max({
                forward_dp[i - 1][j - 1] + p,
                forward_dp[i][j - 1] + gap,
                forward_dp[i - 1][j] + gap
            });
        }

        // Reverse DP (second half)
        if (i >= mid_a) {
            int rev_i = m - i + 1;
            int rev_j = n - j + 1;
            if (rev_i <= 0 || rev_j <= 0) continue;

            int p = (A[m - rev_i] == B[n - rev_j]) ? match : mismatch;
            reverse_dp[rev_i][rev_j] = std::max({
                reverse_dp[rev_i - 1][rev_j - 1] + p,
                reverse_dp[rev_i][rev_j - 1] + gap,
                reverse_dp[rev_i - 1][rev_j] + gap
            });

            // Only compute midpoint if on middle row
            if (i == mid_a) {
                int total = forward_dp[mid_a][j] + reverse_dp[m - mid_a][n - j + 1];
                int current = forward_dp[mid_a][local_midpoint.second] + reverse_dp[m - mid_a][n - local_midpoint.second + 1];
                if (total > current) {
                    local_midpoint = {mid_a, j};
                }
            }
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

    // Initialize DP boundaries
    for (int i = 0; i <= m; ++i) {
        forward_dp[i][0] = i * gap;
        reverse_dp[i][0] = i * gap;
    }
    for (int j = 0; j <= n; ++j) {
        forward_dp[0][j] = j * gap;
        reverse_dp[0][j] = j * gap;
    }

    ThreadPool pool(num_threads);
    std::pair<int, int> global_midpoint = {mid_a, 0};

    for (int d = 1; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        if (num_cells < THRESHOLD) {
            std::pair<int, int> local_mid = {mid_a, 0};
            MMAntidiagonalAuxTp(A, B, match, mismatch, gap, d, i_start, i_end, forward_dp, reverse_dp, mid_a, local_mid);
            int total = forward_dp[mid_a][local_mid.second] + reverse_dp[m - mid_a][n - local_mid.second + 1];
            int best = forward_dp[mid_a][global_midpoint.second] + reverse_dp[m - mid_a][n - global_midpoint.second + 1];
            if (total > best) global_midpoint = local_mid;
            continue;
        }

        std::vector<std::future<std::pair<int, int>>> futures;
        int block_size = std::max(1, num_cells / static_cast<int>(num_threads));
        int current_i = i_start;

        for (size_t t = 0; t < num_threads; ++t) {
            int block_end = (t == num_threads - 1) ? i_end : (current_i + block_size - 1);

            futures.emplace_back(pool.enqueue([=, &A, &B, &forward_dp, &reverse_dp]() mutable {
                std::pair<int, int> local_mid = {mid_a, 0};
                MMAntidiagonalAuxTp(A, B, match, mismatch, gap, d, current_i, block_end, forward_dp, reverse_dp, mid_a, local_mid);
                return local_mid;
            }));

            current_i = block_end + 1;
        }

        for (auto& fut : futures) {
            auto local_mid = fut.get();
            int total = forward_dp[mid_a][local_mid.second] + reverse_dp[m - mid_a][n - local_mid.second + 1];
            int best = forward_dp[mid_a][global_midpoint.second] + reverse_dp[m - mid_a][n - global_midpoint.second + 1];
            if (total > best) global_midpoint = local_mid;
        }
    }

    return {std::move(forward_dp), global_midpoint};
}
