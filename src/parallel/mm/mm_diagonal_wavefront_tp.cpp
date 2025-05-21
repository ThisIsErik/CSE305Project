#include "mm_diagonal_wavefront_tp.h"

void MMAntidiagonalAuxTp(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& forward_dp,
    std::vector<std::vector<int>>& reverse_dp,
    LocalMax& local_max
) {
    int m = A.size();
    int n = B.size();
    int mid_a = m / 2;

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n) {
            continue;
        }

        // Forward DP (top half)
        if (i <= mid_a) {
            int p = (A[i - 1] == B[j - 1]) ? match : mismatch;
            forward_dp[i][j] = std::max({
                forward_dp[i - 1][j - 1] + p,
                forward_dp[i][j - 1] + gap,
                forward_dp[i - 1][j] + gap
            });
            // Track local max for forward DP
            if (forward_dp[i][j] > local_max.val || 
                (forward_dp[i][j] == local_max.val && std::make_pair(i, j) > std::make_pair(local_max.i, local_max.j))) {
                local_max.val = forward_dp[i][j];
                local_max.i = i;
                local_max.j = j;
            }
        }

        // Reverse DP (bottom half)
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

            // Optionally track max on reverse DP if desired, but normally max in forward DP is enough
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

    std::vector<std::vector<int>> forward_dp(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> reverse_dp(m + 1, std::vector<int>(n + 1, 0));

    // Initialize first row and column for forward DP
    for (int i = 0; i <= m; ++i) {
        forward_dp[i][0] = i * gap;
    }
    for (int j = 0; j <= n; ++j) {
        forward_dp[0][j] = j * gap;
    }

    // Initialize first row and column for reverse DP
    for (int i = 0; i <= m; ++i) {
        reverse_dp[i][0] = i * gap;
    }
    for (int j = 0; j <= n; ++j) {
        reverse_dp[0][j] = j * gap;
    }

    ThreadPool pool(num_threads);

    LocalMax global_max;
    global_max.val = INT_MIN;

    for (int d = 1; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        if (num_cells < THRESHOLD) {
            LocalMax local;
            MMAntidiagonalAuxTp(A, B, match, mismatch, gap, d, i_start, i_end, forward_dp, reverse_dp, local);
            if (local.val > global_max.val || 
                (local.val == global_max.val && std::make_pair(local.i, local.j) > std::make_pair(global_max.i, global_max.j))) {
                global_max = local;
            }
            continue;
        }

        std::vector<std::future<LocalMax>> futures;
        int block_size = std::max(1, num_cells / static_cast<int>(num_threads));
        int current_i = i_start;

        for (size_t i = 0; i < num_threads; ++i) {
            int block_end = (i == num_threads - 1) ? i_end : (current_i + block_size - 1);

            futures.push_back(pool.enqueue([=, &A, &B, &forward_dp, &reverse_dp]() mutable {
                LocalMax local;
                MMAntidiagonalAuxTp(A, B, match, mismatch, gap, d, current_i, block_end, forward_dp, reverse_dp, local);
                return local;
            }));

            current_i = block_end + 1;
        }

        for (auto& fut : futures) {
            LocalMax local = fut.get();
            if (local.val > global_max.val || 
                (local.val == global_max.val && std::make_pair(local.i, local.j) > std::make_pair(global_max.i, global_max.j))) {
                global_max = local;
            }
        }
    }

    // Return forward DP and max scoring cell
    MMResult result;
    result.first = std::move(forward_dp);
    result.second = {global_max.i, global_max.j};

    return {forward_dp, {global_max.i, global_max.j}};
}
