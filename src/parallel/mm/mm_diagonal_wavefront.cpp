#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <string>
#include <iostream>
#include "utils/types.h"
#include "thread_pool/thread_pool.h"
#include "mm_diagonal_wavefront_tp.h"

void FillAntiDiagonalMM(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<int>& prev,
    std::vector<int>& curr
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i < 0 || i > m || j < 0 || j > n)
            continue;

        int score = (A[i - 1] == B[j - 1]) ? ma : mi;

        int diag = (i > 0 && j > 0) ? prev[i - 1] + score : INT_MIN;
        int up   = (i > 0) ? prev[i] + g : INT_MIN;
        int left = (j > 0) ? curr[i - 1] + g : INT_MIN;

        curr[i] = std::max({diag, up, left});
    }
}

std::vector<int> NWScore(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
) {
    int m = A.size();
    int n = B.size();

    std::vector<int> prev(m + 1, 0), curr(m + 1, 0);
    for (int i = 0; i <= m; ++i)
        prev[i] = i * g;

    dp::thread_pool pool(num_threads);

    for (int d = 1; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d);

        int num_cells = i_end - i_start + 1;

        if (num_cells < THRESHOLD) {
            FillAntiDiagonalMM(A, B, mi, ma, g, d, i_start, i_end, prev, curr);
        } else {
            int block_size = num_cells / num_threads;
            int current_i = i_start;

            for (size_t t = 0; t < num_threads; ++t) {
                int block_end = (t == num_threads - 1) ? i_end : (current_i + block_size - 1);

                pool.enqueue_detach([=, &A, &B, &prev, &curr]() {
                    FillAntiDiagonalMM(A, B, mi, ma, g, d, current_i, block_end, prev, curr);
                });

                current_i = block_end + 1;
            }

            pool.wait_for_tasks();
        }

        std::swap(prev, curr);
    }

    return prev;
}

MMResult MyersMillerWavefront(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
) {
    int m = A.size();
    int n = B.size();

    if (m == 0) {
        // Return empty DP matrix and midpoint for base case
        std::vector<std::vector<int>> empty_dp(1, std::vector<int>(n + 1));
        for (int j = 0; j <= n; ++j) empty_dp[0][j] = j * g;
        return {std::move(empty_dp), {0, n}};
    }
    if (n == 0) {
        // Return empty DP matrix and midpoint for base case
        std::vector<std::vector<int>> empty_dp(m + 1, std::vector<int>(1));
        for (int i = 0; i <= m; ++i) empty_dp[i][0] = i * g;
        return {std::move(empty_dp), {m, 0}};
    }

    int mid = n / 2;

    std::string B_left = B.substr(0, mid);
    std::string B_right = B.substr(mid);

    std::vector<int> score_l = NWScore(A, B_left, mi, ma, g, num_threads);
    
    std::string A_rev(A.rbegin(), A.rend());
    std::string B_right_rev(B_right.rbegin(), B_right.rend());
    std::vector<int> score_r = NWScore(A_rev, B_right_rev, mi, ma, g, num_threads);

    int max_score = INT_MIN;
    int split = 0;
    for (int i = 0; i <= m; ++i) {
        int total = score_l[i] + score_r[m - i];
        if (total > max_score) {
            max_score = total;
            split = i;
        }
    }

    // Create a minimal DP matrix for compatibility
    std::vector<std::vector<int>> dp_matrix(m + 1, std::vector<int>(n + 1, 0));
    // You might want to fill this properly if needed elsewhere
    
    return {std::move(dp_matrix), {split, mid}};
}