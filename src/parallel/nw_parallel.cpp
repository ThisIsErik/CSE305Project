#include <vector>
#include <thread>
#include <string>
#include <algorithm>
#include <iostream>
#include "../utils/types.h"

void ComputeDiagonalChunk(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& dp
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i < 0 || i > m || j < 0 || j > n) continue;

        if (i == 0) {
            dp[i][j] = j * g;
        } else if (j == 0) {
            dp[i][j] = i * g;
        } else {
            int score = (A[i - 1] == B[j - 1]) ? ma : mi;
            dp[i][j] = std::max({
                dp[i - 1][j - 1] + score,
                dp[i - 1][j] + g,
                dp[i][j - 1] + g
            });
        }
    }
}

std::pair<std::string, std::string> Traceback(
    const std::vector<std::vector<int>>& dp,
    const std::string& A,
    const std::string& B,
    int ma, int mi, int g
) {
    std::string alignedA, alignedB;
    int i = A.size();
    int j = B.size();

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int match_score = (A[i - 1] == B[j - 1]) ? ma : mi;
            if (dp[i][j] == dp[i - 1][j - 1] + match_score) {
                alignedA += A[i - 1];
                alignedB += B[j - 1];
                --i; --j;
                continue;
            }
        }
        if (i > 0 && dp[i][j] == dp[i - 1][j] + g) {
            alignedA += A[i - 1];
            alignedB += '-';
            --i;
        } else {
            alignedA += '-';
            alignedB += B[j - 1];
            --j;
        }
    }

    std::reverse(alignedA.begin(), alignedA.end());
    std::reverse(alignedB.begin(), alignedB.end());
    return {alignedA, alignedB};
}

AlignmentResult NeedlemanWunschWavefront(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
) {
    const int m = A.size();
    const int n = B.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::thread> workers(num_threads - 1);

    for (int i = 0; i <= m; ++i) dp[i][0] = i * g;
    for (int j = 0; j <= n; ++j) dp[0][j] = j * g; //global alignment

    for (int d = 2; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end   = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        if (num_cells <= static_cast<int>(num_threads)) {
            ComputeDiagonalChunk(A, B, mi, ma, g, d, i_start, i_end, dp);
            continue;
        }

        int block_size = num_cells / num_threads;
        int current_i = i_start;

        for (size_t t = 0; t < num_threads - 1; ++t) {
            int block_end = current_i + block_size - 1;
            workers[t] = std::thread(
                ComputeDiagonalChunk,
                std::ref(A), std::ref(B),
                mi, ma, g,
                d, current_i, block_end,
                std::ref(dp)
            );
            current_i = block_end + 1;
        }

        ComputeDiagonalChunk(A, B, mi, ma, g, d, current_i, i_end, dp);

        for (auto& t : workers) {
            if (t.joinable()) t.join();
        }
    }

    int final_score = dp[m][n];
    auto aligned = Traceback(dp, A, B, ma, mi, g);
    return {final_score, aligned};
}