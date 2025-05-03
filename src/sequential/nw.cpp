#include "../sequential/nw.h"
#include <vector>
#include <algorithm>
#include <iostream>


//TODO: make it compatible with amino acids (replace p(i,j) with substitution matrix)
std::vector<std::vector<int>> needleman_wunsh_dp(const std::string& A, const std::string& B, int mi, int ma, int g) {
    size_t m = A.size();
    size_t n = B.size();

    // DP table to store the distances of prefixes
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));


    for (size_t i = 0; i <= m; ++i) dp[i][0] = i*g;
    for (size_t j = 0; j <= n; ++j) dp[0][j] = j*g;

    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int p = (A[i-1] == B[j-1]) ? ma : mi;
            dp[i][j] = std::max({
                dp[i - 1][j-1] + p,
                dp[i][j - 1] + g,
                dp[i - 1][j] + g
            });
        }
    }

    return dp;
}



std::pair<std::string, std::string> needleman_wunsh_traceback(const std::vector<std::vector<int>>& dp, const std::string& A, const std::string& B, int ma, int mi, int g) {
    std::string alignedA;
    std::string alignedB;
    size_t i = A.size();
    size_t j = B.size();

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int p = (A[i - 1] == B[j - 1]) ? ma : mi;
            if (dp[i][j] == dp[i - 1][j - 1] + p) {
                alignedA += A[i - 1];
                alignedB += B[j - 1];
                --i;
                --j;
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
