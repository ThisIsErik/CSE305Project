#include "../sequential/sw.h"
#include <vector>
#include <algorithm>
#include <iostream>


//TODO: make it compatible with amino acids (replace p(i,j) with substitution matrix)
std::pair<std::vector<std::vector<int>>, std::pair<int, int>> smith_waterman_dp(const std::string& A, const std::string& B, int mi, int ma, int g) {
    size_t m = A.size();
    size_t n = B.size();

    size_t maxi, maxj;
    int maxval = 0;

    // DP table to store the distances of prefixes
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));


    for (size_t i = 0; i <= m; ++i) dp[i][0] = 0;
    for (size_t j = 0; j <= n; ++j) dp[0][j] = 0;

    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int p = (A[i-1] == B[j-1]) ? ma : mi;
            dp[i][j] = std::max({
                dp[i - 1][j-1] + p,
                dp[i][j - 1] + g,
                dp[i - 1][j] + g,
                0
            });
            if (dp[i][j] > maxval || (dp[i][j] == maxval && std::make_pair(i, j) > std::make_pair(maxi, maxj))) {
                maxval = dp[i][j];
                maxi = i;
                maxj = j;
            }
        }
    }

    return {dp, {maxi, maxj}};
}


std::pair<std::string, std::string> smith_waterman_traceback(
    const std::vector<std::vector<int>>& dp,
    const std::string& A,
    const std::string& B,
    int ma, int mi, int g,
    std::pair<int, int> start) {

    std::string alignedA, alignedB;
    int i = start.first, j = start.second;

    while (i > 0 && j > 0 && dp[i][j] > 0) {
        int p = (A[i - 1] == B[j - 1]) ? ma : mi;

        if (dp[i][j] == dp[i - 1][j - 1] + p) {
            alignedA += A[i - 1];
            alignedB += B[j - 1];
            --i; --j;
        } else if (dp[i][j] == dp[i - 1][j] + g) {
            alignedA += A[i - 1];
            alignedB += '-';
            --i;
        } else if (dp[i][j] == dp[i][j - 1] + g) {
            alignedA += '-';
            alignedB += B[j - 1];
            --j;
        } else {
            break;  // score is zero
        }
    }

    std::reverse(alignedA.begin(), alignedA.end());
    std::reverse(alignedB.begin(), alignedB.end());
    return {alignedA, alignedB};
}