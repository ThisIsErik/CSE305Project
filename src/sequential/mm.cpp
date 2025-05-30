#include "mm.h"
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <climits>

std::pair<std::string, std::string> myers_miller(
    const std::string& A, 
    const std::string& B, 
    int match, 
    int mismatch, 
    int gap
) {
    return myers_miller_recursive(A, B, 0, A.length(), 0, B.length(), match, mismatch, gap);
}

std::pair<std::string, std::string> myers_miller_recursive(
    const std::string& A, 
    const std::string& B, 
    int start_a, 
    int end_a, 
    int start_b, 
    int end_b, 
    int match, 
    int mismatch, 
    int gap
) {
    int len_a = end_a - start_a;
    int len_b = end_b - start_b;

    if (len_a == 0) {
        return {std::string(len_b, '-'), B.substr(start_b, len_b)};
    } else if (len_b == 0) {
        return {A.substr(start_a, len_a), std::string(len_a, '-')};
    } else if (len_a == 1 || len_b == 1) {
        // Full DP when one sequence is very short
        std::vector<std::vector<int>> dp(len_a + 1, std::vector<int>(len_b + 1, 0));
        for (int i = 0; i <= len_a; ++i) dp[i][0] = i * gap;
        for (int j = 0; j <= len_b; ++j) dp[0][j] = j * gap;

        for (int i = 1; i <= len_a; ++i) {
            for (int j = 1; j <= len_b; ++j) {
                int score = (A[start_a + i - 1] == B[start_b + j - 1]) ? match : mismatch;
                dp[i][j] = std::max({
                    dp[i-1][j-1] + score,
                    dp[i-1][j] + gap,
                    dp[i][j-1] + gap
                });
            }
        }

        std::string alignedA, alignedB;
        int i = len_a, j = len_b;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                int score = (A[start_a + i - 1] == B[start_b + j - 1]) ? match : mismatch;
                if (dp[i][j] == dp[i-1][j-1] + score) {
                    alignedA = A[start_a + i - 1] + alignedA;
                    alignedB = B[start_b + j - 1] + alignedB;
                    --i; --j;
                    continue;
                }
            }

            if (i > 0 && dp[i][j] == dp[i-1][j] + gap) {
                alignedA = A[start_a + i - 1] + alignedA;
                alignedB = '-' + alignedB;
                --i;
            } else {
                alignedA = '-' + alignedA;
                alignedB = B[start_b + j - 1] + alignedB;
                --j;
            }
        }

        return {alignedA, alignedB};
    }

    int mid_a = start_a + len_a / 2;
    int mid_b = find_midpoint(A, B, start_a, mid_a, end_a, start_b, end_b, match, mismatch, gap);

    auto left = myers_miller_recursive(A, B, start_a, mid_a, start_b, mid_b, match, mismatch, gap);
    auto right = myers_miller_recursive(A, B, mid_a, end_a, mid_b, end_b, match, mismatch, gap);

    return {left.first + right.first, left.second + right.second};
}

int find_midpoint(
    const std::string& A, 
    const std::string& B, 
    int start_a, 
    int mid_a, 
    int end_a, 
    int start_b, 
    int end_b, 
    int match, 
    int mismatch, 
    int gap
) {
    std::vector<int> fwd = forward_score(A, B, start_a, mid_a, start_b, end_b, match, mismatch, gap);
    std::vector<int> rev = reverse_score(A, B, mid_a, end_a, start_b, end_b, match, mismatch, gap);

    int best = INT_MIN;
    int mid_b = start_b;

    for (int j = 0; j <= end_b - start_b; ++j) {
        int score = fwd[j] + rev[end_b - start_b - j];
        if (score > best) {
            best = score;
            mid_b = start_b + j;
        }
    }

    return mid_b;
}

std::vector<int> forward_score(
    const std::string& A, 
    const std::string& B, 
    int start_a, 
    int end_a, 
    int start_b, 
    int end_b, 
    int match, 
    int mismatch, 
    int gap
) {
    int len_a = end_a - start_a;
    int len_b = end_b - start_b;

    std::vector<int> current(len_b + 1), previous(len_b + 1);
    for (int j = 0; j <= len_b; ++j) previous[j] = j * gap;

    for (int i = 1; i <= len_a; ++i) {
        current[0] = i * gap;
        for (int j = 1; j <= len_b; ++j) {
            int score = (A[start_a + i - 1] == B[start_b + j - 1]) ? match : mismatch;
            current[j] = std::max({
                previous[j - 1] + score,
                previous[j] + gap,
                current[j - 1] + gap
            });
        }
        std::swap(current, previous);
    }

    return previous;
}

std::vector<int> reverse_score(
    const std::string& A, 
    const std::string& B, 
    int start_a, 
    int end_a, 
    int start_b, 
    int end_b, 
    int match, 
    int mismatch, 
    int gap
) {
    int len_a = end_a - start_a;
    int len_b = end_b - start_b;

    std::vector<int> current(len_b + 1), previous(len_b + 1);
    for (int j = 0; j <= len_b; ++j) previous[j] = j * gap;

    for (int i = 1; i <= len_a; ++i) {
        current[0] = i * gap;
        for (int j = 1; j <= len_b; ++j) {
            int score = (A[end_a - i] == B[end_b - j]) ? match : mismatch;
            current[j] = std::max({
                previous[j - 1] + score,
                previous[j] + gap,
                current[j - 1] + gap
            });
        }
        std::swap(current, previous);
    }

    return previous;
}

// Return DP matrix (for diagnostic purposes) and midpoint (as in SW-style)
std::pair<std::vector<std::vector<int>>, std::pair<int, int>> myers_miller_dp(
    const std::string& A,
    const std::string& B,
    int mismatch,
    int match,
    int gap
) {
    int len_a = A.length();
    int len_b = B.length();

    std::vector<std::vector<int>> dp(len_a + 1, std::vector<int>(len_b + 1));

    for (int i = 0; i <= len_a; ++i) dp[i][0] = i * gap;
    for (int j = 0; j <= len_b; ++j) dp[0][j] = j * gap;

    for (int i = 1; i <= len_a; ++i) {
        for (int j = 1; j <= len_b; ++j) {
            int score = (A[i - 1] == B[j - 1]) ? match : mismatch;
            dp[i][j] = std::max({
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            });
        }
    }

    // Find max score and position
    int max_score = INT_MIN;
    std::pair<int, int> max_pos = {0, 0};
    for (int i = 0; i <= len_a; ++i) {
        for (int j = 0; j <= len_b; ++j) {
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_pos = {i, j};
            }
        }
    }

    return {dp, max_pos};
}

// Traceback using divide-and-conquer (returns aligned strings)
std::pair<std::string, std::string> myers_miller_traceback(
    const std::string& A,
    const std::string& B,
    int mismatch,
    int match,
    int gap
) {
    return myers_miller_recursive(A, B, 0, A.size(), 0, B.size(), match, mismatch, gap);
}
