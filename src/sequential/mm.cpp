#include "mm.h"

int match_score(char a, char b, int ma, int mi) {
    return (a == b) ? ma : mi;
}

// Forward DP: scores from A[0..mid] to B[0..j]
std::vector<int> nw_forward(const std::string& A, const std::string& B, int ma, int mi, int g) {
    int m = A.size(), n = B.size();
    std::vector<int> prev(n+1), curr(n+1);
    for (int j = 0; j <= n; ++j) prev[j] = j * g;
    for (int i = 1; i <= m; ++i) {
        curr[0] = i * g;
        for (int j = 1; j <= n; ++j) {
            int score = match_score(A[i-1], B[j-1], ma, mi);
            curr[j] = std::max({
                prev[j-1] + score,
                prev[j] + g,
                curr[j-1] + g
            });
        }
        std::swap(prev, curr);
    }
    return prev;
}

// Reverse DP: scores from A[mid..end] to B[j..end]
std::vector<int> nw_reverse(const std::string& A, const std::string& B, int ma, int mi, int g) {
    int m = A.size(), n = B.size();
    std::vector<int> prev(n+1), curr(n+1);
    for (int j = 0; j <= n; ++j) prev[j] = (n - j) * g;
    for (int i = 1; i <= m; ++i) {
        curr[n] = i * g;
        for (int j = n - 1; j >= 0; --j) {
            int score = match_score(A[m - i], B[j], ma, mi);
            curr[j] = std::max({
                prev[j+1] + score,
                prev[j] + g,
                curr[j+1] + g
            });
        }
        std::swap(prev, curr);
    }
    return prev;
}

// Recursive divide-and-conquer alignment
std::pair<std::string, std::string> myers_miller_align(
    const std::string& A, const std::string& B, int ma, int mi, int g) {

    if (A.empty()) return {std::string(B.size(), '-'), B};
    if (B.empty()) return {A, std::string(A.size(), '-')};
    if (A.size() == 1 || B.size() == 1) {
        // fall back to classic NW in O(mn)
        auto dp = needleman_wunsh_dp(A, B, mi, ma, g);
        return needleman_wunsh_traceback(dp, A, B, ma, mi, g);
    }

    int mid = A.size() / 2;

    std::vector<int> scoreL = nw_forward(A.substr(0, mid), B, ma, mi, g);
    std::vector<int> scoreR = nw_reverse(A.substr(mid), B, ma, mi, g);

    int split = 0;
    int max_score = std::numeric_limits<int>::min();
    for (int j = 0; j <= B.size(); ++j) {
        int s = scoreL[j] + scoreR[j];
        if (s > max_score) {
            max_score = s;
            split = j;
        }
    }

    auto left = myers_miller_align(A.substr(0, mid), B.substr(0, split), ma, mi, g);
    auto right = myers_miller_align(A.substr(mid), B.substr(split), ma, mi, g);

    return {left.first + right.first, left.second + right.second};
}

// Wrapper to return score and alignment
std::pair<int, std::pair<std::string, std::string>> myers_miller(
    const std::string& A, const std::string& B, int ma, int mi, int g) {
    auto result = myers_miller_align(A, B, ma, mi, g);
    int score = 0;
    for (size_t i = 0; i < result.first.size(); ++i) {
        char a = result.first[i];
        char b = result.second[i];
        if (a == '-' || b == '-') {
            score += g;
        } else {
            score += match_score(a, b, ma, mi);
        }
    }
    return {score, result};
}