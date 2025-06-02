#include "mm.h"
#include <algorithm>
#include <limits>
#include <vector>
#include <iostream>

<<<<<<< Updated upstream
static int score(char a, char b, int match, int mismatch) {
    return (a == b) ? match : mismatch;
}

std::vector<int> nw_score(const std::string& A, const std::string& B, int match, int mismatch, int gap) {
    if (B.empty()) {
        return {static_cast<int>(A.size()) * gap};
=======

struct AlignmentData {
    std::vector<int> score;
    std::vector<int> gap_open;
    std::vector<int> gap_extend;
};

struct AlignmentResult {
    std::string seq_a;
    std::string seq_b;
    int score;
};

AlignmentData forward_pass(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend,
    int start_i,
    int end_i,
    int start_j,
    int end_j) {
    
    int m = end_i - start_i;
    int n = end_j - start_j;
    
    AlignmentData data;
    data.score.resize(n + 1);
    data.gap_open.resize(n + 1);
    data.gap_extend.resize(n + 1);
    
    data.score[0] = 0;
    data.gap_open[0] = std::numeric_limits<int>::min();
    data.gap_extend[0] = std::numeric_limits<int>::min();
    
    for (int j = 1; j <= n; ++j) {
        data.score[j] = gap_open + (j - 1) * gap_extend;
        data.gap_open[j] = std::numeric_limits<int>::min();
        data.gap_extend[j] = data.score[j - 1] + gap_open;
>>>>>>> Stashed changes
    }
    
    std::vector<int> prev(B.size() + 1), curr(B.size() + 1);
    
    for (size_t j = 0; j <= B.size(); ++j) {
        prev[j] = j * gap;
    }
    
    for (size_t i = 1; i <= A.size(); ++i) {
        curr[0] = i * gap;
        for (size_t j = 1; j <= B.size(); ++j) {
            curr[j] = std::max({
                prev[j - 1] + score(A[i - 1], B[j - 1], match, mismatch),
                prev[j] + gap,
                curr[j - 1] + gap
            });
        }
        std::swap(prev, curr);
    }
    
    return prev;
}

std::pair<std::string, std::string> myers_miller(
    const std::string& A,
    const std::string& B,
    int match,
    int mismatch,
    int gap) {
    
    if (A.empty() && B.empty()) {
        return {"", ""};
    }
    if (A.empty()) {
        return {std::string(B.size(), '-'), B};
    }
    if (B.empty()) {
        return {A, std::string(A.size(), '-')};
    }
    
    if (A.size() <= 1 || B.size() <= 1) {
        if (A.size() == 1 && B.size() == 1) {
            if (A[0] == B[0]) {
                return {A, B};
            } else {
                return {A, B};
            }
        } else if (A.size() == 1) {
            return {A + std::string(B.size() - 1, '-'), B};
        } else {
            return {A, std::string(A.size() - 1, '-') + B};
        }
    }
    
    size_t mid = A.size() / 2;
    std::string A_left = A.substr(0, mid);
    std::string A_right = A.substr(mid);
    
<<<<<<< Updated upstream
    size_t split = 0;
    {
        std::vector<int> left_score = nw_score(A_left, B, match, mismatch, gap);
=======
    int m = end_i - start_i;
    int n = end_j - start_j;
    
    AlignmentData data;
    data.score.resize(n + 1);
    data.gap_open.resize(n + 1);
    data.gap_extend.resize(n + 1);
    
    data.score[n] = 0;
    data.gap_open[n] = std::numeric_limits<int>::min();
    data.gap_extend[n] = std::numeric_limits<int>::min();
    
    for (int j = n - 1; j >= 0; --j) {
        data.score[j] = gap_open + (n - j - 1) * gap_extend;
        data.gap_open[j] = std::numeric_limits<int>::min();
        data.gap_extend[j] = data.score[j + 1] + gap_open;
    }
    
    for (int i = m - 1; i >= 0; --i) {
        int diag_score = data.score[n];
        int diag_gap_open = data.gap_open[n];
        int diag_gap_extend = data.gap_extend[n];
>>>>>>> Stashed changes
        
        std::string A_right_rev = A_right;
        std::string B_rev = B;
        std::reverse(A_right_rev.begin(), A_right_rev.end());
        std::reverse(B_rev.begin(), B_rev.end());
        
        std::vector<int> right_score = nw_score(A_right_rev, B_rev, match, mismatch, gap);
        
        int max_score = std::numeric_limits<int>::min();
        for (size_t j = 0; j <= B.size(); ++j) {
            // Fixed boundary checks
            if (j >= left_score.size() || (B.size() - j) > right_score.size() - 1) {
                continue;
            }
            
            int total_score = left_score[j] + right_score[B.size() - j];
            if (total_score > max_score) {
                max_score = total_score;
                split = j;
            }
        }
        
        if (split > B.size()) split = B.size();
    }
    
    auto left = myers_miller(A_left, B.substr(0, split), match, mismatch, gap);
    auto right = myers_miller(A_right, B.substr(split), match, mismatch, gap);
    
    return {left.first + right.first, left.second + right.second};
}

<<<<<<< Updated upstream
int myers_miller_score(
    const std::string& A,
    const std::string& B,
    int match,
    int mismatch,
    int gap) {
    return nw_score(A, B, match, mismatch, gap).back();
}
=======
AlignmentResult align_recursive(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend,
    int start_i,
    int end_i,
    int start_j,
    int end_j) {
    
    if (end_i - start_i == 0) {
        std::string a_gaps(end_j - start_j, '-');
        std::string b_part = B.substr(start_j, end_j - start_j);
        int score = (end_j - start_j > 0) ? gap_open + (end_j - start_j - 1) * gap_extend : 0;
        return {a_gaps, b_part, score};
    }
    
    if (end_j - start_j == 0) {
        std::string a_part = A.substr(start_i, end_i - start_i);
        std::string b_gaps(end_i - start_i, '-');
        int score = (end_i - start_i > 0) ? gap_open + (end_i - start_i - 1) * gap_extend : 0;
        return {a_part, b_gaps, score};
    }
    
    if (end_i - start_i == 1 || end_j - start_j == 1) {
        // Base case: use Needleman-Wunsch for small sequences
        auto dp = needleman_wunsh_dp(
            A.substr(start_i, end_i - start_i),
            B.substr(start_j, end_j - start_j),
            mismatch_score,
            match_score,
            gap_open + gap_extend
        );
        auto alignment = needleman_wunsh_traceback(
            dp,
            A.substr(start_i, end_i - start_i),
            B.substr(start_j, end_j - start_j),
            match_score,
            mismatch_score,
            gap_open + gap_extend
        );
        
        //compute the score for base case alignment
        int score = 0;
        bool in_gap_a = false, in_gap_b = false;
        for (size_t k = 0; k < alignment.first.size(); ++k) {
            char a_char = alignment.first[k];
            char b_char = alignment.second[k];
            
            if (a_char == '-') {
                //gap in sequence A
                if (!in_gap_a) {
                    score += gap_open;
                    in_gap_a = true;
                } else {
                    score += gap_extend;
                }
                in_gap_b = false;
            } else if (b_char == '-') {
                //Gap in sequence B
                if (!in_gap_b) {
                    score += gap_open;
                    in_gap_b = true;
                } else {
                    score += gap_extend;
                }
                in_gap_a = false;
            } else {
                //is a match or mismatch?
                score += (a_char == b_char) ? match_score : mismatch_score;
                in_gap_a = in_gap_b = false;
            }
        }
        
        return {alignment.first, alignment.second, score};
    }
    
    int mid_i = start_i + (end_i - start_i) / 2;
    
    auto forward = forward_pass(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        start_i, mid_i, start_j, end_j
    );
    
    auto reverse = reverse_pass(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        mid_i, end_i, start_j, end_j
    );
    
    int best_j = start_j;
    int best_score = std::numeric_limits<int>::min();
    for (int j = start_j; j <= end_j; ++j) {
        int current_score = forward.score[j - start_j] + reverse.score[j - start_j];
        if (current_score > best_score) {
            best_score = current_score;
            best_j = j;
        }
    }
    
    auto left = align_recursive(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        start_i, mid_i, start_j, best_j
    );
    
    auto right = align_recursive(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        mid_i, end_i, best_j, end_j
    );
    
    return {
        left.seq_a + right.seq_a,
        left.seq_b + right.seq_b,
        left.score + right.score
    };
}

std::tuple<std::string, std::string, int>  myers_miller_align(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend) {
    
    auto result = align_recursive(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        0, A.size(), 0, B.size()
    );
    
    return std::make_tuple(result.seq_a, result.seq_b, result.score);
}
>>>>>>> Stashed changes
