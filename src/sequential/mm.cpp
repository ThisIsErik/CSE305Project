#include "mm.h"
#include <algorithm>
#include <limits>
#include <tuple>

const int MIN = std::numeric_limits<int>::min() / 2;

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
    data.M.resize(n + 1);
    data.Iy.resize(n + 1);
    data.Ix.resize(n + 1);

    data.M[0] = 0;
    data.Iy[0] = MIN;
    data.Ix[0] = MIN;

    for (int j = 1; j <= n; ++j) {
        data.M[j] = gap_open + (j - 1) * gap_extend;
        data.Iy[j] = gap_open + (j - 1) * gap_extend;
        data.Ix[j] = MIN;
    }

    std::vector<int> prev_M = data.M;
    std::vector<int> prev_Ix = data.Ix;
    std::vector<int> prev_Iy = data.Iy;

    for (int i = 1; i <= m; ++i) {
        std::vector<int> curr_M(n + 1);
        std::vector<int> curr_Ix(n + 1);
        std::vector<int> curr_Iy(n + 1);

        curr_Ix[0] = gap_open + (i - 1) * gap_extend;
        curr_Iy[0] = MIN;
        curr_M[0] = curr_Ix[0];

        for (int j = 1; j <= n; ++j) {
            int a_idx = start_i + i - 1;
            int b_idx = start_j + j - 1;
            int score_sub = (A[a_idx] == B[b_idx]) ? match_score : mismatch_score;

            curr_M[j] = std::max({prev_M[j - 1], prev_Ix[j - 1], prev_Iy[j - 1]}) + score_sub;
            curr_Ix[j] = std::max(prev_M[j] + gap_open, prev_Ix[j] + gap_extend);
            curr_Iy[j] = std::max(curr_M[j - 1] + gap_open, curr_Iy[j - 1] + gap_extend);
        }

        prev_M = curr_M;
        prev_Ix = curr_Ix;
        prev_Iy = curr_Iy;
    }

    data.M = prev_M;
    data.Iy = prev_Iy;
    data.Ix = prev_Ix;
    return data;
}

AlignmentData reverse_pass(
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
    data.M.resize(n + 1);
    data.Iy.resize(n + 1);
    data.Ix.resize(n + 1);

    data.M[n] = 0;
    data.Iy[n] = MIN;
    data.Ix[n] = MIN;

    for (int j = n - 1; j >= 0; --j) {
        data.M[j] = gap_open + (n - j - 1) * gap_extend;
        data.Iy[j] = gap_open + (n - j - 1) * gap_extend;
        data.Ix[j] = MIN;
    }

    std::vector<int> next_M = data.M;
    std::vector<int> next_Ix = data.Ix;
    std::vector<int> next_Iy = data.Iy;

    for (int i = m - 1; i >= 0; --i) {
        std::vector<int> curr_M(n + 1);
        std::vector<int> curr_Ix(n + 1);
        std::vector<int> curr_Iy(n + 1);

        curr_Ix[n] = gap_open + (m - i - 1) * gap_extend;
        curr_Iy[n] = MIN;
        curr_M[n] = curr_Ix[n];

        for (int j = n - 1; j >= 0; --j) {
            int a_idx = start_i + i;
            int b_idx = start_j + j;
            int score_sub = (A[a_idx] == B[b_idx]) ? match_score : mismatch_score;

            curr_M[j] = std::max({next_M[j + 1], next_Ix[j + 1], next_Iy[j + 1]}) + score_sub;
            curr_Ix[j] = std::max(next_M[j] + gap_open, next_Ix[j] + gap_extend);
            curr_Iy[j] = std::max(curr_M[j + 1] + gap_open, curr_Iy[j + 1] + gap_extend);
        }

        next_M = curr_M;
        next_Ix = curr_Ix;
        next_Iy = curr_Iy;
    }

    data.M = next_M;
    data.Iy = next_Iy;
    data.Ix = next_Ix;
    return data;
}

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

        int score = 0;
        bool in_gap_a = false, in_gap_b = false;
        for (size_t k = 0; k < alignment.first.size(); ++k) {
            char a_char = alignment.first[k];
            char b_char = alignment.second[k];
            if (a_char == '-') {
                score += in_gap_a ? gap_extend : gap_open;
                in_gap_a = true;
                in_gap_b = false;
            } else if (b_char == '-') {
                score += in_gap_b ? gap_extend : gap_open;
                in_gap_b = true;
                in_gap_a = false;
            } else {
                score += (a_char == b_char) ? match_score : mismatch_score;
                in_gap_a = in_gap_b = false;
            }
        }
        return {alignment.first, alignment.second, score};
    }

    int mid_i = start_i + (end_i - start_i) / 2;
    auto forward = forward_pass(A, B, match_score, mismatch_score, gap_open, gap_extend, start_i, mid_i, start_j, end_j);
    auto reverse = reverse_pass(A, B, match_score, mismatch_score, gap_open, gap_extend, mid_i, end_i, start_j, end_j);

    int best_j = start_j;
    int best_score = MIN;
    for (int j = 0; j <= end_j - start_j; ++j) {
        int current_score = forward.M[j] + reverse.M[j];
        if (current_score > best_score) {
            best_score = current_score;
            best_j = start_j + j;
        }
    }

    auto left = align_recursive(A, B, match_score, mismatch_score, gap_open, gap_extend, start_i, mid_i, start_j, best_j);
    auto right = align_recursive(A, B, match_score, mismatch_score, gap_open, gap_extend, mid_i, end_i, best_j, end_j);

    return {
        left.seq_a + right.seq_a,
        left.seq_b + right.seq_b,
        left.score + right.score
    };
}

std::tuple<std::string, std::string, int> myers_miller_align(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend) {

    auto result = align_recursive(A, B, match_score, mismatch_score, gap_open, gap_extend, 0, A.size(), 0, B.size());
    return std::make_tuple(result.seq_a, result.seq_b, result.score);
}
