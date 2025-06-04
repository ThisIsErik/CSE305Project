#include "mm.h"
#include <iostream>
#include <climits>

std::pair<std::vector<int>, std::vector<int>> MyersMillerAligner::forward_pass(
    const std::string& A, const std::string& B, 
    
    int rows = end_i - start_i + 1;
    int cols = end_j - start_j + 1;
    
    std::vector<int> curr(cols, INT_MIN);
    std::vector<int> curr_gap(cols, INT_MIN);
    std::vector<int> prev(cols, INT_MIN);
    
    curr[0] = 0; //init the first row
    for (int j = 1; j < cols; j++) {
        curr[j] = curr[j-1] + gap_score;
        curr_gap[j] = curr[j];
    }
    
    for (int i = 1; i < rows; i++) { //filling the dp table
        prev = curr;
        curr[0] = prev[0] + gap_score;
        curr_gap[0] = curr[0];
        
        for (int j = 1; j < cols; j++) {
            char a_char = A[start_i + i - 1];
            char b_char = B[start_j + j - 1];
            int score = (a_char == b_char) ? match_score : mismatch_score;

            int diagonal = prev[j-1] + score; //get reg score
            int left = curr[j-1] + gap_score;
            int up = prev[j] + gap_score;
            
            curr[j] = std::max({diagonal, left, up});
            
            curr_gap[j] = std::max(curr_gap[j-1] + gap_score, curr[j-1] + gap_score); //get gap extension score
        }
    }
    
    return {curr, curr_gap};
}

std::pair<std::vector<int>, std::vector<int>> MyersMillerAligner::backward_pass(
    const std::string& A, const std::string& B, 
    int start_i, int end_i, int start_j, int end_j) {
    
    int rows = end_i - start_i + 1;
    int cols = end_j - start_j + 1;
    
    std::vector<int> curr(cols, INT_MIN);
    std::vector<int> curr_gap(cols, INT_MIN);
    std::vector<int> prev(cols, INT_MIN); //same as before
    
    curr[cols-1] = 0; //now we initialize the last row and work backwards
    for (int j = cols-2; j >= 0; j--) {
        curr[j] = curr[j+1] + gap_score;
        curr_gap[j] = curr[j];
    }
    
    for (int i = rows-2; i >= 0; i--) {
        prev = curr;
        curr[cols-1] = prev[cols-1] + gap_score;
        curr_gap[cols-1] = curr[cols-1];
        
        for (int j = cols-2; j >= 0; j--) {
            char a_char = A[start_i + i];
            char b_char = B[start_j + j];
            int score = (a_char == b_char) ? match_score : mismatch_score;

            int diagonal = prev[j+1] + score;
            int right = curr[j+1] + gap_score;
            int down = prev[j] + gap_score;
            
            curr[j] = std::max({diagonal, right, down});

            curr_gap[j] = std::max(curr_gap[j+1] + gap_score, curr[j+1] + gap_score);
        }
    }
    
    return {curr, curr_gap};
}

std::pair<int, bool> MyersMillerAligner::find_midpoint(
    const std::vector<int>& forward_scores,
    const std::vector<int>& forward_gap_scores,
    const std::vector<int>& backward_scores,
    const std::vector<int>& backward_gap_scores,
    int cols) {
    
    int max_score = INT_MIN;
    int best_j = 0;
    bool use_gap = false;
    
    for (int j = 0; j < cols; j++) {
        int score1 = forward_scores[j] + backward_scores[j]; //type 1 is reg alignment
        int score2 = forward_gap_scores[j] + backward_gap_scores[j]; //type 2 gap alginemnt
        
        if (score1 > max_score) {
            max_score = score1;
            best_j = j;
            use_gap = false;
        }
        
        if (score2 > max_score) {
            max_score = score2;
            best_j = j;
            use_gap = true;
        }
    }
    
    return {best_j, use_gap};
}

std::vector<AlignOp> MyersMillerAligner::align_recursive(
    const std::string& A, const std::string& B,
    int start_i, int end_i, int start_j, int end_j) {
    
    std::vector<AlignOp> result;

    if (start_i >= end_i && start_j >= end_j) { //base case
        return result; //empty alignment
    }
    
    if (start_i >= end_i) { //only insertions
        for (int j = start_j; j < end_j; j++) {
            result.push_back(AlignOp(AlignOp::INSERT, 0, B[j]));
        }
        return result;
    }
    
    if (start_j >= end_j) { //only deletions
        for (int i = start_i; i < end_i; i++) {
            result.push_back(AlignOp(AlignOp::DELETE, A[i], 0));
        }
        return result;
    }
    
    if (end_i - start_i == 1) {
        char a_char = A[start_i];
        
        int best_score = INT_MIN; //looking for best alignment is A only has one char
        int best_j = start_j;
        bool best_is_match = false;
        
        //either delete A[start_i] and insert all of B
        int delete_score = gap_score + (end_j - start_j) * gap_score;
        if (delete_score > best_score) {
            best_score = delete_score;
            best_j = start_j;
            best_is_match = false;
        }
        
        //either try aligning A[start_i] with each B[j]
        for (int j = start_j; j < end_j; j++) {
            char b_char = B[j];
            int align_score = (a_char == b_char) ? match_score : mismatch_score;
            align_score += (j - start_j) * gap_score; //insertions before
            align_score += (end_j - j - 1) * gap_score; //insertions after
            
            if (align_score > best_score) {
                best_score = align_score;
                best_j = j;
                best_is_match = true;
            }
        }
        if (best_is_match) {
            for (int j = start_j; j < best_j; j++) {
                result.push_back(AlignOp(AlignOp::INSERT, 0, B[j]));
            }
            if (a_char == B[best_j]) { //match or mismatch
                result.push_back(AlignOp(AlignOp::MATCH, a_char, B[best_j]));
            } else {
                result.push_back(AlignOp(AlignOp::MISMATCH, a_char, B[best_j]));
            }
            for (int j = best_j + 1; j < end_j; j++) {
                result.push_back(AlignOp(AlignOp::INSERT, 0, B[j]));
            }
        } else {
            result.push_back(AlignOp(AlignOp::DELETE, a_char, 0));
            for (int j = start_j; j < end_j; j++) {
                result.push_back(AlignOp(AlignOp::INSERT, 0, B[j]));
            }
        }
        
        return result;
    }
    int mid_i = start_i + (end_i - start_i) / 2;
    auto forward_result = forward_pass(A, B, start_i, mid_i, start_j, end_j);
    auto backward_result = backward_pass(A, B, mid_i, end_i, start_j, end_j);
    
    auto midpoint = find_midpoint(forward_result.first, forward_result.second,
                                 backward_result.first, backward_result.second,
                                 end_j - start_j + 1);
    
    int mid_j = start_j + midpoint.first;
    
    auto left_ops = align_recursive(A, B, start_i, mid_i, start_j, mid_j);
    auto right_ops = align_recursive(A, B, mid_i, end_i, mid_j, end_j);
    
    result.insert(result.end(), left_ops.begin(), left_ops.end());
    result.insert(result.end(), right_ops.begin(), right_ops.end());
    
    return result;
}

std::pair<std::string, std::string> MyersMillerAligner::ops_to_strings(const std::vector<AlignOp>& ops) {
    std::string aligned_A, aligned_B;
    
    for (const auto& op : ops) {
        switch (op.type) {
            case AlignOp::MATCH:
            case AlignOp::MISMATCH:
                aligned_A += op.a_char;
                aligned_B += op.b_char;
                break;
            case AlignOp::DELETE:
                aligned_A += op.a_char;
                aligned_B += '-';
                break;
            case AlignOp::INSERT:
                aligned_A += '-';
                aligned_B += op.b_char;
                break;
        }
    }
    
    return {aligned_A, aligned_B};
}

int MyersMillerAligner::calculate_score(const std::vector<AlignOp>& ops) {
    int score = 0;
    
    for (const auto& op : ops) {
        switch (op.type) {
            case AlignOp::MATCH:
                score += match_score;
                break;
            case AlignOp::MISMATCH:
                score += mismatch_score;
                break;
            case AlignOp::DELETE:
            case AlignOp::INSERT:
                score += gap_score;
                break;
        }
    }
    
    return score;
}

MMResult MyersMillerAligner::align(const std::string& A, const std::string& B) {
    if (A.empty() && B.empty()) {
        return {{"", ""}, 0};
    }
    
    if (A.empty()) {
        std::string gaps(B.length(), '-');
        return {{gaps, B}, static_cast<int>(B.length()) * gap_score};
    }
    
    if (B.empty()) {
        std::string gaps(A.length(), '-');
        return {{A, gaps}, static_cast<int>(A.length()) * gap_score};
    }
    
    auto ops = align_recursive(A, B, 0, A.length(), 0, B.length());
    auto aligned_strings = ops_to_strings(ops);
    int score = calculate_score(ops);
    
    return {aligned_strings, score};
}

MMResult myers_miller_align(const std::string& A, const std::string& B, 
                           int match, int mismatch, int gap) {
    MyersMillerAligner aligner(match, mismatch, gap);
    return aligner.align(A, B);
} 