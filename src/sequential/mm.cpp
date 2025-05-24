#include "mm.h"
#include <algorithm>
#include <limits>
#include <vector>
#include <iostream>

static int score(char a, char b, int match, int mismatch) {
    return (a == b) ? match : mismatch;
}

std::vector<int> nw_score(const std::string& A, const std::string& B, int match, int mismatch, int gap) {
    if (B.empty()) {
        return {static_cast<int>(A.size()) * gap};
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
    
    size_t split = 0;
    {
        std::vector<int> left_score = nw_score(A_left, B, match, mismatch, gap);
        
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

int myers_miller_score(
    const std::string& A,
    const std::string& B,
    int match,
    int mismatch,
    int gap) {
    return nw_score(A, B, match, mismatch, gap).back();
}
