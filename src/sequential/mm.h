#ifndef MM_H
#define MM_H

#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <climits>
#include "../utils/types.h"

struct AlignOp {
    enum Type { MATCH, MISMATCH, INSERT, DELETE };
    Type type;
    char a_char;
    char b_char;
    
    AlignOp(Type t, char a = 0, char b = 0) : type(t), a_char(a), b_char(b) {}
};

class MyersMillerAligner {
private:
    int match_score;
    int mismatch_score;
    int gap_score;

    std::pair<std::vector<int>, std::vector<int>> forward_pass(
        const std::string& A, const std::string& B, 
        int start_i, int end_i, int start_j, int end_j);

    std::pair<std::vector<int>, std::vector<int>> backward_pass(
        const std::string& A, const std::string& B, 
        int start_i, int end_i, int start_j, int end_j);

    std::pair<int, bool> find_midpoint(
        const std::vector<int>& forward_scores,
        const std::vector<int>& forward_gap_scores,
        const std::vector<int>& backward_scores,
        const std::vector<int>& backward_gap_scores,
        int cols);

    std::vector<AlignOp> align_recursive(
        const std::string& A, const std::string& B,
        int start_i, int end_i, int start_j, int end_j);

    std::pair<std::string, std::string> ops_to_strings(const std::vector<AlignOp>& ops);

    int calculate_score(const std::vector<AlignOp>& ops);

public:
    MyersMillerAligner(int match = 1, int mismatch = -1, int gap = -2)
        : match_score(match), mismatch_score(mismatch), gap_score(gap) {}

    MMResult align(const std::string& A, const std::string& B);
};

MMResult myers_miller_align(const std::string& A, const std::string& B, 
                           int match, int mismatch, int gap);

#endif // MM_H
