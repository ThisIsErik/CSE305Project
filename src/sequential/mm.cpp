#include "mm.h"
#include <algorithm>
#include <limits>
#include <tuple>

namespace {

struct AlignmentData {
    std::vector<int> score;
    std::vector<int> gap_open;
    std::vector<int> gap_extend;
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
    
    // Initialize first row
    data.score[0] = 0;
    data.gap_open[0] = std::numeric_limits<int>::min();
    data.gap_extend[0] = std::numeric_limits<int>::min();
    
    for (int j = 1; j <= n; ++j) {
        data.score[j] = gap_open + (j - 1) * gap_extend;
        data.gap_open[j] = std::numeric_limits<int>::min();
        data.gap_extend[j] = data.score[j - 1] + gap_open;
    }
    
    for (int i = 1; i <= m; ++i) {
        int diag_score = data.score[0];
        int diag_gap_open = data.gap_open[0];
        int diag_gap_extend = data.gap_extend[0];
        
        data.score[0] = gap_open + (i - 1) * gap_extend;
        data.gap_open[0] = std::numeric_limits<int>::min();
        data.gap_extend[0] = data.score[0];
        
        for (int j = 1; j <= n; ++j) {
            int match_mismatch = (A[start_i + i - 1] == B[start_j + j - 1]) ? 
                                match_score : mismatch_score;
            
            int new_score = std::max({
                diag_score + match_mismatch,
                diag_gap_open,
                diag_gap_extend
            });
            
            int new_gap_open = std::max(
                data.score[j] + gap_open,
                data.gap_extend[j] + gap_extend
            );
            
            int new_gap_extend = std::max(
                data.score[j - 1] + gap_open,
                data.gap_extend[j - 1] + gap_extend
            );
            
            diag_score = data.score[j];
            diag_gap_open = data.gap_open[j];
            diag_gap_extend = data.gap_extend[j];
            
            data.score[j] = new_score;
            data.gap_open[j] = new_gap_open;
            data.gap_extend[j] = new_gap_extend;
        }
    }
    
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
    data.score.resize(n + 1);
    data.gap_open.resize(n + 1);
    data.gap_extend.resize(n + 1);
    
    // Initialize last row
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
        
        data.score[n] = gap_open + (m - i - 1) * gap_extend;
        data.gap_open[n] = std::numeric_limits<int>::min();
        data.gap_extend[n] = data.score[n];
        
        for (int j = n - 1; j >= 0; --j) {
            int match_mismatch = (A[start_i + i] == B[start_j + j]) ? 
                                match_score : mismatch_score;
            
            int new_score = std::max({
                diag_score + match_mismatch,
                diag_gap_open,
                diag_gap_extend
            });
            
            int new_gap_open = std::max(
                data.score[j] + gap_open,
                data.gap_extend[j] + gap_extend
            );
            
            int new_gap_extend = std::max(
                data.score[j + 1] + gap_open,
                data.gap_extend[j + 1] + gap_extend
            );
            
            diag_score = data.score[j];
            diag_gap_open = data.gap_open[j];
            diag_gap_extend = data.gap_extend[j];
            
            data.score[j] = new_score;
            data.gap_open[j] = new_gap_open;
            data.gap_extend[j] = new_gap_extend;
        }
    }
    
    return data;
}

std::pair<std::string, std::string> align_recursive(
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
        return {a_gaps, b_part};
    }
    
    if (end_j - start_j == 0) {
        std::string a_part = A.substr(start_i, end_i - start_i);
        std::string b_gaps(end_i - start_i, '-');
        return {a_part, b_gaps};
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
        return needleman_wunsh_traceback(
            dp,
            A.substr(start_i, end_i - start_i),
            B.substr(start_j, end_j - start_j),
            match_score,
            mismatch_score,
            gap_open + gap_extend
        );
    }
    
    int mid_i = start_i + (end_i - start_i) / 2;
    
    // Forward pass to middle row
    auto forward = forward_pass(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        start_i, mid_i, start_j, end_j
    );
    
    // Reverse pass from middle row
    auto reverse = reverse_pass(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        mid_i, end_i, start_j, end_j
    );
    
    // Find the best midpoint
    int best_j = start_j;
    int best_score = std::numeric_limits<int>::min();
    for (int j = start_j; j <= end_j; ++j) {
        int current_score = forward.score[j - start_j] + reverse.score[j - start_j];
        if (current_score > best_score) {
            best_score = current_score;
            best_j = j;
        }
    }
    
    // Recursively align the two halves
    auto left = align_recursive(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        start_i, mid_i, start_j, best_j
    );
    
    auto right = align_recursive(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        mid_i, end_i, best_j, end_j
    );
    
    return {
        left.first + right.first,
        left.second + right.second
    };
}

} // namespace

std::pair<std::string, std::string> myers_miller_align(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend) {
    
    return align_recursive(
        A, B, match_score, mismatch_score, gap_open, gap_extend,
        0, A.size(), 0, B.size()
    );
}