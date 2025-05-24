#include "mm_diagonal_score_only.h"

void MMAntidiagonalAux_ScoreOnly(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    int diagonal,
    int start_i,
    int end_i,
    const std::vector<int>& prev_diag,
    const std::vector<int>& prev_prev_diag,
    int wall_case,
    std::vector<int>& curr_diag,
    std::pair<int, int>& midpoint,
    int mid_a
) {
    int m = A.size();
    int n = B.size();

    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n)
            continue;
        
        int diag_idx = i - start_i; 
        int diag_left, diag_above;
        int diag_above_left;

        if (wall_case == 0) { // left wall
            diag_left = diag_idx;
            diag_above = diag_idx - 1;
            diag_above_left = diag_idx - 1;
        } else if (wall_case == 1) { // transition
            diag_left = diag_idx + 1;
            diag_above = diag_idx;
            diag_above_left = diag_idx;
        } else { // right wall
            diag_left = diag_idx + 1;
            diag_above = diag_idx;
            diag_above_left = diag_idx + 1;
        }

        int match_score = (A[i - 1] == B[j - 1]) ? match : mismatch;

        int val = std::max({
            (diag_above_left >= 0 && static_cast<size_t>(diag_above_left) < prev_prev_diag.size() 
                ? prev_prev_diag[diag_above_left] + match_score 
                : match_score),
            (diag_left >= 0 && static_cast<size_t>(diag_left) < prev_diag.size()
                ? prev_diag[diag_left] + gap 
                : gap),
            (diag_above >= 0 && static_cast<size_t>(diag_above) < prev_diag.size()
                ? prev_diag[diag_above] + gap 
                : gap)
        });

        curr_diag[diag_idx] = val;

        // If we are at the middle row update midpoint
        if (i == mid_a) {
            if (midpoint.second == 0 || val > curr_diag[midpoint.second - start_i]) {
                midpoint.second = j;
            }
        }
    }
}

std::tuple<int, int, int> MyersMillerWavefront_ScoreOnly(
    const std::string& A,
    const std::string& B,
    int match, int mismatch, int gap,
    size_t num_threads
) {
    const int m = A.size();
    const int n = B.size();
    int mid_a = m / 2;

    std::pair<int, int> midpoint = {mid_a, 0};
    std::vector<std::pair<int, int>> thread_midpoints(num_threads, {mid_a, 0});
    std::vector<std::thread> workers(num_threads - 1);

    std::vector<int> prev_diag_forward, prev_prev_diag_forward, curr_diag_forward;
    std::vector<int> prev_diag_reverse, prev_prev_diag_reverse, curr_diag_reverse;

    // Forward pass from start to middle
    for (int d = 1; d <= mid_a + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(mid_a, d - 1);
        int len = i_end - i_start + 1;
        
        int wall_case;
        if (i_start == 1) wall_case = 0;
        else if (i_start == 2) wall_case = 1;
        else wall_case = 2;

        prev_prev_diag_forward = std::move(prev_diag_forward);
        prev_diag_forward = std::move(curr_diag_forward);
        curr_diag_forward.assign(len, 0);

        if (len < static_cast<int>(num_threads * 2)) {
            MMAntidiagonalAux_ScoreOnly(A, B, match, mismatch, gap, d, i_start, i_end,
                                   prev_diag_forward, prev_prev_diag_forward,
                                   wall_case,
                                   curr_diag_forward, midpoint, mid_a);
            continue;
        }

        int block_size = len / num_threads;
        int current_i = i_start;

        for (size_t t = 0; t < num_threads - 1; ++t) {
            int block_end = current_i + block_size;
            workers[t] = std::thread(
                MMAntidiagonalAux_ScoreOnly,
                std::ref(A), std::ref(B),
                match, mismatch, gap, d,
                current_i, block_end,
                std::ref(prev_diag_forward),
                std::ref(prev_prev_diag_forward),
                wall_case,
                std::ref(curr_diag_forward),
                std::ref(thread_midpoints[t]),
                mid_a
            );
            current_i = block_end + 1;
        }

        MMAntidiagonalAux_ScoreOnly(A, B, match, mismatch, gap, d, current_i, i_end,
                               prev_diag_forward, prev_prev_diag_forward,
                               wall_case,
                               curr_diag_forward, thread_midpoints[num_threads - 1], mid_a);

        for (auto& t : workers)
            t.join();

        // Combine midpoints
        for (const auto& mp : thread_midpoints) {
            if (mp.second != 0 && (midpoint.second == 0 || 
                curr_diag_forward[mp.second - i_start] > curr_diag_forward[midpoint.second - i_start])) {
                midpoint = mp;
            }
        }
    }

    // Reverse pass from end to middle
    std::string rev_A = A;
    std::string rev_B = B;
    std::reverse(rev_A.begin(), rev_A.end());
    std::reverse(rev_B.begin(), rev_B.end());
    
    for (int d = 1; d <= (m - mid_a) + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m - mid_a, d - 1);
        int len = i_end - i_start + 1;
        
        int wall_case;
        if (i_start == 1) wall_case = 0;
        else if (i_start == 2) wall_case = 1;
        else wall_case = 2;

        prev_prev_diag_reverse = std::move(prev_diag_reverse);
        prev_diag_reverse = std::move(curr_diag_reverse);
        curr_diag_reverse.assign(len, 0);

        if (len < static_cast<int>(num_threads * 2)) {
            MMAntidiagonalAux_ScoreOnly(rev_A, rev_B, match, mismatch, gap, d, i_start, i_end,
                                   prev_diag_reverse, prev_prev_diag_reverse,
                                   wall_case,
                                   curr_diag_reverse, midpoint, m - mid_a);
            continue;
        }

        int block_size = len / num_threads;
        int current_i = i_start;

        for (size_t t = 0; t < num_threads - 1; ++t) {
            int block_end = current_i + block_size;
            workers[t] = std::thread(
                MMAntidiagonalAux_ScoreOnly,
                std::ref(rev_A), std::ref(rev_B),
                match, mismatch, gap, d,
                current_i, block_end,
                std::ref(prev_diag_reverse),
                std::ref(prev_prev_diag_reverse),
                wall_case,
                std::ref(curr_diag_reverse),
                std::ref(thread_midpoints[t]),
                m - mid_a
            );
            current_i = block_end + 1;
        }

        MMAntidiagonalAux_ScoreOnly(rev_A, rev_B, match, mismatch, gap, d, current_i, i_end,
                               prev_diag_reverse, prev_prev_diag_reverse,
                               wall_case,
                               curr_diag_reverse, thread_midpoints[num_threads - 1], m - mid_a);

        for (auto& t : workers)
            t.join();
    }

    // Combine the forward and reverse passes to find the optimal midpoint j
    int best_j = 0;
    int best_score = INT_MIN;

    for (int j = 1; j <= n; ++j) {
        int forward_idx = j - 1;
        int reverse_idx = n - j;
        
        if (forward_idx < 0 || forward_idx >= static_cast<int>(curr_diag_forward.size()) ||
            reverse_idx < 0 || reverse_idx >= static_cast<int>(curr_diag_reverse.size()))
            continue;
            
        int combined_score = curr_diag_forward[forward_idx] + curr_diag_reverse[reverse_idx];
        
        if (combined_score > best_score) {
            best_score = combined_score;
            best_j = j;
        }
    }

    return {best_score, mid_a, best_j};
}