#include "nw_parallel.h"
#include <thread>
#include <algorithm>
#include <future>
#include <mutex>

std::vector<std::vector<std::pair<int, int>>> get_antidiagonal_indices(int m, int n) {
    std::vector<std::vector<std::pair<int, int>>> antidiagonals;
    int total_diagonals = m + n + 1;
    
    for (int d = 0; d < total_diagonals; ++d) {
        std::vector<std::pair<int, int>> diagonal;
        
        //for diagonal d we need i + j = d
        //constraints: 0 <= i <= m, 0 <= j <= n
        int start_i = std::max(0, d - n);
        int end_i = std::min(m, d);
        
        for (int i = start_i; i <= end_i; ++i) {
            int j = d - i;
            if (j >= 0 && j <= n) {
                diagonal.push_back({i, j});
            }
        }
        
        if (!diagonal.empty()) {
            antidiagonals.push_back(diagonal);
        }
    }
    
    return antidiagonals;
}

void compute_diagonal_parallel(std::vector<std::vector<int>>& dp,
                             const std::vector<std::pair<int, int>>& diagonal,
                             const std::string& A, const std::string& B,
                             int mi, int ma, int g) {
    
    const int chunk_size = 10;
    const int num_threads = std::thread::hardware_concurrency();
    
    if (diagonal.size() <= chunk_size || num_threads <= 1) {
        for (const auto& cell : diagonal) {
            int i = cell.first;
            int j = cell.second;
            
            if (i == 0 && j == 0) {
                dp[i][j] = 0;
            } else if (i == 0) {
                dp[i][j] = j * g;
            } else if (j == 0) {
                dp[i][j] = i * g;
            } else {
                int match_score = (A[i-1] == B[j-1]) ? ma : mi;
                dp[i][j] = std::max({
                    dp[i-1][j-1] + match_score,
                    dp[i-1][j] + g,
                    dp[i][j-1] + g
                });
            }
        }
        return;
    }

    std::vector<std::future<void>> futures;
    int cells_per_chunk = std::max(1, static_cast<int>(diagonal.size()) / num_threads);
    
    for (size_t start = 0; start < diagonal.size(); start += cells_per_chunk) {
        size_t end = std::min(start + cells_per_chunk, diagonal.size());
        
        futures.push_back(std::async(std::launch::async, [&, start, end]() {
            for (size_t idx = start; idx < end; ++idx) {
                int i = diagonal[idx].first;
                int j = diagonal[idx].second;
                
                if (i == 0 && j == 0) {
                    dp[i][j] = 0;
                } else if (i == 0) {
                    dp[i][j] = j * g;
                } else if (j == 0) {
                    dp[i][j] = i * g;
                } else {
                    int match_score = (A[i-1] == B[j-1]) ? ma : mi;
                    dp[i][j] = std::max({
                        dp[i-1][j-1] + match_score,
                        dp[i-1][j] + g,
                        dp[i][j-1] + g
                    });
                }
            }
        }));
    }

    for (auto& future : futures) {
        future.wait();
    }
}

std::pair<std::string, std::string> parallel_traceback(const std::vector<std::vector<int>>& dp,
                                                     const std::string& A, const std::string& B,
                                                     int ma, int mi, int g) {
    std::string alignedA;
    std::string alignedB;
    size_t i = A.size();
    size_t j = B.size();
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int match_score = (A[i-1] == B[j-1]) ? ma : mi;
            if (dp[i][j] == dp[i-1][j-1] + match_score) {
                alignedA += A[i-1];
                alignedB += B[j-1];
                --i;
                --j;
                continue;
            }
        }
        
        if (i > 0 && dp[i][j] == dp[i-1][j] + g) {
            alignedA += A[i-1];
            alignedB += '-';
            --i;
        } else {
            alignedA += '-';
            alignedB += B[j-1];
            --j;
        }
    }
    
    std::reverse(alignedA.begin(), alignedA.end());
    std::reverse(alignedB.begin(), alignedB.end());
    return {alignedA, alignedB};
}

AlignmentResult parallel_needleman_wunsch(const std::string& A, const std::string& B, 
                                        int mi, int ma, int g) {
    size_t m = A.size();
    size_t n = B.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    auto antidiagonals = get_antidiagonal_indices(m, n);
    for (const auto& diagonal : antidiagonals) {
        compute_diagonal_parallel(dp, diagonal, A, B, mi, ma, g);
    }
    
    int final_score = dp[m][n];
    auto aligned_sequences = parallel_traceback(dp, A, B, ma, mi, g);
    
    return {final_score, aligned_sequences};
}