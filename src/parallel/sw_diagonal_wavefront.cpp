#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>

typedef std::pair<std::vector<std::vector<int>>, std::pair<int, int>> SWResult;

void AntiDiagonalAux(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& dp,
    int& max_val,
    std::pair<int, int>& max_pos,
    std::mutex& max_mutex
) {
    int m = A.size();
    int n = B.size();
    for (int i = start_i; i <= end_i; ++i) {
        int j = diagonal - i;
        if (i <= 0 || i > m || j <= 0 || j > n)
            continue;
        int p = (A[i-1] == B[j-1]) ? ma : mi;
        int val = std::max({
            dp[i-1][j-1] + p,
            dp[i][j-1] + g,
            dp[i-1][j] + g,
            0
        });
        dp[i][j] = val;

        std::lock_guard<std::mutex> lock(max_mutex);
        if (val > max_val) {
                max_val = val;
                max_pos = {i, j};
        }
    }
}

SWResult SmithWatermanWavefront(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
) {
    const int m = A.size();
    const int n = B.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    int max_val = 0;
    std::pair<int, int> max_pos = {0, 0};
    std::mutex max_mutex;

    //For each antidiagonal
    for (int d = 1; d <= m + n; ++d) {
        //Bounds for the curr antidiagonal
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        int num_cells = i_end - i_start + 1;

        // if (num_cells < num_threads * 2) {
        //     AntiDiagonalAux(A, B, mi, ma, g, d, i_start, i_end, dp, max_val, max_pos, max_mutex);
        //     continue;
        // }

                // Split work across threads
        int block_size = num_cells / num_threads;
        std::vector<std::thread> workers(num_threads - 1);

        int current_i = i_start;
        for (size_t i = 0; i < num_threads - 1; ++i) {
            int block_end = current_i + block_size;
            workers[i] = std::thread(
                AntiDiagonalAux,
                std::ref(A), std::ref(B),
                mi, ma, g,
                d,
                current_i, block_end,
                std::ref(dp),
                std::ref(max_val), std::ref(max_pos),
                std::ref(max_mutex)
            );
            current_i = block_end + 1;
        }

        AntiDiagonalAux(A, B, mi, ma, g, d, current_i, i_end, dp, max_val, max_pos, max_mutex);

        //Sync all workers of this antidiagonal as next one depends on these values
        for (auto& t : workers) {
            t.join();
        }
    }

    return {dp, max_pos};
}