#ifndef PERF_TEST
#define PERF_TEST

#include <string>
#include <tuple>
#include <vector>
#include "utils/types.h"
#include "utils/random_dna.h"
#include <chrono>
#include <functional>

//Functions and their inputs:
//GPU:
//(score, i, j) CUDAlign(seq A, seq B, (mi,ma,g))
//CPU:
//* (score, i, j) SmithWatermanWavefront_ScoreOnly(seq A, seq B, (mi,ma,g), num_threads)
//* (DP, (i,j)) SmithWatermanWavefrontTp(seq A, seq B, (mi,ma,g), num_threads
//* (DP, (i,j)) SmithWatermanWavefront(seq A, seq B, (mi,ma,g), num_threads
//Sequential:
//* (DP, (i,j)) smith_waterman_dp(seq A, seq B, (mi,ma,g));

template <typename Func, typename... Args>
auto function_test_timer(Func func, Args&&... args)
{
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    auto result = std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return std::make_tuple(result, elapsed.count());
}

template <typename Func, typename... Args>
auto function_test_random(Func func, size_t lenA, size_t lenB, int mi, int ma, int g, double similarity, Args&&... args) {
    std::string A = generate_random_dna(lenA);
    std::string B = generate_similar_dna(lenB, similarity, A);
    return function_test_timer(func, A, B, mi, ma, g, std::forward<Args>(args)...);
}

template <typename Func>
std::vector<double> function_test_threads(Func func, std::vector<int> num_threads = {1,2,4,8,16}) {
    std::string A = generate_random_dna(1<<14);
    std::string B = generate_random_dna(1<<14);
    std::vector<double> times;
    for(int i: num_threads){
        auto [result, time] = function_test_timer(func, A, B, -1, 1, -2, i);
        times.push_back(time);
    }
    return times;
}


#endif // PERF_TEST