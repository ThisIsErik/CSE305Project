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

// Timer for functions with 6 args (A, B, mi, ma, g, num_threads)
template <typename Func>
auto function_test_timer(Func func,
                        const std::string& A,
                        const std::string& B,
                        int mi,
                        int ma,
                        int g,
                        size_t num_threads)
    -> std::enable_if_t<std::is_invocable_v<Func, const std::string&, const std::string&, int, int, int, size_t>, double>
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func(A, B, mi, ma, g, num_threads);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

// Timer for functions with 5 args (A, B, mi, ma, g)
template <typename Func>
auto function_test_timer(Func func,
                        const std::string& A,
                        const std::string& B,
                        int mi,
                        int ma,
                        int g,
                        size_t /*placehold for num_threads*/)
    -> std::enable_if_t<!std::is_invocable_v<Func, const std::string&, const std::string&, int, int, int, size_t> &&
                        std::is_invocable_v<Func, const std::string&, const std::string&, int, int, int>, double>
{
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func(A, B, mi, ma, g);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
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
        double time = function_test_timer(func, A, B, -1, 1, -2, i);
        times.push_back(time);
    }
    return times;
}

template <typename Func>
std::vector<double> function_test_size(Func func,
                                     std::vector<int> size_a = {1<<10, 1<<11, 1<<12, 1<<13, 1<<14, 1<<15}, 
                                     std::vector<int> size_b = {1<<10, 1<<11, 1<<12, 1<<13, 1<<14, 1<<15},
                                     size_t num_threads = 8) {
    std::vector<double> times;
    for(int i=0; i<size_a.size(); ++i){
        std::string A = generate_random_dna(size_a[i]);
        std::string B = generate_random_dna(size_b[i]);
        double time = function_test_timer(func, A, B, -1, 1, -2, num_threads);
        times.push_back(time);
    }
    return times;
}

template <typename Func>
std::vector<double> function_test_similarity(Func func,
                                    std::vector<double> similarity = {0,0.2,0.4,0.6,0.8,1.0},
                                    size_t num_threads = 8) {
    std::vector<double> times;
    size_t length = 1<<14;
    for(double i: similarity){
        std::string A = generate_random_dna(length);
        std::string B = generate_similar_dna(length, i, A);
        double time = function_test_timer(func, A, B, -1, 1, -2, num_threads);
        times.push_back(time);
    }
    return times;
}


#endif // PERF_TEST