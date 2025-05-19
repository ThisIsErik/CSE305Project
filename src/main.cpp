#include <iostream>
#include "utils/random_dna.h"
//#include "utils/timing_utils.h"
#include "sequential/nw.h"
#include "sequential/sw.h"
#include "parallel/sw_parallel_database.h"
// #include "parallel/sw_diagonal_wavefront_tp.h"

#include "parallel/sw_diagonal_wavefront.h"
#include "parallel/sw_diagonal_score_only.h"
#include "tests.h"
#include <vector>
#include <iomanip>
#include <chrono>


int main() {
    std::string A = "TAGC";
    std::string B = "TAGTC";

    std::cout << "Sequence A: " << A.substr(0,10) <<  ((A.size()>10)?"...\n":"\n");
    std::cout << "Sequence B: " << B.substr(0,10) <<  ((B.size()>10)?"...\n":"\n");

    std::vector<std::vector<int>> nw_dp = needleman_wunsh_dp(A, B, -1, 1, -2);
    std::pair<std::string, std::string> nw_aligned = needleman_wunsh_traceback(nw_dp, A, B, 1, -1, -2);
    std::cout << "Needleman-Wunsh Score: " << nw_dp[A.size()][B.size()] << "\n";
    std::cout << "Aligned A: " << nw_aligned.first << "\n";
    std::cout << "Aligned B: " << nw_aligned.second << "\n";

    std::pair<std::vector<std::vector<int>>, std::pair<int,int>> sw_dp = smith_waterman_dp(A, B, -1, 1, -2);
    std::pair<std::string, std::string> aligned = smith_waterman_traceback(sw_dp.first, A, B, 1, -1, -2, sw_dp.second);
    std::cout << "Smith-Waterman Score: " << sw_dp.first[sw_dp.second.first][sw_dp.second.second] << "\n";
    std::cout << "Aligned A: " << aligned.first << "\n";
    std::cout << "Aligned B: " << aligned.second << "\n";

    //First parallel implementation. Create a dictionary of 1000 sequences of bases
    std::string query = generate_random_dna(1000);
    std::vector<std::string> database;
    database.reserve(1000);
    double length_step = static_cast<double>(1000 - 500) / (1000 - 1); //create sequences of length 1000,999,...,500
    for (int i = 0; i < 1000; ++i) {
        int current_length = 1000 - static_cast<int>(i * length_step);
        database.push_back(generate_random_dna(current_length));
    }
    std::cout << "Database generated "<< "\n";

    //#################  Sanity check that diagonal wavefront works  #######################
    int succ = 0;
    int succ_scoreonly = 0;
    for (size_t i = 13; i < 16; ++i) {
        for(size_t j=0; j < 1; ++j){
            std::string refA = generate_random_dna(1<<i);
            std::string refB = generate_random_dna(1<<i);
            std::cout << "DNA sequences generated."<< "\n";

            auto start_seq = std::chrono::high_resolution_clock::now();
            auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
            auto end_seq = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq = end_seq - start_seq;
            std::cout << "Sequential implementation finished. (score = " << sw_seq.first[sw_seq.second.first][sw_seq.second.second] << ")\n";

            // printAntidiagonals(sw_seq.first);

            // auto start_par = std::chrono::high_resolution_clock::now();
            // auto sw_par = SmithWatermanWavefrontTp(refA, refB, -1, 1, -2, 8);
            // auto end_par = std::chrono::high_resolution_clock::now();
            // std::chrono::duration<double> duration_par = end_par - start_par;
            // std::cout << "Parallel implementation finished."<< "\n";

            auto start_par_scoreonly = std::chrono::high_resolution_clock::now();
            auto [score, pos_i, pos_j] = SmithWatermanWavefront_ScoreOnly(refA, refB, -1, 1, -2, 1);
            auto end_par_scoreonly = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par_scoreonly = end_par_scoreonly - start_par_scoreonly;
            std::cout << "Parallel (score only) implementation finished. (score = " << score << ")\n";

            // succ = CheckSWWavefront(sw_seq, sw_par);
            succ_scoreonly = CheckSWWavefront_ScoreOnly(sw_seq, std::make_tuple(score, pos_i, pos_j));
            if (succ == -1) {
                std::cerr << "Test failed for sequences of length " << i << "\n";
                break;
            }

            double t_seq = duration_seq.count();
            double t_par = duration_par_scoreonly.count();
            double speedup = t_seq / t_par;

            std::cout << "Length: " << (1<<i) << "\n";
            std::cout << "Sequential time: " << t_seq << " sec\n";
            std::cout << "Parallel time: " << t_par << " sec\n";
            // std::cout << "Parallel time (score only): " << t_par_scoreonly << " sec\n";
            std::cout << "Speedup:  " << speedup << "x\n\n";
        }

    }

    return 0;
}