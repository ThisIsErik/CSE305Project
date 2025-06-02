#include <iostream>
#include "utils/random_dna.h"
//#include "utils/timing_utils.h"
#include "sequential/nw.h"
#include "sequential/sw.h"
#include "sequential/mm.h"
#include "parallel/sw_parallel_database.h"
#include "parallel/sw_diagonal_wavefront_tp.h"

#include "parallel/sw_diagonal_wavefront.h"
#include "parallel/sw_diagonal_score_only.h"
<<<<<<< Updated upstream
<<<<<<< Updated upstream
=======
=======
>>>>>>> Stashed changes

#include "parallel/mm/mm_parallel_database.h"
#include "parallel/mm/mm_diagonal_wavefront_tp.h"
#include "parallel/mm/mm_diagonal_wavefront.h"
#include "parallel/mm/mm_diagonal_score_only.h"

>>>>>>> Stashed changes
#include "tests.h"
#include <vector>
#include <iomanip>
#include <chrono>

#ifdef USE_CUDA
#include "gpu/CUDAlign.cuh"
#endif


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

<<<<<<< Updated upstream
<<<<<<< Updated upstream
    // Myers-Miller
    std::pair<std::string, std::string> mm_aligned = myers_miller(A, B, 1, -1, -2);
    int mm_score = myers_miller_score(A, B, 1, -1, -2);
    std::cout << "Myers-Miller Score: " << mm_score << "\n";
    std::cout << "Aligned A: " << mm_aligned.first << "\n";
    std::cout << "Aligned B: " << mm_aligned.second << "\n";
=======
=======
>>>>>>> Stashed changes
    // Modified MM section to match the implementation
    std::tuple<std::string, std::string, int> mm_aligned = myers_miller_align(A, B, -1, 1, -2, -1);
    std::cout << "Myers-Miller Score: " << std::get<2>(mm_aligned)<< "\n";
    std::cout << "Aligned A: " << std::get<0>(mm_aligned) << "\n";
    std::cout << "Aligned B: " << std::get<1>(mm_aligned) << "\n";
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

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
    int succ_par = 0;
    int succ_cuda = 0;
<<<<<<< Updated upstream
    int succ_scoreonly = 0;
    for (size_t length_seq= 16; length_seq < 17; ++length_seq) {
        for(size_t repetitions = 0; repetitions < 1; ++repetitions) {
=======
    int succ_scoreonly_sw = 0;
    for (size_t length_seq= 10; length_seq < 15; ++length_seq) {
        for(size_t repetitions = 0; repetitions < 3; ++repetitions) {
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
            std::string refA = generate_random_dna(1<<length_seq);
            std::string refB = generate_random_dna(1<<length_seq);
            std::cout << "DNA sequences generated."<< "\n";

            auto start_seq = std::chrono::high_resolution_clock::now();
            auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
            auto end_seq = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq = end_seq - start_seq;
            std::cout << "Sequential implementation finished. (score = " << sw_seq.first[sw_seq.second.first][sw_seq.second.second] << ")\n";

<<<<<<< Updated upstream
            auto start_par = std::chrono::high_resolution_clock::now();
=======
            auto start_seq_mm = std::chrono::high_resolution_clock::now();
            auto mm_seq = myers_miller_align(refA, refB, -1, 1, -2, -2);
            auto end_seq_mm = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_mm = end_seq_mm - start_seq_mm;
            std::cout << "Sequential MM implementation finished. (score = " << std::get<2>(mm_seq) << ")\n";

            auto sw_start_par = std::chrono::high_resolution_clock::now();
>>>>>>> Stashed changes
            auto sw_par = SmithWatermanWavefrontTp(refA, refB, -1, 1, -2, 8);
            auto end_par = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par = end_par - start_par;
            std::cout << "Parallel implementation finished."<< "\n";

            #ifdef USE_CUDA
            // call CUDA-related functions
            auto start_cuda = std::chrono::high_resolution_clock::now();
            auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
            auto end_cuda = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
            std::cout << "CUDA implementation finished."<< "\n";
            #endif

            auto start_par_scoreonly = std::chrono::high_resolution_clock::now();
            auto sw_par_scoreonly = SmithWatermanWavefront_ScoreOnly(refA, refB, -1, 1, -2, 8);
            auto end_par_scoreonly = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par_scoreonly = end_par_scoreonly - start_par_scoreonly;
            std::cout << "Parallel (score only) implementation finished. << \n";
            
<<<<<<< Updated upstream
<<<<<<< Updated upstream
=======
=======
>>>>>>> Stashed changes
            // For MM verification, we need to compare with NW since both are global alignment
            auto nw_seq = needleman_wunsh_dp(refA, refB, -1, 1, -2);
            auto nw_aligned = needleman_wunsh_traceback(nw_seq, refA, refB, 1, -1, -2);
            
            if (nw_aligned.first == std::get<0>(mm_aligned) && nw_aligned.second == std::get<1>(mm_aligned)) {
                std::cout << "MM alignment matches NW alignment (correct)\n";
            } else {
                std::cerr << "MM alignment differs from NW alignment!\n";
            }

            // Rest of the original checks remain the same
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
            succ_par = Check_Matrix_Matrix(sw_seq, sw_par);
            if (succ_par == -1) {
                std::cerr << "Test failed for parallel implementation" << length_seq << "\n";
            }
            succ_scoreonly = Check_Matrix_Score(sw_seq, sw_par_scoreonly);
            if (succ_scoreonly == -1) {
                std::cerr << "Test failed for sequences of length " << length_seq << "\n";
                break;
            }

<<<<<<< Updated upstream
<<<<<<< Updated upstream
            double t_par = duration_par.count();
            double t_seq = duration_seq.count();
            double t_par_scoreonly = duration_par_scoreonly.count();
            double speedup_par = t_seq / t_par;
            double speedup_par_scoreonly = t_seq / t_par_scoreonly;
=======
=======
>>>>>>> Stashed changes
            double t_par_sw = sw_duration_par.count();
            double t_par_mm = mm_duration_par.count();
            double t_seq_sw = duration_seq_sw.count();
            double t_seq_mm = duration_seq_mm.count();
            double t_par_scoreonly_sw = duration_par_scoreonly_sw.count();
            double t_par_scoreonly_mm = duration_par_scoreonly_mm.count();
            double speedup_par_sw = t_seq_sw / t_par_sw;
            double speedup_par_mm = t_seq_mm / t_par_mm;
            double speedup_par_scoreonly_sw = t_seq_sw / t_par_scoreonly_sw;
            double speedup_par_scoreonly_mm = t_seq_mm / t_par_scoreonly_mm;
>>>>>>> Stashed changes

            #ifdef USE_CUDA
            succ_cuda = Check_Matrix_Score(sw_seq, cuda_result);
            if (succ_cuda == -1) {
                std::cerr << "Test failed for cuda implementation" << length_seq << "\n";
            }
            double t_cuda = duration_cuda.count();
            double speedup_cuda = t_seq/ t_cuda;
            #endif

            std::cout << "Length: " << (1<<length_seq) << "\n";
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            std::cout << "Sequential time: " << t_seq << " sec\n";
            std::cout << "Parallel time: " << t_par << " sec\n";
=======
=======
>>>>>>> Stashed changes
            std::cout << "Sequential SW time: " << t_seq_sw << " sec\n";
            std::cout << "Parallel SW time: " << t_par_sw << " sec\n";
            std::cout << "Sequential MM time: " << t_seq_mm << " sec\n";
            std::cout << "Parallel MM time: " << t_par_mm << " sec\n";

>>>>>>> Stashed changes
            #ifdef USE_CUDA
            std::cout << "Cuda time: " << t_cuda << " sec\n";
            #endif
            std::cout << "Parallel time (score only): " << t_par_scoreonly << " sec\n";
            std::cout << "Speedup parallel - sequential:  " << speedup_par << "x\n";
            std::cout << "Speedup parallel (score only) - sequential:  " << speedup_par_scoreonly << "x\n";
            #ifdef USE_CUDA
            std::cout << "Speedup cuda - sequential:  " << speedup_cuda << "x\n";
            #endif
    }
<<<<<<< Updated upstream
}
=======

    // Rest of the original file remains exactly the same
    /////////// Time performance tests SW ///////////////
    auto timing1 = function_test_threads(SmithWatermanWavefrontTp, {1,2,4,8,16});
    for(auto i: timing1){
        std::cout << i << " ";
    }
    std::cout << "\n";

    auto timingscoreonly1 = function_test_size(SmithWatermanWavefront_ScoreOnly, {1<<10,1<<11,1<<12,1<<13,1<<14},{1<<10,1<<11,1<<12,1<<13,1<<14}, 8);
    auto timingseq1 = function_test_size(smith_waterman_dp, {1<<10,1<<11,1<<12,1<<13,1<<14},{1<<10,1<<11,1<<12,1<<13,1<<14});
    for(size_t i=0; i<timingseq1.size(); ++i){
        std::cout << timingseq1[i]/timingscoreonly1[i] << " ";
    }
    std::cout << "\n";

    /////////// Time performance tests MM ///////////////
    auto timing2 = function_test_threads(MyersMillerWavefrontTp, {1,2,4,8,16});
    for(auto i: timing2){
        std::cout << i << " ";
    }
    std::cout << "\n";

    auto timingscoreonly2 = function_test_size(MyersMillerWavefront_ScoreOnly, {1<<10,1<<11,1<<12,1<<13,1<<14},{1<<10,1<<11,1<<12,1<<13,1<<14}, 8);
    auto timingseq2 = function_test_size(needleman_wunsh_dp, {1<<10,1<<11,1<<12,1<<13,1<<14},{1<<10,1<<11,1<<12,1<<13,1<<14});
    for(size_t i=0; i<timingseq2.size(); ++i){
        std::cout << timingseq2[i]/timingscoreonly2[i] << " ";
    }
    std::cout << "\n";
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

    return 0;
}