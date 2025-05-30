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
#include "tests.h"
#include "utils/performance_test.h"
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
    int succ_scoreonly = 0;
    for (size_t length_seq= 16; length_seq < 18; ++length_seq) {
        for(size_t repetitions = 0; repetitions < 1; ++repetitions) {
            std::string refA = generate_random_dna(1<<length_seq);
            std::string refB = generate_random_dna(1<<length_seq);
            std::cout << "DNA sequences generated."<< "\n";

            auto start_seq_sw = std::chrono::high_resolution_clock::now();
            auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
            auto end_seq_sw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_sw = end_seq_sw - start_seq_sw;
            std::cout << "Sequential SW implementation finished. (score = " << sw_seq.first[sw_seq.second.first][sw_seq.second.second] << ")\n";

            auto start_seq_mm = std::chrono::high_resolution_clock::now();
            auto mm_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
            auto end_seq_mm = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_mm = end_seq_mm - start_seq_mm;
            std::cout << "Sequential MM implementation finished. (score = " << mm_seq.first[mm_seq.second.first][mm_seq.second.second] << ")\n";

            auto start_par = std::chrono::high_resolution_clock::now();
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
            

            int succ_seq = Check_Matrix_Matrix(sw_seq, mm_seq);
            if (succ_seq == -1) {
                std::cerr << "Test failed for MM implementation" << length_seq << "\n";
            }
            succ_par = Check_Matrix_Matrix(sw_seq, sw_par);
            if (succ_par == -1) {
                std::cerr << "Test failed for parallel implementation" << length_seq << "\n";
            }
            succ_scoreonly = Check_Matrix_Score(sw_seq, sw_par_scoreonly);
            if (succ_scoreonly == -1) {
                std::cerr << "Test failed for sequences of length " << length_seq << "\n";
                break;
            }

            double t_par = duration_par.count();
            double t_seq_sw = duration_seq_sw.count();
            double t_seq_mm = duration_seq_mm.count();
            double t_par_scoreonly = duration_par_scoreonly.count();
            double speedup_par = t_seq_sw / t_par;
            double speedup_par_scoreonly = t_seq_sw / t_par_scoreonly;

            #ifdef USE_CUDA
            succ_cuda = Check_Matrix_Score(sw_seq, cuda_result);
            if (succ_cuda == -1) {
                std::cerr << "Test failed for cuda implementation" << length_seq << "\n";
            }
            double t_cuda = duration_cuda.count();
            double speedup_cuda = t_seq_sw/ t_cuda;
            #endif

            std::cout << "Length: " << (1<<length_seq) << "\n";
            std::cout << "Sequential time: " << t_seq_sw << " sec\n";
            std::cout << "Parallel time: " << t_par << " sec\n";
            #ifdef USE_CUDA
            std::cout << "Cuda time: " << t_cuda << " sec\n";
            #endif
            std::cout << "Parallel time (score only): " << t_par_scoreonly << " sec\n";
            std::cout << "Speedup parallel - sequential:  " << speedup_par << "x\n";
            std::cout << "Speedup parallel (score only) - sequential:  " << speedup_par_scoreonly << "x\n";
            #ifdef USE_CUDA
            std::cout << "Speedup cuda - sequential:  " << speedup_cuda << "x\n";
            #endif
            std::cout << "\n";
        }
    }

    /////////// Time performance tests ///////////////
    // std::string refA = generate_random_dna(1<<15);
    // std::string refB = generate_random_dna(1<<15);
    // auto start_seq_sw = std::chrono::high_resolution_clock::now();
    // auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
    // auto end_seq_sw = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> seq_timing_temp = end_seq_sw - start_seq_sw;
    // double seq_timing = seq_timing_temp.count();
    // std::cout << "Sequential: " << seq_timing << "\n";

    // auto timing = function_test_threads(SmithWatermanWavefront, {1,2,4,8,16});
    // std::cout << "Regular: ";
    // for(auto i: timing){
    //     std::cout << i <<"(" << seq_timing/i<< ") ";
    // }
    // std::cout << "\n\n";

    // timing = function_test_threads(SmithWatermanWavefrontTp, {1,2,4,8,16});
    // std::cout << "Threadpool: ";
    // for(auto i: timing){
    //     std::cout << i <<"(" << seq_timing/i<< ") ";
    // }
    // std::cout << "\n\n";

    // timing = function_test_threads(SmithWatermanWavefront_ScoreOnly, {1,2,4,8,16});
    // std::cout << "Score only: ";
    // for(auto i: timing){
    //     std::cout << i <<"(" << seq_timing/i<< ") ";
    // }
    // std::cout << "\n\n";

    

    // auto timingscoreonly = function_test_size(SmithWatermanWavefront_ScoreOnly, {1<<10,1<<11,1<<12,1<<13,1<<14},{1<<10,1<<11,1<<12,1<<13,1<<14}, 8);
    // auto timingseq = function_test_size(smith_waterman_dp, {1<<10,1<<11,1<<12,1<<13,1<<14},{1<<10,1<<11,1<<12,1<<13,1<<14});
    // for(int i=0; i<timingseq.size(); ++i){
    //     std::cout << timingseq[i]/timingscoreonly[i] << " ";
    // }
    // std::cout << "\n";

    
    // std::vector<int> sizes = {1<<11, 1<<12, 1<<13, 1<<14, 1<<15, 1<<16, 1<<17};
    // std::vector<double> seq_timings;
    // for(int i: sizes){
    //     std::string refA = generate_random_dna(i);
    //     std::string refB = generate_random_dna(i);  
    //     auto start_seq_sw = std::chrono::high_resolution_clock::now();
    //     auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
    //     auto end_seq_sw = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> seq_timing_temp = end_seq_sw - start_seq_sw;
    //     double seq_timing = seq_timing_temp.count();
    //     seq_timings.push_back(seq_timing);
    // }
    // std::cout << "Seq (SW): ";
    // for(auto i: seq_timings){
    //     std::cout << i << " ";
    // }
    // std::cout << "\n\n";

    // auto timing = function_test_size(SmithWatermanWavefront, sizes, sizes, 8);
    // std::cout << "Standard parallel: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_size(SmithWatermanWavefrontTp, sizes, sizes, 8);
    // std::cout << "Parallel TP: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_size(SmithWatermanWavefront_ScoreOnly, sizes, sizes, 8);
    // std::cout << "Parallel score only: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // std::vector<double> gpu_timings;
    // #ifdef USE_CUDA
    // for(int i: sizes){
    //     std::string refA = generate_random_dna(i);
    //     std::string refB = generate_random_dna(i);  
    //     auto start_cuda = std::chrono::high_resolution_clock::now();
    //     auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
    //     auto end_cuda = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
    //     gpu_timings.push_back(duration_cuda.count());
    // }
    // #endif
    // std::cout << "GPU: ";
    // for(int i=0; i<gpu_timings.size(); ++i){
    //     std::cout << " " << gpu_timings[i] << " (" << seq_timings[i]/gpu_timings[i] << ") |";
    // }
    // std::cout << "\n\n";
    
    

    return 0;
}