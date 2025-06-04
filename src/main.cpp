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
#include "parallel/sw_diagonal_score_only_tp.h"

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
    std::pair<std::string, std::string> sw_aligned = smith_waterman_traceback(sw_dp.first, A, B, 1, -1, -2, sw_dp.second);
    std::cout << "Smith-Waterman Score: " << sw_dp.first[sw_dp.second.first][sw_dp.second.second] << "\n";
    std::cout << "Aligned A: " << sw_aligned.first << "\n";
    std::cout << "Aligned B: " << sw_aligned.second << "\n";

    std::pair<int, std::pair<std::string, std::string>> mm_aligned = myers_miller(A, B, 1, -1, -2);
    std::cout << "Myers-Miller Alignment:\n";
    std::cout << "Aligned A: " << mm_aligned.second.first << "\n";
    std::cout << "Aligned B: " << mm_aligned.second.second << "\n";

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
    int succ_scoreonly_sw = 0;
    
    for (size_t length_seq= 10; length_seq < 13; ++length_seq) {
        for(size_t repetitions = 0; repetitions < 0; ++repetitions) {
            std::string refA = generate_random_dna(1<<length_seq);
            std::string refB = generate_random_dna(1<<length_seq);
            std::cout << "DNA sequences generated."<< "\n";

            auto start_seq_sw = std::chrono::high_resolution_clock::now();
            auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
            auto end_seq_sw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_sw = end_seq_sw - start_seq_sw;
            std::cout << "Sequential SW implementation finished. (score = " << sw_seq.first[sw_seq.second.first][sw_seq.second.second] << ")\n";

            auto simiA = generate_random_dna(1<<length_seq);
            auto simiB = generate_similar_dna(1<<length_seq, 0.5, simiA);
            auto start_seq_nw = std::chrono::high_resolution_clock::now();
            auto nw_dp = needleman_wunsh_dp(simiA, simiB, -1, 1, -2);
            auto nw_alignment = needleman_wunsh_traceback(nw_dp, simiA, simiB, 1, -1, -2);
            int nw_score = 0;
            for (size_t i = 0; i < nw_alignment.first.size(); ++i) {
                char a = nw_alignment.first[i];
                char b = nw_alignment.second[i];
                if (a == '-' || b == '-') nw_score += -2;
                else nw_score += (a == b ? 1 : -1);
            }
            auto end_seq_nw = std::chrono::high_resolution_clock::now();
            std::cout << "Sequential Needle-Wunsch implementation finished. (score = " << nw_score << ").\n";

            auto start_seq_mm = std::chrono::high_resolution_clock::now();
            auto [mm_score, mm_alignment] = myers_miller(simiA, simiB, 1, -1, -2);
            auto end_seq_mm = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_mm = end_seq_mm - start_seq_mm;
            std::cout << "Sequential MM implementation finished. (score = " << mm_score <<").\n";

            auto sw_start_par = std::chrono::high_resolution_clock::now();
            auto sw_par = SmithWatermanWavefrontTp(refA, refB, -1, 1, -2, 8);
            auto sw_end_par = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> sw_duration_par = sw_end_par - sw_start_par;
            std::cout << "Parallel SW implementation finished."<< "\n";

            #ifdef USE_CUDA
            // call CUDA-related functions
            auto start_cuda = std::chrono::high_resolution_clock::now();
            auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
            auto end_cuda = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
            std::cout << "CUDA implementation finished."<< "\n";
            #endif

            auto start_par_scoreonly_sw = std::chrono::high_resolution_clock::now();
            auto sw_par_scoreonly = SmithWatermanWavefront_ScoreOnly(refA, refB, -1, 1, -2, 8);
            auto end_par_scoreonly_sw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par_scoreonly_sw = end_par_scoreonly_sw - start_par_scoreonly_sw;
            std::cout << "Parallel SW (score only) implementation finished. << \n";
            

            succ_par = Check_Matrix_Matrix(sw_seq, sw_par);
            if (succ_par == -1) {
                std::cerr << "Test failed for SW parallel implementation" << length_seq << "\n";
            }
            succ_scoreonly_sw = Check_Matrix_Score(sw_seq, sw_par_scoreonly);
            if (succ_scoreonly_sw == -1) {
                std::cerr << "Test failed for SW sequences of length " << length_seq << "\n";
                break;
            }
            double t_par_sw = sw_duration_par.count();
            double t_seq_sw = duration_seq_sw.count();
            double t_seq_mm = duration_seq_mm.count();
            double t_par_scoreonly_sw = duration_par_scoreonly_sw.count();
            double speedup_par_sw = t_seq_sw / t_par_sw;
            double speedup_par_scoreonly_sw = t_seq_sw / t_par_scoreonly_sw;

            #ifdef USE_CUDA
            succ_cuda = Check_Matrix_Score(sw_seq, cuda_result);
            if (succ_cuda == -1) {
                std::cerr << "Test failed for cuda implementation" << length_seq << "\n";
            }
            double t_cuda = duration_cuda.count();
            double speedup_cuda = t_seq_sw/ t_cuda;
            #endif

            std::cout << "Length: " << (1<<length_seq) << "\n";
            std::cout << "Sequential MM time: " << t_seq_mm << " sec\n";
            std::cout << "Parallel SW time: " << t_par_sw << " sec\n";

            #ifdef USE_CUDA
            std::cout << "Cuda time: " << t_cuda << " sec\n";
            #endif
            std::cout << "Parallel SW time (score only): " << t_par_scoreonly_sw << " sec\n";
            std::cout << "Speedup parallel SW - sequential:  " << speedup_par_sw << "x\n";
            std::cout << "Speedup parallel SW (score only) - sequential:  " << speedup_par_scoreonly_sw << "x\n";

            #ifdef USE_CUDA
            std::cout << "Speedup cuda - sequential:  " << speedup_cuda << "x\n";
            #endif
            std::cout << "\n";
        }
    }

    ///////////////////////// Num threads comparison ////////////////////////////////////////////////////
    // std::vector<int> num_threads = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    // std::string refA = generate_random_dna(1<<15);
    // std::string refB = generate_random_dna(1<<15);
    // auto start_seq_sw = std::chrono::high_resolution_clock::now();
    // auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
    // auto end_seq_sw = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> seq_timing_temp_sw = end_seq_sw - start_seq_sw;
    // double seq_timing_sw = seq_timing_temp_sw.count();
    // std::cout << "Sequential (sw): " << seq_timing_sw << "\n";

    // // For when mm works:
    // // auto start_seq_mm = std::chrono::high_resolution_clock::now();
    // // auto mm_seq = myers_miller_dp(refA, refB, -1, 1, -2);
    // // auto end_seq_mm = std::chrono::high_resolution_clock::now();
    // // std::chrono::duration<double> seq_timing_temp_mm = end_seq_mm - start_seq_sw;
    // // double seq_timing_mm = seq_timing_temp_mm.count();
    // // std::cout << "Sequential (mm): " << seq_timing_mm << "\n";

    // auto timing = function_test_threads(SmithWatermanWavefront, num_threads);
    // std::cout << "Regular: ";
    // for(auto i: timing){
    //     std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_threads(SmithWatermanWavefrontTp, num_threads);
    // std::cout << "Threadpool: ";
    // for(auto i: timing){
    //     std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_threads(SmithWatermanWavefront_ScoreOnly, num_threads);
    // std::cout << "Score only: ";
    // for(auto i: timing){
    //     std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_threads(SmithWatermanWavefront_ScoreOnly_Tp, num_threads);
    // std::cout << "Score only + TP: ";
    // for(auto i: timing){
    //     std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    // }
    // std::cout << "\n\n";

    //////////////////////////////SIZE COMPARISON///////////////////////////////////////////////////////////////////
    // std::vector<int> sizes = {1<<11, 1<<12, 1<<13, 1<<14, 1<<15, 1<<16, 1<<17};

    // std::vector<double> seq_timings;
    // std::vector<double> gpu_timings;
    // for(int i: sizes){
    //     std::string refA = generate_random_dna(i);
    //     std::string refB = generate_random_dna(i);  
        
    //     auto start_seq_sw = std::chrono::high_resolution_clock::now();
    //     auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
    //     auto end_seq_sw = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> seq_timing_temp = end_seq_sw - start_seq_sw;
    //     seq_timings.push_back(seq_timing_temp.count());
    //     int res_seq = sw_seq.first[sw_seq.second.first][sw_seq.second.second];

        // #ifdef USE_CUDA
        // auto start_cuda = std::chrono::high_resolution_clock::now();
        // auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
        // auto end_cuda = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
        // gpu_timings.push_back(duration_cuda.count());
        // int res_gpu = std::get<0>(cuda_result);
        // if(res_gpu!=res_seq){
        //     std::cout << "Incorrect results: " << res_seq << " / " << res_gpu << "\n";
        // }
        // #endif
    // }
    // std::cout << "Seq (SW): ";
    // for(auto i: seq_timings){
    //     std::cout << i << " ";
    // }
    // std::cout << "\n\n";

    // #ifdef USE_CUDA
    // std::cout << "GPU: ";
    // for(int i=0; i<gpu_timings.size(); ++i){
    //     std::cout << " " << gpu_timings[i] << " (" << seq_timings[i]/gpu_timings[i] << ") |";
    // }
    // std::cout << "\n\n";
    // #endif

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

    // auto timing = function_test_size(SmithWatermanWavefront_ScoreOnly_Tp, sizes, sizes, 8);
    // std::cout << "Parallel score only + Tp: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    


    //////////////////////// Similarity comparison //////////////////////////////////////////////////////////////
    // std::vector<double> similarities = {0,0.2,0.4,0.6,0.8,1.0};
    // std::vector<double> seq_timings;
    // int size = 1<<14;
    // for(double i: similarities){
    //     std::string refA = generate_random_dna(size);
    //     std::string refB = generate_similar_dna(size, i, refA);  
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

    // auto timing = function_test_similarity(SmithWatermanWavefront, similarities, 8);
    // std::cout << "Standard parallel: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_similarity(SmithWatermanWavefrontTp, similarities, 8);
    // std::cout << "Parallel TP: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_similarity(SmithWatermanWavefront_ScoreOnly, similarities, 8);
    // std::cout << "Parallel score only: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // timing = function_test_similarity(SmithWatermanWavefront_ScoreOnly_Tp, similarities, 8);
    // std::cout << "Parallel score only + Tp: ";
    // for(int i=0; i<timing.size(); ++i){
    //     std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    // }
    // std::cout << "\n\n";

    // std::vector<double> gpu_timings;
    // #ifdef USE_CUDA
    // for(double i: similarities){
    //     std::string refA = generate_random_dna(size);
    //     std::string refB = generate_similar_dna(size, i, refA);  
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