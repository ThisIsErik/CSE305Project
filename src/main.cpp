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
#include "parallel/nw_parallel.h"

#include "tests.h"
#include "utils/performance_test.h"
#include <vector>
#include <iomanip>
#include <chrono>

#ifdef USE_CUDA
#include "gpu/CUDAlign.cuh"
#endif


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./run <mode>\n";
        std::cout << "Available modes: corectness, thread, size, similarity, memory\n";
        return 1;
    }
    std::string mode = argv[1];

    if (mode == "corectness") {
    //#################  Sanity check that sequential baselines works  #######################
    std::string A = "TAGC";
    std::string B = "TAGTC";

    std::cout << "Sequence A: " << A.substr(0,10) <<  ((A.size()>10)?"...\n":"\n");
    std::cout << "Sequence B: " << B.substr(0,10) <<  ((B.size()>10)?"...\n":"\n");

    std::vector<std::vector<int>> nw_dp = needleman_wunsh_dp(A, B, -1, 1, -2);
    std::pair<std::string, std::string> nw_aligned = needleman_wunsh_traceback(nw_dp, A, B, 1, -1, -2);
    std::cout << "Needleman-Wunsh Score: " << nw_dp[A.size()][B.size()] << "\n";
    std::cout << "Aligned A: " << nw_aligned.first << "\n";
    std::cout << "Aligned B: " << nw_aligned.second << "\n";

    std::pair<int, std::pair<std::string, std::string>> mm_aligned = myers_miller(A, B, 1, -1, -2);
    std::cout << "Myers-Miller Alignment:\n";
    std::cout << "Aligned A: " << mm_aligned.second.first << "\n";
    std::cout << "Aligned B: " << mm_aligned.second.second << "\n";

    std::pair<std::vector<std::vector<int>>, std::pair<int,int>> sw_dp = smith_waterman_dp(A, B, -1, 1, -2);
    std::pair<std::string, std::string> sw_aligned = smith_waterman_traceback(sw_dp.first, A, B, 1, -1, -2, sw_dp.second);
    std::cout << "Smith-Waterman Score: " << sw_dp.first[sw_dp.second.first][sw_dp.second.second] << "\n";
    std::cout << "Aligned A: " << sw_aligned.first << "\n";
    std::cout << "Aligned B: " << sw_aligned.second << "\n";

    //#################  Check that paralllel algorithms work  #######################
    int succ_par = 0;
    int succ_cuda = 0;
    int succ_scoreonly_sw = 0;
    
    for (size_t length_seq= 10; length_seq < 14; ++length_seq) {
        for(size_t repetitions = 0; repetitions < 1; ++repetitions) {
            std::string refA = generate_random_dna(1<<length_seq);
            std::string refB = generate_random_dna(1<<length_seq);
            std::cout << "DNA sequences generated."<< "\n";

            auto start_seq_sw = std::chrono::high_resolution_clock::now();
            auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
            auto end_seq_sw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_sw = end_seq_sw - start_seq_sw;
            std::cout << "Sequential SW implementation finished. (score = " << sw_seq.first[sw_seq.second.first][sw_seq.second.second] << ")\n";

            auto sw_start_par = std::chrono::high_resolution_clock::now();
            auto sw_par = SmithWatermanWavefrontTp(refA, refB, -1, 1, -2, 8);
            auto sw_end_par = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> sw_duration_par = sw_end_par - sw_start_par;
            std::cout << "Parallel SW (TP) implementation finished. (score = " << sw_par.first[sw_par.second.first][sw_par.second.second] << ")\n";

            auto start_par_scoreonly_sw = std::chrono::high_resolution_clock::now();
            auto sw_par_scoreonly = SmithWatermanWavefront_ScoreOnly(refA, refB, -1, 1, -2, 8);
            auto end_par_scoreonly_sw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par_scoreonly_sw = end_par_scoreonly_sw - start_par_scoreonly_sw;
            std::cout << "Parallel SW (score only) implementation finished. (score = " << std::get<0>(sw_par_scoreonly) << ")\n";

            auto start_par_scoreonly_tp_sw = std::chrono::high_resolution_clock::now();
            auto sw_par_scoreonly_tp = SmithWatermanWavefront_ScoreOnly_Tp(refA, refB, -1, 1, -2, 8);
            auto end_par_scoreonly_tp_sw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par_scoreonly_tp_sw = end_par_scoreonly_tp_sw - start_par_scoreonly_tp_sw;
            std::cout << "Parallel SW (score only + TP) implementation finished. (score = " << std::get<0>(sw_par_scoreonly_tp) << ")\n";

            auto simiA = generate_random_dna(1<<length_seq);
            auto simiB = generate_similar_dna(1<<(length_seq-1), 0.99, simiA);
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
            std::chrono::duration<double> duration_seq_nw = end_seq_nw - start_seq_nw;
            std::cout << "Sequential Needle-Wunsch implementation finished. (score = " << nw_score << ").\n";

            auto start_par_nw = std::chrono::high_resolution_clock::now();
            auto [par_nw_score, par_nw_alignment] = NeedlemanWunschWavefront(simiA, simiB, -1, 1, -2, 8);
            auto end_par_nw = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_par_nw = end_par_nw - start_par_nw;
            std::cout << "Parallel Needle-Wunsch implementation finished. (score = " << par_nw_score <<").\n";

            auto start_seq_mm = std::chrono::high_resolution_clock::now();
            auto [mm_score, mm_alignment] = myers_miller(simiA, simiB, 1, -1, -2);
            auto end_seq_mm = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_seq_mm = end_seq_mm - start_seq_mm;
            std::cout << "Sequential MM implementation finished. (score = " << mm_score <<").\n";

            #ifdef USE_CUDA
            // call CUDA-related functions
            auto start_cuda = std::chrono::high_resolution_clock::now();
            auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
            auto end_cuda = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
            std::cout << "CUDA implementation finished."<< "\n";
            #endif
            
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
            double t_par_nw = duration_par_nw.count();
            double t_seq_nw = duration_seq_nw.count();
            double t_par_scoreonly_sw = duration_par_scoreonly_sw.count();
            double speedup_par_sw = t_seq_sw / t_par_sw;
            double speedup_par_scoreonly_sw = t_seq_sw / t_par_scoreonly_sw;
            double speedup_par_nw = t_seq_nw / t_par_nw;

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
            std::cout << "Speedup parallel Needle-Wunsch: " << speedup_par_nw << "x\n";
            std::cout << "Speedup parallel SW - sequential:  " << speedup_par_sw << "x\n";
            std::cout << "Speedup parallel SW (score only) - sequential:  " << speedup_par_scoreonly_sw << "x\n";

            #ifdef USE_CUDA
            std::cout << "Speedup cuda - sequential:  " << speedup_cuda << "x\n";
            #endif
            std::cout << "\n";
            }
        }
    }
    if (mode == "thread") {
    /////////////////////// Number threads comparison ///////////////////////////////////////
    std::vector<int> num_threads = {2,3,4};
    std::string refA = generate_random_dna(1<<15);
    std::string refB = generate_random_dna(1<<15);
    auto start_seq_sw = std::chrono::high_resolution_clock::now();
    auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
    auto end_seq_sw = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> seq_timing_temp_sw = end_seq_sw - start_seq_sw;
    double seq_timing_sw = seq_timing_temp_sw.count();
    std::cout << "Sequential (sw): " << seq_timing_sw << "\n";

    // Local alignment
    auto timing = function_test_threads(SmithWatermanWavefront, num_threads);
    std::cout << "Regular: ";
    for(auto i: timing){
        std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    }
    std::cout << "\n\n";

    timing = function_test_threads(SmithWatermanWavefrontTp, num_threads);
    std::cout << "Threadpool: ";
    for(auto i: timing){
        std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    }
    std::cout << "\n\n";

    timing = function_test_threads(SmithWatermanWavefront_ScoreOnly, num_threads);
    std::cout << "Score only: ";
    for(auto i: timing){
        std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    }
    std::cout << "\n\n";

    timing = function_test_threads(SmithWatermanWavefront_ScoreOnly_Tp, num_threads);
    std::cout << "Score only + TP: ";
    for(auto i: timing){
        std::cout << i <<" (" << seq_timing_sw/i<< ") |";
    }
    std::cout << "\n\n";

    // Global alignment:
    std::string simiA = generate_random_dna(1<<15);
    std::string simiB = generate_similar_dna(1<<15, 0.9, simiA);

    auto start_seq_nw = std::chrono::high_resolution_clock::now();
    auto nw_dp_threads = needleman_wunsh_dp(simiA, simiB, -1, 1, -2);
    auto nw_alignment_threads = needleman_wunsh_traceback(nw_dp_threads, simiA, simiB, 1, -1, -2);
    int nw_score = 0;
    for (size_t i = 0; i < nw_alignment_threads.first.size(); ++i) {
        char a = nw_alignment_threads.first[i];
        char b = nw_alignment_threads.second[i];
        if (a == '-' || b == '-') nw_score += -2;
        else nw_score += (a == b ? 1 : -1);
    }
    auto end_seq_nw = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_seq_nw = end_seq_nw - start_seq_nw;
    std::cout << "Reference sequential nw time: " << duration_seq_nw.count() << "\n";

    auto start_seq_mm = std::chrono::high_resolution_clock::now();
    auto [mm_score, mm_alignment] = myers_miller(simiA, simiB, 1, -1, -2);
    auto end_seq_mm = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_seq_mm = end_seq_mm - start_seq_mm;
    std::cout << "Reference sequential mm time: " << duration_seq_mm.count() << "\n";

    auto timing_global = function_test_threads(NeedlemanWunschWavefront, num_threads);
    std::cout << "Parallel NW (wrt NW)";
    for(auto i: timing_global){
        std::cout << i <<" (" << duration_seq_nw.count()/i<< ") |";
    }
    std::cout << "\n\n";
    std::cout << "Parallel NW (wrt MM)";
    for(auto i: timing_global){
        std::cout << i <<" (" << duration_seq_mm.count()/i<< ") |";
    }
    std::cout << "\n\n";
    }

    if (mode == "size") {
    //////////////////////////// Size comparison //////////////////////////////////////
    std::vector<int> sizes = {1<<11, 1<<12, 1<<13};

    // Local Alignment
    std::vector<double> seq_timings;
    std::vector<double> gpu_timings;
    for(int i: sizes){
        std::string refA = generate_random_dna(i);
        std::string refB = generate_random_dna(i);  
        
        auto start_seq_sw = std::chrono::high_resolution_clock::now();
        auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
        auto end_seq_sw = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq_timing_temp = end_seq_sw - start_seq_sw;
        seq_timings.push_back(seq_timing_temp.count());
        int res_seq = sw_seq.first[sw_seq.second.first][sw_seq.second.second];

        #ifdef USE_CUDA
        auto start_cuda = std::chrono::high_resolution_clock::now();
        auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
        auto end_cuda = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
        gpu_timings.push_back(duration_cuda.count());
        int res_gpu = std::get<0>(cuda_result);
        if(res_gpu!=res_seq){
            std::cout << "Incorrect results: " << res_seq << " / " << res_gpu << "\n";
        }
        #endif
    }
    std::cout << "Seq (SW): ";
    for(auto i: seq_timings){
        std::cout << i << " ";
    }
    std::cout << "\n\n";

    #ifdef USE_CUDA
    std::cout << "GPU: ";
    for(int i=0; i<gpu_timings.size(); ++i){
        std::cout << " " << gpu_timings[i] << " (" << seq_timings[i]/gpu_timings[i] << ") |";
    }
    std::cout << "\n\n";
    #endif

    auto timing = function_test_size(SmithWatermanWavefront, sizes, sizes, 8);
    std::cout << "Standard parallel: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    timing = function_test_size(SmithWatermanWavefrontTp, sizes, sizes, 8);
    std::cout << "Parallel TP: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    timing = function_test_size(SmithWatermanWavefront_ScoreOnly, sizes, sizes, 8);
    std::cout << "Parallel score only: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    timing = function_test_size(SmithWatermanWavefront_ScoreOnly_Tp, sizes, sizes, 8);
    std::cout << "Parallel score only + Tp: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    // Global Alignment
    std::vector<double> seq_nw_timings;
    std::vector<double> seq_mm_timings;
    for(int i: sizes){
        std::string simiA = generate_random_dna(i);
        std::string simiB = generate_random_dna(i);  
        
        auto start_seq_nw = std::chrono::high_resolution_clock::now();
        auto nw_dp_threads = needleman_wunsh_dp(simiA, simiB, -1, 1, -2);
        auto nw_alignment_threads = needleman_wunsh_traceback(nw_dp_threads, simiA, simiB, 1, -1, -2);
        int nw_score = 0;
        for (size_t i = 0; i < nw_alignment_threads.first.size(); ++i) {
            char a = nw_alignment_threads.first[i];
            char b = nw_alignment_threads.second[i];
            if (a == '-' || b == '-') nw_score += -2;
            else nw_score += (a == b ? 1 : -1);
        }
        auto end_seq_nw = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_seq_nw = end_seq_nw - start_seq_nw;
        seq_nw_timings.push_back(duration_seq_nw.count());

        auto start_seq_mm = std::chrono::high_resolution_clock::now();
        auto [mm_score, mm_alignment] = myers_miller(simiA, simiB, 1, -1, -2);
        auto end_seq_mm = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_seq_mm = end_seq_mm - start_seq_mm;
        seq_mm_timings.push_back(duration_seq_mm.count());

    }
    std::cout << "Seq (NW): ";
    for(auto i: seq_nw_timings){
        std::cout << i << " ";
    }
    std::cout << "\n\n";

    std::cout << "Seq (MM): ";
    for(auto i: seq_mm_timings){
        std::cout << i << " ";
    }
    std::cout << "\n\n";

    auto timing_global = function_test_size(NeedlemanWunschWavefront, sizes, sizes, 8);
    std::cout << "NW parallel wrt NW";
    for(int i=0; i<timing_global.size(); ++i){
        std::cout << " " << timing_global[i] << " (" << seq_nw_timings[i]/timing_global[i] << ") |";
    }
    std::cout << "\n\n";
    std::cout << "NW parallel wrt MM";
    for(int i=0; i<timing_global.size(); ++i){
        std::cout << " " << timing_global[i] << " (" << seq_mm_timings[i]/timing_global[i] << ") |";
    }
    std::cout << "\n\n";
    }

    if (mode == "similarity") {
     ////////////////////// Similarity comparison /////////////////////////////////////////
    std::vector<double> similarities = {0,0.2,0.8};
    std::vector<double> seq_timings;
    int size = 1<<14;
    for(double i: similarities){
        std::string refA = generate_random_dna(size);
        std::string refB = generate_similar_dna(size, i, refA);  
        auto start_seq_sw = std::chrono::high_resolution_clock::now();
        auto sw_seq = smith_waterman_dp(refA, refB, -1, 1, -2);
        auto end_seq_sw = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq_timing_temp = end_seq_sw - start_seq_sw;
        double seq_timing = seq_timing_temp.count();
        seq_timings.push_back(seq_timing);
    }
    std::cout << "Seq (SW): ";
    for(auto i: seq_timings){
        std::cout << i << " ";
    }
    std::cout << "\n\n";

    auto timing = function_test_similarity(SmithWatermanWavefront, similarities, 8);
    std::cout << "Standard parallel: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    timing = function_test_similarity(SmithWatermanWavefrontTp, similarities, 8);
    std::cout << "Parallel TP: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    timing = function_test_similarity(SmithWatermanWavefront_ScoreOnly, similarities, 8);
    std::cout << "Parallel score only: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    timing = function_test_similarity(SmithWatermanWavefront_ScoreOnly_Tp, similarities, 8);
    std::cout << "Parallel score only + Tp: ";
    for(int i=0; i<timing.size(); ++i){
        std::cout << " " << timing[i] << " (" << seq_timings[i]/timing[i] << ") |";
    }
    std::cout << "\n\n";

    std::vector<double> gpu_timings;
    #ifdef USE_CUDA
    for(double i: similarities){
        std::string refA = generate_random_dna(size);
        std::string refB = generate_similar_dna(size, i, refA);  
        auto start_cuda = std::chrono::high_resolution_clock::now();
        auto cuda_result = CUDAlign(refA, refB, -1, 1, -2);
        auto end_cuda = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_cuda = end_cuda - start_cuda;
        gpu_timings.push_back(duration_cuda.count());
    }
    #endif
    std::cout << "GPU: ";
    for(int i=0; i<gpu_timings.size(); ++i){
        std::cout << " " << gpu_timings[i] << " (" << seq_timings[i]/gpu_timings[i] << ") |";
    }
    std::cout << "\n\n";
    }
    if (mode == "memory") {
        size_t length_seq = 14;
        auto refA = generate_random_dna(1<<length_seq);
        auto refB = generate_random_dna(1<<length_seq);
        auto nw_dp = needleman_wunsh_dp(refA, refB, -1, 1, -2);
        print_memory_usage(); 
    }
    return 0;
}