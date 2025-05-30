#include <iostream>
#include <vector>
#include <iomanip>
#include <tuple>
#include "utils/types.h"

int Check_Matrix_Matrix(
    const SWResult& sw_seq,
    const SWResult& sw_par
) {
    const auto& dp_seq = sw_seq.first;
    const auto& dp_par = sw_par.first;
    const auto& pos_seq = sw_seq.second;
    const auto& pos_par = sw_par.second;

    if (dp_seq != dp_par) {
        std::cerr << "Fail: DP matrices different\n";
        for (size_t i = 0; i < dp_seq.size(); ++i) {
            if (dp_seq[i] != dp_par[i]) {
                for (size_t j = 0; j < dp_seq[i].size(); ++j) {
                    if (dp_seq[i][j] != dp_par[i][j]) {
                        std::cerr << "  Mismatch at (" << i << "," << j << "): "
                                  << "Seq=" << dp_seq[i][j] << ", "
                                  << "Par=" << dp_par[i][j] << "\n";
                        return -1;
                    }
                }
            }
        }
        return -1;
    }

    // if (pos_seq != pos_par) {
    //      std::cerr << "Fail: Max positions different\n"
    //                << "  Seq: (" << pos_seq.first << "," << pos_seq.second << ")\n"
    //                << "  Par: (" << pos_par.first << "," << pos_par.second << ")\n";
    //      return -1;
    //  } 
     
    else if (dp_seq[pos_seq.first][pos_seq.second] != dp_par[pos_par.first][pos_par.second]) {
        std::cerr << "Fail: Max values different\n"
                  << "  Seq val: " << dp_seq[pos_seq.first][pos_seq.second] << "\n"
                  << "  Par val: " << dp_par[pos_par.first][pos_par.second] << "\n";
        return -1;
    }
    return 0;
}

int Check_Matrix_Score(
    const SWResult& sw_full,
    const SWResultScore& sw_score_only
) {
    const auto& dp = sw_full.first;
    const auto& pos = sw_full.second;
    int full_score = dp[pos.first][pos.second];

    const int score = std::get<0>(sw_score_only);
    const int pos_i = std::get<1>(sw_score_only);
    const int pos_j = std::get<2>(sw_score_only);

    if (full_score != score) {
        std::cerr << "Fail: Scores differ\n"
                  << "  Full matrix score: " << full_score << "\n"
                  << "  Score-only result: " << score << "\n";
        return -1;
    }

    // if (pos.first != pos_i || pos.second != pos_j) {
    //     std::cerr << "Fail: Max positions differ\n"
    //               << "  Full matrix position: (" << pos.first << ", " << pos.second << ")\n"
    //               << "  Score-only position: (" << pos_i << ", " << pos_j << ")\n";
    //     return -1;
    // }

    return 0;
}


int Check_Score_Score(
    const SWResultScore& result1,
    const SWResultScore& result2
) {

    const int score1 = std::get<0>(result1);
    const int pos_i1 = std::get<1>(result1);
    const int pos_j1 = std::get<2>(result1);

    const int score2 = std::get<0>(result2);
    const int pos_i2 = std::get<1>(result2);
    const int pos_j2 = std::get<2>(result2);

    if (score1 != score2) {
        std::cerr << "Fail: Scores differ\n"
                  << "  Score of first input: " << score1 << "\n"
                  << "  Score of second input: " << score2 << "\n";
        return -1;
    }

    if (pos_i1 != pos_i2 || pos_j1 != pos_j2) {
        std::cerr << "Fail: Max positions differ\n"
                  << "  First input position: (" << pos_i1 << ", " << pos_j1 << ")\n"
                  << "  Second input position: (" << pos_i2 << ", " << pos_j2 << ")\n";
        return -1;
    }

    return 0;
}