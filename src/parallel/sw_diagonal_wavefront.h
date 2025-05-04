#ifndef SW_WF_H
#define SW_WF_H

#include <vector>
#include <string>
#include <utility>
#include <mutex>

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
);

SWResult SmithWatermanWavefront(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
);

#endif 