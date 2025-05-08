#ifndef SW_WF_H
#define SW_WF_H

#include <vector>
#include <string>
#include <utility>
#include <mutex>

typedef std::pair<std::vector<std::vector<int>>, std::pair<int, int>> SWResult;

struct LocalMax {
    int val = 0;
    int i = 0;
    int j = 0;
};

void AntiDiagonalAux(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    int diagonal,
    int start_i,
    int end_i,
    std::vector<std::vector<int>>& dp,
    LocalMax& local_max
);

SWResult SmithWatermanWavefront(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g,
    size_t num_threads
);

#endif 