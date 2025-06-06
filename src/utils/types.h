#ifndef LOCAL_MAX_H
#define LOCAL_MAX_H

#include <vector>   // for std::vector
#include <utility>  // for std::pair
#include <string>

struct LocalMax {
    int val = 0;
    int i = 0;
    int j = 0;
};
typedef std::pair<std::vector<std::vector<int>>, std::pair<int, int>> SWResult;
typedef std::tuple<int, int, int> SWResultScore;

typedef std::tuple<int, int, int> MMResultScore;
typedef std::pair<int, std::pair<std::string, std::string>> AlignmentResult;

constexpr int THRESHOLD = 2000; //TUNE HERE 

constexpr int SCORE_THRESHOLD = 100; // For general scoring
constexpr int PARALLEL_THRESHOLD = 2000; // For parallel processing decisions

#endif // LOCAL_MAX_H