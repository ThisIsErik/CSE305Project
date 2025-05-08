#ifndef LOCAL_MAX_H
#define LOCAL_MAX_H

#include <vector>   // for std::vector
#include <utility>  // for std::pair

struct LocalMax {
    int val = 0;
    int i = 0;
    int j = 0;
};
typedef std::pair<std::vector<std::vector<int>>, std::pair<int, int>> SWResult;

#endif // LOCAL_MAX_H