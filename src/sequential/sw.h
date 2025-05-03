#ifndef SW_H
#define SW_H

#include <string>

std::pair<std::vector<std::vector<int>>, std::pair<int, int>> smith_waterman_dp(const std::string& A, const std::string& B, int mi, int ma, int g);
std::pair<std::string, std::string> smith_waterman_traceback(
    const std::vector<std::vector<int>>& dp,
    const std::string& A,
    const std::string& B,
    int ma, int mi, int g,
    std::pair<int, int> start);

#endif
