#ifndef NW_H
#define NW_H

#include <vector>
#include <string>

std::vector<std::vector<int>> needleman_wunsh_dp(const std::string& seqA, const std::string& seqB, int mi, int ma, int g);
std::pair<std::string, std::string> needleman_wunsh_traceback(const std::vector<std::vector<int>>& dp, const std::string& A, const std::string& B, int ma, int mi, int g);

#endif
