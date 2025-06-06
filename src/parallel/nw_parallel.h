#ifndef NW_PARALLEL_H
#define NW_PARALLEL_H

#include "../utils/types.h"
#include <string>
#include <vector>
#include <utility>


AlignmentResult parallel_needleman_wunsch(const std::string& A, const std::string& B, 
                                        int mi, int ma, int g);


std::vector<std::vector<std::pair<int, int>>> get_antidiagonal_indices(int m, int n);


void compute_diagonal_parallel(std::vector<std::vector<int>>& dp,
                             const std::vector<std::pair<int, int>>& diagonal,
                             const std::string& A, const std::string& B,
                             int mi, int ma, int g);

std::pair<std::string, std::string> parallel_traceback(const std::vector<std::vector<int>>& dp,
                                                     const std::string& A, const std::string& B,
                                                     int ma, int mi, int g);

#endif // PARALLEL_NW_H