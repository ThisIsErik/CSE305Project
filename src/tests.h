#ifndef TESTS_H
#define TESTS_H

#include <vector>
#include <utility>
#include <iomanip>

typedef std::pair<std::vector<std::vector<int>>, std::pair<int, int>> SWResult;

int CheckSWWavefront(const SWResult& sw_seq, const SWResult& sw_par);
int CheckSWWavefront_ScoreOnly(const SWResult& sw_seq, const std::tuple<int, int, int>& sw_par);


template <typename T>
void printMatrix(const std::vector<std::vector<T>>& matrix, 
                 std::pair<int, int> highlight = {-1, -1}) {
    const int cell_width = 5;

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            if (highlight.first == static_cast<int>(i) && highlight.second == static_cast<int>(j)) {
                std::cout << std::setw(cell_width - 1) << "*" << matrix[i][j] << "*";
            } else {
                std::cout << std::setw(cell_width) << matrix[i][j];
            }
        }
        std::cout << "\n";
    }
}



#endif