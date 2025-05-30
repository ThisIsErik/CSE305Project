#ifndef TESTS_H
#define TESTS_H

#include <vector>
#include <utility>
#include <iomanip>
#include "utils/types.h"

int Check_Matrix_Matrix(const SWResult& sw_seq, const SWResult& sw_par);
int Check_Matrix_Score(const SWResult& sw_seq, const SWResultScore& sw_score_only);
int Check_Score_Score(const SWResultScore& result1, const SWResultScore& result2);

template <typename T>
void printMatrix(const std::vector<std::vector<T>>& matrix, 
                 std::pair<int, int> highlight = {-1, -1}) {
    const int cell_width = 3;

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


void printAntidiagonals(const std::vector<std::vector<int>>& matrix) {
    int totalRows = matrix.size();
    int totalCols = matrix[0].size();
    int rows = totalRows - 1;
    int cols = totalCols - 1;

    for (int s = 0; s < rows + cols - 1; ++s) {
        std::vector<int> antidiagonal;
        for (int i = 0; i < rows; ++i) {
            int j = s - i;
            if (j >= 0 && j < cols) {
                antidiagonal.push_back(matrix[i + 1][j + 1]);  
            }
        }

        for (int val : antidiagonal) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
}



#endif