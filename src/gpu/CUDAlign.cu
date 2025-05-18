#include <algorithm>
#include <iostream>
#include <chrono>
#include <math.h>
#include "CUDAlign.cuh"
#include <cuda_runtime.h>

__device__ __host__ inline int get_index(int i, int j, int n) {
    return i * (n + 1) + j;
}

__global__
void CUDAlignAux(int* dp, 
    const char* A,
    const char* B,
    int external_diagonal, 
    int mi, int ma, int g, 
    int rows, int columns, 
    int BIG_ROWS, int BIG_COLUMNS,
    int m, int n,
    int alpha) {
    if (blockIdx.x <= external_diagonal && external_diagonal - blockIdx.x < BIG_ROWS){
        // in big grid, this block will treat the big cell indexed G_(external_diagonal-blockIdx, blockIdx)
        int topLeft_i = (external_diagonal-blockIdx.x)*rows;
        int topLeft_j = blockIdx.x * columns;
        for (int internal_diagonal = 0; internal_diagonal < rows/alpha + columns - 1; internal_diagonal++){
            // inside the big gird cell, we have rows x columns individual cells
            // our thread will take care of lines threadIdx*alpha to (threadIdx+1)*alpha - 1

            // on this internal diagonal, the current thread will take care of alpha rows on this column
            int row_thread = threadIdx.x * alpha;
            int column_thread = internal_diagonal - threadIdx.x; 
            if (column_thread >= 0 && column_thread < columns){
                // check column index is valid

                for (int offset = 0; offset < alpha; offset++){
                    // offset will say which one out of the alpha rows we process
                    int i = topLeft_i + row_thread + offset; // i = topLeft_i index of big gridcell + which rows these thread processes + offset
                    int j = topLeft_j + column_thread; // j = topLeft_j index of big gridcell + which column this thread processes for this antidiagonal

                    int p = (A[i] == B[j]) ? ma : mi;
                    int val = max(max(max(
                    dp[get_index(i, j, n)] + p,
                    dp[get_index(i + 1, j, n)] + g),
                    dp[get_index(i, j + 1, n)] + g),
                    0);
                    dp[get_index(i + 1, j + 1, n)] = val;
                }
            }
            __syncthreads(); // sync all threads within the block before moving on to the next internal diagonal
        }
    }
}

void CUDAlign(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g
) {
    const size_t BLOCKS_NUM = 32;
    const size_t THREADS_PER_BLOCK = 32;
    const size_t ROWS_PER_THREAD = 4;

    const int m = static_cast<int>(A.size());
    const int n = static_cast<int>(B.size());
    int* dp = new int[(m + 1) * (n + 1)]();

    // int max_val = 0;
    // std::pair<int, int> max_pos = {0, 0};

    const size_t BIG_COLUMNS = BLOCKS_NUM; // B
    const size_t BIG_ROWS = m/(THREADS_PER_BLOCK * ROWS_PER_THREAD); // m/(alpha*T)
    const size_t columns = n/BLOCKS_NUM; // C = n/B
    const size_t rows = THREADS_PER_BLOCK * ROWS_PER_THREAD; // R = alpha*T
    const size_t alpha = ROWS_PER_THREAD; // alpha

    // moving the data to device 
    int* dpd;
    cudaMalloc(&dpd, (m + 1) * (n + 1) * sizeof(int));
    cudaMemcpy(dpd, dp, (m + 1) * (n + 1) * sizeof(int), cudaMemcpyHostToDevice);

    char* A_dev;
    char* B_dev;
    cudaMalloc(&A_dev, m * sizeof(char));
    cudaMalloc(&B_dev, n * sizeof(char));
    cudaMemcpy(A_dev, A.c_str(), m * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(B_dev, B.c_str(), n * sizeof(char), cudaMemcpyHostToDevice);

    // computing on GPU
    for (int external_diagonal = 0; external_diagonal < BIG_COLUMNS + BIG_ROWS - 1; external_diagonal++) {
        CUDAlignAux<<<BLOCKS_NUM, THREADS_PER_BLOCK>>>(
            dpd,
            A_dev, B_dev,
            external_diagonal,
            mi, ma, g,
            rows, columns,
            BIG_ROWS, BIG_COLUMNS,
            m, n,
            alpha
        );
    }

    // copying the result back
    cudaMemcpy(dp, dpd, (m + 1) * (n + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(dpd);
    cudaFree(A_dev);
    cudaFree(B_dev);
}

// //----------------------------------------------------

// int main(int argc, char* argv[]) {
//     // setting the random seed to get the same result each time
//     srand(42);

//     // taking as input, which algo to run
//     int alg_ind = std::stoi(argv[1]);

//     // Generating data
//     size_t N = 1 << 27;
//     double* x = (double*) malloc(N * sizeof(double));
//     double* y = (double*) malloc(N * sizeof(double));
//     for (size_t i = 0; i < N; ++i) {
//           x[i] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
//           y[i] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
//     }
 
//     // Warming up the kernel
//     auto start = std::chrono::steady_clock::now();
//     warmup<<<1,1>>>();
//     auto finish = std::chrono::steady_clock::now();
//     auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count(); 
//     std::cout << "Warm-up time: " << elapsed << std::endl;
 

//     // Allocating the result
//     double* result = (double*) malloc(N * sizeof(double));
//     start = std::chrono::steady_clock::now();
//     switch (alg_ind) {
//         case 0: 
//             Add(x, y, result, N);
//             break;
//         case 1:
//             AddGPU(x, y, result, N);
//             break;
//         case 2:
//             AddGPU2(x, y, result, N);
//             break;
//     }
//     finish = std::chrono::steady_clock::now();
//     elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count(); 
//     std::cout << "Elapsed time: " << elapsed << std::endl;
    
//     delete[] x;
//     delete[] y;
//     delete[] result;
//     return 0;
// }
