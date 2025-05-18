#include <algorithm>
#include <iostream>
#include <chrono>
#include <math.h>
#include "CUDAlign.cuh"
#include <cuda_runtime.h>
#include "utils/types.h"

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
    int alpha,
    int* max_val,
    long long* max_pos) {
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

                    int old_max = atomicMax(max_val, val);
                    long long pos_code = (static_cast<long long>(i + 1) << 32) | (j + 1);

                    if (val > old_max) {
                        *max_pos = pos_code;
                    }
                    else if (val == old_max) {
                        atomicMax(reinterpret_cast<unsigned long long*>(max_pos),
                                static_cast<unsigned long long>(pos_code));
                    }
                    
                }
            }
            __syncthreads(); // sync all threads within the block before moving on to the next internal diagonal
        }
    }
}

SWResultScore CUDAlign(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g
) {
    const size_t BLOCKS_NUM = 32;
    const size_t THREADS_PER_BLOCK = 32;
    const size_t ROWS_PER_THREAD = 4;

    const int m = static_cast<int>(A.size());
    const int n = static_cast<int>(B.size());

    int host_max_val = 0;
    long long host_max_pos = 0;

    int* dev_max_val;
    long long* dev_max_pos;
    cudaMalloc(&dev_max_val, sizeof(int));
    cudaMalloc(&dev_max_pos, sizeof(long long));
    cudaMemcpy(dev_max_val, &host_max_val, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_max_pos, &host_max_pos, sizeof(long long), cudaMemcpyHostToDevice);

    const size_t BIG_COLUMNS = BLOCKS_NUM; // B
    const size_t BIG_ROWS = m/(THREADS_PER_BLOCK * ROWS_PER_THREAD); // m/(alpha*T)
    const size_t columns = n/BLOCKS_NUM; // C = n/B
    const size_t rows = THREADS_PER_BLOCK * ROWS_PER_THREAD; // R = alpha*T
    const size_t alpha = ROWS_PER_THREAD; // alpha

    // moving the data to device 
    int* dpd;
    cudaMalloc(&dpd, (m + 1) * (n + 1) * sizeof(int));
    cudaMemset(dpd, 0, (m + 1) * (n + 1) * sizeof(int));

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
            alpha,
            dev_max_val,
            dev_max_pos
        );
    }

    // copying the result back
    cudaMemcpy(&host_max_val, dev_max_val, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_max_pos, dev_max_pos, sizeof(long long), cudaMemcpyDeviceToHost);
    int max_i = static_cast<int>(host_max_pos >> 32);
    int max_j = static_cast<int>(host_max_pos & 0xFFFFFFFF);
    // std::cout << "Max value: " << host_max_val << " at (" << max_i << ", " << max_j << ")\n";

    cudaFree(dpd);
    cudaFree(A_dev);
    cudaFree(B_dev);
    cudaFree(dev_max_val);
    cudaFree(dev_max_pos);

    return {host_max_val, max_i, max_j};
}
