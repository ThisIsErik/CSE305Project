#include <algorithm>
#include <iostream>
#include <chrono>
#include <math.h>
#include "CUDAlign.cuh"
#include <cuda_runtime.h>
#include "utils/types.h"

__device__ __host__ inline size_t get_index(size_t i, size_t j, size_t n) {
    return i * (n + 1ULL) + j;
}

__global__
void CUDAlignAux(unsigned short* dp, 
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
        int topLeft_i = (external_diagonal - blockIdx.x) * rows;
        int topLeft_j = blockIdx.x * columns;

        for (int internal_diagonal = 0; internal_diagonal < rows / alpha + columns - 1; internal_diagonal++) {
            int row_thread = threadIdx.x * alpha;
            int column_thread = internal_diagonal - threadIdx.x;

            if (column_thread >= 0 && column_thread < columns) {
                for (int offset = 0; offset < alpha; offset++) {
                    int i = topLeft_i + row_thread + offset;
                    int j = topLeft_j + column_thread;

                    if (i < 0 || j < 0 || i >= m || j >= n || i + 1 > m || j + 1 > n) {
                        printf("[OOB] Thread (%d,%d): i=%d j=%d (m=%d, n=%d)\n", blockIdx.x, threadIdx.x, i, j, m, n);
                        continue;
                    }

                    int p = (A[i] == B[j]) ? ma : mi;

                    size_t write_index = ((size_t)(i + 1)) * (n + 1ULL) + (j + 1);
                    if (write_index >= (size_t)(m + 1) * (n + 1)) {
                        printf("[ERROR] Invalid dp index at (%d,%d): index=%zu\n", i + 1, j + 1, write_index);
                        continue;
                    }

                    if (i == 0 && j == 0 && threadIdx.x == 0 && blockIdx.x == 0) {
                        printf("[TRACE] Kernel start: A[0]=%c B[0]=%c n=%d\n", A[0], B[0], n);
                    }

                    unsigned short val = max(max(max(
                        (int)dp[get_index(i, j, n)] + p,
                        (int)dp[get_index(i + 1, j, n)] + g),
                        (int)dp[get_index(i, j + 1, n)] + g),
                        0);

                    dp[write_index] = val;

                    int old_max = atomicMax(max_val, static_cast<int>(val));
                    long long pos_code = (static_cast<long long>(i + 1) << 32) | (j + 1);

                    if (val > old_max) {
                        *max_pos = pos_code;
                    } else if (val == old_max) {
                        atomicMax(reinterpret_cast<unsigned long long*>(max_pos),
                                  static_cast<unsigned long long>(pos_code));
                    }
                }
            }
            __syncthreads();
        }
    }
}

SWResultScore CUDAlign(const std::string& A, const std::string& B, int mi, int ma, int g) {
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

    const size_t rows = THREADS_PER_BLOCK * ROWS_PER_THREAD;
    const size_t columns = (n + BLOCKS_NUM - 1) / BLOCKS_NUM;
    const size_t BIG_COLUMNS = BLOCKS_NUM;
    const size_t BIG_ROWS = (m + rows - 1) / rows;
    const size_t alpha = ROWS_PER_THREAD;

    unsigned short* dpd;
    cudaError_t err = cudaMalloc(&dpd, (m + 1) * (n + 1) * sizeof(unsigned short));
    if (err != cudaSuccess) {
        std::cerr << "CUDA malloc failed: " << cudaGetErrorString(err) << std::endl;
        return {0, 0, 0};
    }
    cudaMemset(dpd, 0, (m + 1) * (n + 1) * sizeof(unsigned short));

    char* A_dev;
    char* B_dev;
    cudaMalloc(&A_dev, m * sizeof(char));
    cudaMalloc(&B_dev, n * sizeof(char));
    cudaMemcpy(A_dev, A.c_str(), m * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(B_dev, B.c_str(), n * sizeof(char), cudaMemcpyHostToDevice);

    for (int external_diagonal = 0; external_diagonal < BIG_COLUMNS + BIG_ROWS - 1; external_diagonal++) {
        CUDAlignAux<<<BLOCKS_NUM, THREADS_PER_BLOCK>>>(
            dpd, A_dev, B_dev, external_diagonal,
            mi, ma, g,
            rows, columns,
            BIG_ROWS, BIG_COLUMNS,
            m, n, alpha,
            dev_max_val, dev_max_pos);

        cudaError_t kernel_err = cudaGetLastError();
        if (kernel_err != cudaSuccess) {
            std::cerr << "Kernel launch failed: " << cudaGetErrorString(kernel_err) << std::endl;
            break;
        }
        cudaDeviceSynchronize();
    }

    cudaMemcpy(&host_max_val, dev_max_val, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_max_pos, dev_max_pos, sizeof(long long), cudaMemcpyDeviceToHost);
    int max_i = static_cast<int>(host_max_pos >> 32);
    int max_j = static_cast<int>(host_max_pos & 0xFFFFFFFF);

    cudaFree(dpd);
    cudaFree(A_dev);
    cudaFree(B_dev);
    cudaFree(dev_max_val);
    cudaFree(dev_max_pos);

    return {host_max_val, max_i, max_j};
}
