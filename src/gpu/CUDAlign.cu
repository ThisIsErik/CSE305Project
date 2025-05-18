#include <algorithm>
#include <iostream>
#include <chrono>
#include <math.h>
#include "CUDAlign.cuh"
#include <cuda_runtime.h>


__global__
void CUDAlignAux(int* dp, int external_diagonal, int mi, int ma, int g) {
    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
}

void CUDAlign(
    const std::string& A,
    const std::string& B,
    int mi, int ma, int g
) {
    const size_t BLOCKS_NUM = 32;
    const size_t THREADS_PER_BLOCK = 32;
    const size_t ROWS_PER_THREAD = 4;

    const int m = A.size();
    const int n = B.size();
    int* dp = new int[(m + 1) * (n + 1)]();  

    int max_val = 0;
    std::pair<int, int> max_pos = {0, 0};

    const size_t BIG_COLUMNS = BLOCKS_NUM;
    const size_t BIG_ROWS = m/(THREADS_PER_BLOCK * ROWS_PER_THREAD);

    // moving the data to device 
    int* dpd;
    cudaMalloc(&dpd, (m + 1) * (n + 1) * sizeof(int));
    cudaMemcpy(dpd, dp, (m + 1) * (n + 1) * sizeof(int), cudaMemcpyHostToDevice);

    // computing on GPU
    for (int external_diagonal=0;external_diagonal<BIG_COLUMNS + BIG_ROWS - 1; external_diagonal++){
        CUDAlignAux<<<BLOCKS_NUM, THREADS_PER_BLOCK>>>(dpd, external_diagonal, mi, ma, g);
    }

    // copying the result back
    cudaMemcpy(dpd, dp, (m + 1) * (n + 1) * sizeof(int), cudaMemcpyDeviceToHost);
  
    // Free memory
    cudaFree(dpd);
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
