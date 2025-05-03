// #ifndef TIMING_UTILS_H
// #define TIMING_UTILS_H

// #include <chrono>
// #include <iostream>

// // Measures and prints how long a function takes to run
// inline int time_levenshtein(const std::string& A, const std::string& B, int (*algo)(const std::string&, const std::string&)) {
//     auto start = std::chrono::high_resolution_clock::now();

//     int result = algo(A, B);  // run algorithm

//     auto end = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration<double>(end - start);

//     std::cout << "Time taken: " << duration.count() << " seconds\n";
//     return result;
// }

// #endif
