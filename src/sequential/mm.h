#ifndef MM_H
#define MM_H

#include <string>
#include <utility>

<<<<<<< Updated upstream
<<<<<<< Updated upstream
std::pair<std::string, std::string> myers_miller(
=======
std::tuple<std::string, std::string, int> myers_miller_align(
>>>>>>> Stashed changes
=======
std::tuple<std::string, std::string, int> myers_miller_align(
>>>>>>> Stashed changes
    const std::string& A,
    const std::string& B,
    int match,
    int mismatch,
    int gap);

int myers_miller_score(
    const std::string& A,
    const std::string& B,
    int match,
    int mismatch,
    int gap);

#endif // MM_H
