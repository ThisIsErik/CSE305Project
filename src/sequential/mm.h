#ifndef MM_H
#define MM_H

#include <string>
#include <utility>

std::pair<std::string, std::string> myers_miller(
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
