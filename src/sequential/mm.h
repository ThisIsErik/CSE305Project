#ifndef MM_H
#define MM_H

#include <vector>
#include <string>
#include <utility>
#include "nw.h"

std::pair<std::string, std::string> myers_miller_align(
    const std::string& A,
    const std::string& B,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend);

#endif // MM_H
