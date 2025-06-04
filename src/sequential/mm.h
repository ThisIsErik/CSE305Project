// mm.h
#ifndef MM_H
#define MM_H

#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <algorithm>
#include <limits>
#include "../utils/types.h"
#include "nw.h"

std::pair<int, std::pair<std::string, std::string>> myers_miller(
    const std::string& A, const std::string& B, int ma, int mi, int g);

#endif // MM_H