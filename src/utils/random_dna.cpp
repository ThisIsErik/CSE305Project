#include "../utils/random_dna.h"
#include <random>

std::string generate_random_dna(size_t length) {
    const std::string bases = "ATCG";
    std::string result;
    result.reserve(length);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 3);

    for (size_t i = 0; i < length; ++i) {
        result += bases[dist(gen)];
    }

    return result;
}
