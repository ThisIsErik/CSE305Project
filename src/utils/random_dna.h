#ifndef RANDOM_DNA_H
#define RANDOM_DNA_H

#include <string>

std::string generate_random_dna(size_t length);
std::string generate_similar_dna(size_t length, double similarity, const std::string& reference = "");

#endif
