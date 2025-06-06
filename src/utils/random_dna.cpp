#include "../utils/random_dna.h"
#include <random>
#include <iostream>

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


std::string generate_similar_dna(size_t length, double similarity, const std::string& reference) {
   static const char nucleotides[] = {'A', 'C', 'G', 'T'};
   static std::random_device rd;
   static std::mt19937 gen(rd());
   static std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
   static std::uniform_int_distribution<> base_dist(0, 3);
   static std::uniform_int_distribution<> diff_dist(1, 3);
  
   std::string sequence; //the new sequecne we are making
   sequence.reserve(length);
  
    if (reference.empty() || similarity == 0.0) {
        for (size_t i = 0; i < length; ++i) {
            sequence.push_back(nucleotides[base_dist(gen)]);
        }
        return sequence;
    }   

   //If the ref is longer, then want some random indices where the similarity is satisfied
   std::vector<size_t> indices;
   if (reference.length() > length) {
       indices.reserve(reference.length());
       for (size_t i = 0; i < reference.length(); ++i) {
           indices.push_back(i);
       }
       std::shuffle(indices.begin(), indices.end(), gen); //mix up all of the indices
       indices.resize(length); //choose the ones you will use for the new sequence
       std::sort(indices.begin(), indices.end());  // keep the original order: now stores the indices which we will be similar to
       for(size_t i = 0; i<length; ++i){
           if(prob_dist(gen) < similarity){
               sequence.push_back(reference[indices[i]]);
           }
           else{
               size_t temp_index = 0;
                   while(nucleotides[temp_index]!=reference[indices[i]]) temp_index++;
                   sequence.push_back(nucleotides[(temp_index+diff_dist(gen))%4]);
           }
       }
   }


   else { //if the ref is shorter then it should be the opposite
       indices.reserve(length);
       for (size_t i = 0; i < length; ++i) {
           indices.push_back(i);
       }
       std::shuffle(indices.begin(), indices.end(), gen);
       indices.resize(reference.length());
       std::sort(indices.begin(), indices.end());
       int curr_index = 0;
       for(size_t i = 0; i<length; ++i){
           if(i==indices[curr_index]){ //need to do it (this is a copy index)
               if(prob_dist(gen) < similarity){
                   sequence.push_back(reference[curr_index]);
               }
               else{
                   size_t temp_index = 0;
                   while(nucleotides[temp_index]!=reference[curr_index]) temp_index++;
                   sequence.push_back(nucleotides[(temp_index+diff_dist(gen))%4]);
               }
               curr_index++;
           }
           else{
               sequence.push_back(nucleotides[base_dist(gen)]);
           }
       }
   }
   return sequence;
}
