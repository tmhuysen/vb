#include "bits.hpp"

namespace vb {

size_t next_bitset_permutation(const size_t occupation_num_vec){
    size_t intermediate_vec = occupation_num_vec;

    // t gets occupation_num_vec's least significant 0 bits set to 1
    size_t intermediate_vec_two = intermediate_vec | (intermediate_vec - 1);
    // Next set to 1 the most significant bit to change and set to 0 the least significant ones,
    // and add the necessary 1 bits.
    size_t target_vec = (intermediate_vec_two + 1) | (((~intermediate_vec_two & (intermediate_vec_two+1)) - 1) >> (__builtin_ctzl(intermediate_vec) + 1));

    return  target_vec;
}
size_t smallest_bitset(size_t num_occupation) {
    size_t occupation_num_vec = 1;
    for (size_t i = 1; i < num_occupation; i++) {
        occupation_num_vec <<= 1;
        occupation_num_vec += 1;
    }
    return occupation_num_vec;

}

}