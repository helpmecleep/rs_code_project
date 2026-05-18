#include <iostream>
#include <vector>
#include <stdexcept>
#include "../../Inc/gfOps.h"
#include "../../Inc/polyOps.h"


std::vector<uint32_t> computeSyndromes(int m, int t, const std::vector<uint32_t>& received);

int main() { 
    return 0;
}


/**
 * This function computes the syndromes for a received codeword in a BCH code.
 * This will compute 2t syndromes, where t is the error correction capability of the code. 
 * The syndromes are computed using the received codeword 
 * and the Galois Field elements corresponding to the powers of alpha.
 */
std::vector<uint32_t> computeSyndromes(int m, int t, const std::vector<uint32_t>& received) {
    std::vector<uint32_t> syndromes(2 * t, 0U);
    for (int s = 1; s <= 2*t; s++) {
        uint32_t syndromeVal = 0U;
        for (size_t i = 0; i < received.size(); i++) {
            int degree = (int)(received.size() - 1 - i); // actual degree of this term
            int exponent = (s * degree) % (m - 1);
            syndromeVal ^= multiplyGFElement(m, received[i], retrieveGFElement(m, exponent));
        }
        syndromes[s - 1] = syndromeVal;
    }
    return syndromes;
}


