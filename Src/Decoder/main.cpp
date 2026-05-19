#include <iostream>
#include <vector>
#include <stdexcept>
#include "../../Inc/gfOps.h"
#include "../../Inc/polyOps.h"

struct ErrorPolys {
    // lambda is the error locator 
    std::vector<uint32_t> lambda;
    // omega is the error evaluator
    std::vector<uint32_t> omega;
};

std::vector<uint32_t> computeSyndromes(int m, int t, const std::vector<uint32_t>& received);
ErrorPolys computeErrorPolys(int m, int t, const std::vector<uint32_t>& syndromes);
std::vector<int> chienSearch (int m, int t, int n, const std::vector<uint32_t>& lambda);

int main() {
    int m = 8, t = 2;
    
    // Valid codeword from encoder
    std::vector<uint32_t> received = {5U, 6U, 7U, 3U, 4U, 2U};
    received[2] ^= 3U; // inject error at position 2
    std::cout << "Received: ";
    for (auto c : received) std::cout << c << " ";
    std::cout << std::endl;
    
    // Compute syndromes
    std::vector<uint32_t> syndromes = computeSyndromes(m, t, received);
    std::cout << "Syndromes: ";
    for (auto s : syndromes) std::cout << s << " ";
    std::cout << std::endl;
    
    // Get lambda and omega
    ErrorPolys ep = computeErrorPolys(m, t, syndromes);
    std::cout << "Lambda: ";
    for (auto c : ep.lambda) std::cout << c << " ";
    std::cout << std::endl;
    
    // Chien search
    std::vector<int> errorPositions = chienSearch(m, t, received.size(), ep.lambda);
    std::cout << "Error positions (expected 2): ";
    for (auto p : errorPositions) std::cout << p << " ";
    std::cout << std::endl;
    
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
    for (int s = 1; s <= 2 * t; s++) {
        uint32_t syndromeVal = 0U;
        for (size_t i = 0; i < received.size(); i++) {
            int degree = (int)(received.size() - 1 - i); // actual degree of this term
            int exponent = (s * degree) % (m - 1);
            syndromeVal ^= multiplyGFElement(m, received[i], retrieveGFElement(m, exponent));
        }
        syndromes[s-1] = syndromeVal;
    }
    return syndromes;
}



ErrorPolys computeErrorPolys(int m, int t, const std::vector<uint32_t>& syndromes) {
    std::vector<uint32_t> syndromePoly(2 * t + 1, 0U);
    for (int i = 0; i <= 2 * t; i++) {
        syndromePoly[i + 1] = syndromes[i];
    }
    std::vector<uint32_t> x2t(2 * t + 1, 0U);
    x2t[2*t] = static_cast<uint32_t>(retrieveGFElement(m, 0));
    GCDStruct result = extendedGCDPoly(m, t, x2t, syndromePoly);
    return {result.y, result.x};
}

std::vector<int> chienSearch(int m, int t, int n, const std::vector<uint32_t>& lambda) {
    std::vector<int> errorPositions;
    for (int i = 0; i < (m - 1); i++) {
        uint32_t result = 0U;
        for (int j = 0; j < lambda.size(); j++) {
            int exponent = ((m - 1 - i) * j) % (m - 1);
            result ^= multiplyGFElement(m, lambda[j], retrieveGFElement(m, exponent));
        }
        if (result == 0U) errorPositions.push_back(n - 1 - (int)i);
    }
    return errorPositions;
}