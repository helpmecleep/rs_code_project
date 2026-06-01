#include <iostream>
#include <vector>
#include <stdexcept>
#include "../../Inc/gfOps.h"
#include "../../Inc/polyOps.h"


// Struct to hold omega and lambda 
struct ErrorPolys {
    // lambda is the error locator 
    std::vector<uint32_t> lambda;
    // omega is the error evaluator
    std::vector<uint32_t> omega;
};

// Forward declaration of functions
std::vector<uint32_t> computeSyndromes(int m, int t, const std::vector<uint32_t>& received);
ErrorPolys computeErrorPolys(int m, int t, const std::vector<uint32_t>& syndromes);
std::vector<int> findRoots (int m, int t, int n, const std::vector<uint32_t>& lambda);
std::vector<uint32_t> findErrorMag(int m, int n,
    const std::vector<uint32_t>& omega,const std::vector<uint32_t>& lambda,
    const std::vector<int>& errorPositions);
std::vector<uint32_t> decodeRS(int m, int t, const std::vector<uint32_t>& received);



int main() {
    // Test parameters
    int m = 8, t = 2;
    std::vector<uint32_t> original = {5U, 6U, 7U, 3U, 4U, 2U};
    std::vector<uint32_t> received = original;
    
    // inject error at position 2
    received[2] ^= 5U;

    // Print original and received codewords
    std::cout << "Original: ";
    for (auto c : original) std::cout << c << " ";
    std::cout << std::endl;

    // Print received codeword
    std::cout << "Received: ";
    for (auto c : received) std::cout << c << " ";
    std::cout << std::endl;
    
    // Compute syndromes and error locations and magnitudes
    std::vector<uint32_t> syndromes = computeSyndromes(m, t, received);
    ErrorPolys ep = computeErrorPolys(m, t, syndromes);
    std::vector<int> errorPositions = findRoots(m, t, received.size(), ep.lambda);
    std::vector<uint32_t> errorMags = findErrorMag(m, received.size(), ep.omega, ep.lambda, errorPositions);

    // Print error positions and magnitudes
    std::cout << "Error positions: ";
    for (auto p : errorPositions) std::cout << p << " ";
    std::cout << std::endl;
    std::cout << "Error magnitudes: ";
    for (auto v : errorMags) std::cout << v << " ";
    std::cout << std::endl;

    // Print corrected codeword
    std::vector<uint32_t> corrected = decodeRS(m, t, received);
    std::cout << "Corrected: ";
    for (auto c : corrected) std::cout << c << " ";
    std::cout << std::endl;
    
    return 0;
}

/**
 * @author Saleh D.
 * @date 2026-05-20
 * @brief This function computes the syndromes for a received codeword in a BCH code.
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



/**
 * @author Saleh D.
 * @date 2026-05-20
 * @brief This function computes the error location and magnitude polynomials (lambda and omega)
 * using the Extended Euclidean Algorithm.
 * It takes the syndrome polynomial and x^(2t) as inputs and computes the GCD, which gives omega,
 * 
 */
ErrorPolys computeErrorPolys(int m, int t, const std::vector<uint32_t>& syndromes) {
    std::vector<uint32_t> syndromePoly(2 * t, 0U);
    for (int i = 0; i < 2 * t; i++) {
        syndromePoly[i] = syndromes[i];
    }
    // compute the GCD and retrieve lambda and omega
    std::vector<uint32_t> x2t(2 * t + 1, 0U);
    x2t[2*t] = static_cast<uint32_t>(retrieveGFElement(m, 0));
    GCDStruct result = extendedGCDPoly(m, t, x2t, syndromePoly);
    std::vector<uint32_t> lambda = result.y;
    std::vector<uint32_t> omega = result.gcd;

    // Normalize so lambda[0] == 1 when possible.
    if (!lambda.empty() && lambda[0] != 0U && lambda[0] != 1U) {
        uint32_t inv = inverseGFElement(m, lambda[0]);
        for (size_t i = 0; i < lambda.size(); ++i) {
            lambda[i] = multiplyGFElement(m, lambda[i], inv);
        }
        for (size_t i = 0; i < omega.size(); ++i) {
            omega[i] = multiplyGFElement(m, omega[i], inv);
        }
    }

    return {lambda, omega};
}


/**
 * @author Saleh D.
 * @date 2026-05-20
 * @brief This function performs the Chien Search to find the error positions.
 * It evaluates the error locator polynomial lambda at the inverses of the field elements.
 * If lambda evaluates to zero at a particular field element, it indicates an error at that position.
 * The positions are returned in terms of the original codeword indexing.   
 */
std::vector<int> findRoots(int m, int t, int n, const std::vector<uint32_t>& lambda) {
    // set up error position vectors
    std::vector<int> errorPositions;
    for (int i = 0; i < (m - 1); i++) {
        // evaluate lambda at alpha^(-i)
        uint32_t result = 0U;
        for (int j = 0; j < lambda.size(); j++) {
            int exponent = ((m - 1 - i) * j) % (m - 1);
            result ^= multiplyGFElement(m, lambda[j], retrieveGFElement(m, exponent));
        }
        if (result == 0U) errorPositions.push_back(n - 1 - (int)i);
    }
    return errorPositions;
}
/**
 * @author Saleh D.
 * @date 2026-05-20
 * @brief This function computes the error magnitudes using Forney's formula.
 */
std::vector<uint32_t> findErrorMag(int m, int n,
    const std::vector<uint32_t>& omega,const std::vector<uint32_t>& lambda,
    const std::vector<int>& errorPositions) {

    // Set up error magnitude and lambda derivative vectors
    std::vector<uint32_t> errorMags(errorPositions.size(), 0U);
    std::vector<uint32_t> lambdaDeriv = polyDerivative(lambda);

    // loop through the error positions and compute the error magnitude for each using Forney's formula
    for (size_t k = 0; k < errorPositions.size(); ++k) {
        int pos = errorPositions[k];
        int i = (n - 1) - pos;
        int xExp = i % (m - 1);
        uint32_t Xinv = retrieveGFElement(m, xExp);

        // evaluate omega and lambda derivative at Xinv
        uint32_t omegaEval = 0U;
        for (size_t j = 0; j < omega.size(); ++j) {
            int exponent = (xExp * (int)j) % (m - 1);
            omegaEval ^= multiplyGFElement(m, omega[j], retrieveGFElement(m, exponent));
        }

        uint32_t lambdaDerivEval = 0U;
        for (size_t j = 0; j < lambdaDeriv.size(); ++j) {
            int exponent = (xExp * (int)j) % (m - 1);
            lambdaDerivEval ^= multiplyGFElement(m, lambdaDeriv[j], retrieveGFElement(m, exponent));
        }

        // If lambda derivative evaluates to zer set error magnitude to 0
        if (lambdaDerivEval == 0U) {
            errorMags[k] = 0U;
            continue;
        }
        errorMags[k] = divideGFElement(m, omegaEval, lambdaDerivEval);
    }
    return errorMags;
}

/**
 * @author Saleh D.
 * @date 2026-05-20
 * @brief This function calls the respective other functions to perform the full Reed-Solomon decoding process.
 * It computes the syndromes, error locator and evaluator polynomials, error positions, and error magnitudes, 
 * and then applies the corrections to the received codeword to produce the corrected codeword.
 */
std::vector<uint32_t> decodeRS(int m, int t, const std::vector<uint32_t>& received) {
    std::vector<uint32_t> syndromes = computeSyndromes(m, t, received);
    ErrorPolys ep = computeErrorPolys(m, t, syndromes);
    std::vector<int> errorPositions = findRoots(m, t, received.size(), ep.lambda);
    std::vector<uint32_t> errorMags = findErrorMag(m, received.size(), ep.omega, ep.lambda, errorPositions);

    std::vector<uint32_t> corrected = received;
    for (size_t i = 0; i < errorPositions.size(); ++i) {
        int pos = errorPositions[i];
        corrected[pos] ^= errorMags[i];
    }
    return corrected;
}