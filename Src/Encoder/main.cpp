#include <iostream>
#include "../../Inc/gfOps.h"
#include "../../Inc/polyOps.h"  

std::vector<uint32_t> generatorPolynomial(int m, int t);
std::vector<uint32_t> generateCode(int m, int t, const std::vector<uint32_t>& message, const std::vector<uint32_t>& generator);

int main() {
    int m = 8, t = 2;
    
    // Build generator poly
    std::vector<uint32_t> g = generatorPolynomial(m, t);
    std::cout << "Generator poly: ";
    for (auto c : g) std::cout << c << " ";
    std::cout << std::endl;

    // Simple message: [α^0, α^1] = [4, 2] in GF(8)
    std::vector<uint32_t> message = {4U, 2U};
    
    // Encode
    std::vector<uint32_t> codeword = generateCode(m, t, message, g);
    std::cout << "Codeword: ";
    for (auto c : codeword) std::cout << c << " ";
    std::cout << std::endl;

    // Verify: codeword mod g should be all zeros
    std::vector<uint32_t> check = divideGFPoly(m, codeword, g);
    std::cout << "Remainder (should be all zeros): ";
    for (auto c : check) std::cout << c << " ";
    std::cout << std::endl;
}

/**
 * @author: Saleh D.
 * @date: 2026-05-11
 * Construct the generator polynomial for the Reed-Solomon code.
 * The generator polynomial is defined as:
 * g(x) = (x - α^1)(x - α^2)...(x - α^t)
 * where α is a primitive element of the finite field GF(m) and t is the number of errors the code can correct.
 */
std::vector<uint32_t> generatorPolynomial(int m, int t) {
    // Start with constant term
    // Retrieve one in the respective GF(m)
    uint32_t oneTerm = retrieveGFElement(m, 0);
    std::vector<uint32_t> generatorPoly = {oneTerm};
    for (int i = 1; i <= 2*t; i++) {
        uint32_t alphaTerm = retrieveGFElement(m, i);
        std::vector<uint32_t> factor = {oneTerm, alphaTerm};
        generatorPoly = multiplyGFPoly(m, generatorPoly, factor);
    }
    return generatorPoly;
}


std::vector<uint32_t> generateCode(int m, int t, const std::vector<uint32_t>& message, const std::vector<uint32_t>& generator) {
    // Shift message by 2*t positions (multiply by x^(2*t))
    std::vector<uint32_t> shiftedMessage(message.size() + 2*t, 0U);
    for (size_t i = 0; i < message.size(); i++) {
        shiftedMessage[i + 2*t] = message[i]; // shift to high-degree end
    }
    // Calculate the remainder of the division of the shifted message by the generator polynomial
    std::vector<uint32_t> remainder = divideGFPoly(m, shiftedMessage, generator);
    remainder.resize(2*t, 0U);
    // The codeword is the shifted message plus the remainder
    std::vector<uint32_t> codeword = shiftedMessage;
    for (size_t i = 0; i < remainder.size(); i++) {
        codeword[i] ^= remainder[i];
    }
    return codeword;
}