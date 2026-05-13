#pragma once
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include "gfOps.h"


uint32_t multiplyPolyElement(uint32_t multiplicand, uint32_t multiplier) {
    if (multiplicand == 0U || multiplier == 0U) return 0U;
    uint32_t product = 0U;
    while (multiplier != 0U) {
        if (multiplier & 1U) {
            product ^= multiplicand;
        }
        multiplicand <<= 1U;
        multiplier >>= 1U;
    }
    return product;
}

int polyDegree(uint32_t p) {
    int degree = 0;
    while (p >>= 1U) {
        degree++;
    }
    return degree;
}

uint32_t polyMod(uint32_t dividend, uint32_t divisor) {
    if (divisor == 0U) throw std::invalid_argument("Divisor cannot be zero");
    uint32_t remainder = dividend;
    
    int degDivisor = polyDegree(divisor);
    int degRemainder = polyDegree(remainder);

    while (remainder != 0U && degRemainder >= degDivisor) {
        int shift = degRemainder - degDivisor;
        remainder ^= (divisor << shift);
        degRemainder = polyDegree(remainder);
    }
    return remainder;
}

std::vector<uint32_t> multiplyGFPoly(int m, const std::vector<uint32_t>& multiplicand, const std::vector<uint32_t>& multiplier) {
    std::vector<uint32_t> result(multiplicand.size() + multiplier.size() - 1, 0U);

    for (size_t i = 0; i < multiplicand.size(); i++) {
        for (size_t j = 0; j < multiplier.size(); j++) {
            result[i + j] ^= multiplyGFElement(m, multiplicand[i], multiplier[j]);
        }
    }
    
    return result;
}


std::vector<uint32_t> divideGFPoly (int m, const std::vector<uint32_t>& dividend, const std::vector<uint32_t>& divisor) {
    if (divisor.empty() || (divisor.size() == 1 && divisor[0] == 0U)) {
        throw std::invalid_argument("Divisor cannot be zero polynomial");
    }
    
    std::vector<uint32_t> quotient(dividend.size(), 0U);
    std::vector<uint32_t> remainder = dividend;

    int degDivisor = divisor.size() - 1;
    int degRemainder = remainder.size() - 1;

    while (degRemainder >= degDivisor) {
        int shift = degRemainder - degDivisor;
        quotient[shift] = divideGFElement(m, remainder[degRemainder], divisor[degDivisor]);        
        for (size_t i = 0; i < divisor.size(); i++) {
            remainder[shift + i] ^= multiplyGFElement(m, quotient[shift], divisor[i]);
        }
        degRemainder = remainder.size() - 1;
        while (degRemainder >= 0 && remainder[degRemainder] == 0U) {
            degRemainder--;
        }
    }
    
    return remainder;
}