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

