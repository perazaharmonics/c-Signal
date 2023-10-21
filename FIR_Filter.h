#pragma once
#ifndef FIR3_FILTER_H
#define FIR3_FILTER_H

#include <vector>

// Type for signal and coefficient data
typedef double word;

// FIR3 Filter Variables Structure
struct fir3Vars {
    std::vector<word> Aout;  // Output array
    std::vector<word> Ain;   // Input array
    word b0;  // Filter coefficient
    word b1;  // Filter coefficient
    word b2;  // Filter coefficient
    word s1;  // Filter state
    word s2;  // Filter state
};

// FIR3 Filter Function Declaration
void fir3(fir3Vars* a);

#endif // FIR3_FILTER_H