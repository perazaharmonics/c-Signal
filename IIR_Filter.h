#pragma once    // Include guard

#include "Iir.h"
#include <vector>
#include "Stk.h"
#include <iostream>

using namespace stk;

class IIRFilter {
public:
    IIRFilter(const std::vector<StkFloat>& numerator, const std::vector<StkFloat>& denominator);
    void applyFilter(StkFrames& frames);
    void printOutput(const StkFrames& frames) const;

private:
    Iir filter_;
    StkFrames output_;
};
