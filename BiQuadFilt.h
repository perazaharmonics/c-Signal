#pragma once    // Include guard

#include "BiQuad.h"
#include "Noise.h"
#include <iostream>

using namespace stk;

class BiQuadFilt {
public:
    BiQuadFilt();
    void setResonance(StkFloat frequency, StkFloat radius, bool normalize = false);
    void applyFilter(Noise& noiseSource, StkFrames& frames);
    void printOutput(const StkFrames& frames) const;

private:
    BiQuad biquad_;
    StkFrames output_;
};
