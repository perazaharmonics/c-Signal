#include "BiQuadFilt.h"

BiQuadFilt::BiQuadFilt()
    : output_(20, 1)  // Adjust the size and channels as needed
{
    output_[0] = 1.0;  // Initial value, adjust as needed
}

void BiQuadFilt::setResonance(StkFloat frequency, StkFloat radius, bool normalize) {
    biquad_.setResonance(frequency, radius, normalize);
}

void BiQuadFilt::applyFilter(Noise& noiseSource, StkFrames& frames) {
    biquad_.tick(noiseSource.tick(frames));
}

void BiQuadFilt::printOutput(const StkFrames& frames) const {
    for (unsigned int i = 0; i < frames.size(); i++) {
        std::cout << "i = " << i << " : output = " << frames[i] << std::endl;
    }
}
