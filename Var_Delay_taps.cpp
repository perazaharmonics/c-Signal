#include <iostream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>  // For std::setprecision and std::fixed
#include "FileWvIn.h"
#include "FileWvOut.h"
#include "FileLoop.h"
#include "Stk.h"
#include "MVarDelayLine.h"

using namespace stk;

int main()
{
    // Set the global sample rate before creating class instances.
    Stk::setSampleRate(44100.0);

    FileLoop input;
    FileWvOut output;

    // Load the sine wave file.
    try {
        input.openFile("Test_Signal.wav");
    }
    catch (StkError&) {
        std::cerr << "Could not load the input file." << std::endl;
        exit(1);
    }

    // Open a 16-bit per-channel, one-channel WAV formatted output file
    try {
        output.openFile("BB_Anal2.wav", 1, FileWrite::FILE_WAV, Stk::STK_SINT16);
    }
    catch (StkError&) {
        std::cerr << "Could not open the output file." << std::endl;
        exit(1);
    }

    // Set precision for console output
    std::cout << std::fixed << std::setprecision(5);

    setdelay(100);  // Example: setting a delay of 100 samples

    // Run the oscillator for 40000 samples, writing to the output file
    for (int i = 0; i < 40000; i++) {
        double sample = input.tick();
        double delayedSample = delayline(sample);
        output.tick(delayedSample);

        // Pretty print the original and delayed samples to the console
        std::cout << "Sample " << std::setw(5) << i
            << ": Original = " << std::setw(10) << sample
            << ", Delayed = " << std::setw(10) << delayedSample << std::endl;
    }

    return 0;
}
