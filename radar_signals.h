#include <vector>
#include <string>
#include <cmath>    // Include cmath for mathematical functions
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <stdio.h>
// Define a class called radar_signals

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class radar_signals {
public:
    static radar_signals& instance() {
        static radar_signals instance;
        return instance;
    }

    // Define a function called quadratic_chirp that takes in a vector of doubles called t, and three doubles called f0, f1, and T
    std::vector<double> quadratic_chirp(const std::vector<double>& t, double f0, double f1, double T) const {
        // Create an empty vector called chirp
        std::vector<double> chirp;
        // Reserve space in the chirp vector for the size of the t vector
        chirp.reserve(t.size());
        // Calculate the beta value
        const double beta = (f1 - f0) / (T * T);
        // Loop through each value in the t vector
        for (const double& t_ : t)
        {
            // Calculate the value of the chirp waveform at the current time
            chirp.push_back(std::cos(2 * M_PI * (f0 * t_ + beta * (t_ * t_ * t_) / 3)));
        }
        // Return the chirp waveform
        return chirp;
    }

    // Define a function called exponential_chirp that takes in a vector of doubles called t, and three doubles called f0, f1, and T
    std::vector<double> exponential_chirp(const std::vector<double>& t, const double f0, const double f1, const double T) const {
        // Create an empty vector called chirp
        std::vector<double> chirp;
        // Reserve space in the chirp vector for the size of the t vector
        chirp.reserve(t.size());
        // Calculate the beta value
        const double beta = (f1 - f0) / (T * T);
        // Loop through each value in the t vector
        for (const double& t_ : t)
        {
            // Calculate the value of the chirp waveform at the current time
            chirp.push_back(std::cos(2 * M_PI * (f0 * t_ + beta * (t_ * t_ * t_) / 3)));
        }
        // Return the chirp waveform
        return chirp;
    }

    // Define a function called bpsk_modulate that takes in a vector of doubles called data, a vector of doubles called t, and two doubles called carrier_frequency and modulation_index
    std::vector<double> bpsk_modulate(const std::vector<double>& data, const std::vector<double>& t, double carrier_frequency, double modulation_index) const {
        std::vector<double> bpsk_waveform;
        bpsk_waveform.reserve(t.size());

        // Correct the implementation according to BPSK modulation technique
        for (size_t i = 0; i < data.size(); ++i)
        {
            double carrier = std::cos(2 * 3.14159265358979323846 * carrier_frequency * t[i]);
            double phase_shift = (data[i] == 0) ? 3.14159265358979323846 : 0;  // Phase shift depending on the data bit
            bpsk_waveform.push_back(carrier * std::cos(phase_shift));
        }
        return bpsk_waveform;
    }

private:
    radar_signals() {}
};
