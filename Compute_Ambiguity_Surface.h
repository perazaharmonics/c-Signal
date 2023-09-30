// Include necessary libraries
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// Use the standard namespace
using namespace std;

// Define M_PI if it is not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Define a complex number type
typedef complex<double> Complex;

// Define a class to compute the ambiguity surface
class ComputeAmbiguitySurface {
public:
    // Define a static instance of the class
    static ComputeAmbiguitySurface& instance() {
        static ComputeAmbiguitySurface instance;
        return instance;
    }

    // Define a function to compute the ambiguity surface
    void compute(vector<Complex>& s, vector<double>& t, vector<double>& tau_range, vector<double>& fd_range, vector<vector<Complex>>& ambiguity_surface) {
        // Get the number of tau values, fd values, and time values
        int num_tau = tau_range.size();
        int num_fd = fd_range.size();
        int num_t = t.size();

        // Create vectors to hold the doppler signal and shifted signal
        vector<Complex> doppler_signal(num_t);
        vector<Complex> shifted_signal(num_t);

        // Resize the ambiguity surface vector
        ambiguity_surface.resize(num_fd, vector<Complex>(num_tau, 0));

        // Loop over all tau and fd values
        for (int i = 0; i < num_tau; i++) {
            for (int j = 0; j < num_fd; j++) {
                // Shift the signal by the current fd value
                for (int k = 0; k < num_t; k++) {
                    shifted_signal[k] = s[k] * exp(-1.0 * Complex(0, 1) * 2.0 * M_PI * fd_range[j] * t[k]);
                }

                // Compute the ambiguity function
                Complex sum = 0;
                for (int k = 0; k < num_t; k++) {
                    sum += conj(shifted_signal[k]) * s[k];
                }

                // Store the result in the ambiguity surface vector
                ambiguity_surface[j][i] = sum;
            }
        }
    }

private:
    // Define a private constructor and copy constructor to prevent instantiation
    ComputeAmbiguitySurface() {}
    ComputeAmbiguitySurface(const ComputeAmbiguitySurface&) = delete;
    ComputeAmbiguitySurface& operator=(const ComputeAmbiguitySurface&) = delete;
};