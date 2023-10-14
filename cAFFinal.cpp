#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <random>
#include <omp.h>
#include <immintrin.h> // for AVX intrinsics
#include "computeAF.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI

using namespace std;

int main() {
    // Constants
    const int fs = 65536;  // sampling frequency

    // Generate a random pseudo-noise signal
    vector<double> pnSignal(fs, 0); // initialize with 0
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);

#pragma omp parallel for
    for (int i = 0; i < fs; ++i) {
        pnSignal[i] = dis(gen) * 2 - 1; // converting 0,1 to -1,1 for BPSK
    }

    vector<double> tau_values;
    vector<double> f_values;

    double tau_start = 0;     // start value of time delay
    double tau_end = 1;       // end value of time delay
    double tau_step = 1;   // step size of time delay

    double f_start = -0.5;    // start value of normalized frequency
    double f_end = 0.5;       // end value of normalized frequency
    double f_step = 1.0 / fs; // step size of normalized frequency

    // fill tau_values with values from tau_start to tau_end with a step size of tau_step
    for (double tau = tau_start; tau <= tau_end; tau += tau_step) {
        tau_values.push_back(tau);
    }

    // fill f_values with values from f_start to f_end with a step size of f_step
    for (double f = f_start; f <= f_end; f += f_step) {
        f_values.push_back(f);
    }

    AFMatrix af = computeAF(pnSignal, tau_values, f_values);

    // Print or process the AF matrix as required.
    return 0;
}
