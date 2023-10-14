/* Here is the explanation for the code above:
1. The first part of the code is to initialize the matrix that will be used to store the values of the Ambiguity Function (AF) for each combination of tau and f.
2. The next part of the code is to calculate the AF for each combination of tau and f. This part is parallelized using OpenMP.
3. The final part of the code is to print the result for each combination of tau and f. This part is not parallelized.
4. The code then stops the timer and prints the duration of the calculation. */


#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

using namespace std;
using namespace std::chrono;

typedef vector<double> Signal;
typedef vector<vector<double>> AFMatrix;

AFMatrix computeAF(const Signal& x, const vector<double>& tau_values, const vector<double>& f_values) {
    int tau_length = tau_values.size();
    int f_length = f_values.size();
    AFMatrix AF(tau_length, vector<double>(f_length, 0));

    // Start the timer
    auto start = system_clock::now();
    time_t start_time = system_clock::to_time_t(start);
    char str[26];
    ctime_s(str, sizeof str, &start_time);
    cout << "Starting the timer at: " << str; // Printing the start time

#pragma omp parallel for collapse(2)
    for (int m = 0; m < tau_length; m++) {
        for (int n = 0; n < f_length; n++) {
            double tau = tau_values[m];
            double f = f_values[n];
            double result = 0;

#pragma omp simd reduction(+:result)
            for (size_t t = 0; t < x.size(); t++) {
                double time_val = t + tau;
                if (time_val >= 0 && time_val < x.size()) {
                    result += x[t] * x[static_cast<int>(time_val)];
                }
            }

            AF[m][n] = result;

            if (m % 10000 == 0 && n % 10000 == 0) {
                cout << "For Time Delay " << tau << ", Doppler Shift " << f << ", AF = " << result << endl;
            }
        }
    }

    // Stop the timer and calculate the duration
    auto stop = system_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    // Print the duration
    cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;

    return AF;
}
