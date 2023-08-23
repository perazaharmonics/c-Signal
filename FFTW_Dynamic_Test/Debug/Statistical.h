#ifndef STATISTICAL_H
#define STATISTICAL_H
#endif // STATISTICAL_H

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <valarray>
#include <numeric>    // For std::accumulate
// #include <fftw3.h>
#include <cstdint>
#include "transforms.h"
#include "Audiohandler.h"
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Statistical
    {
public:
private:
    template <typename T>
    T squaredDifference(const T& a, const T& b) {
        T diff = a - b;
        return diff * diff; // Generic implementation
    }

    double squaredDifference(const std::complex<double>& a, const std::complex<double>& b) {
        std::complex<double> diff = a - b;
        return std::norm(diff); // Specialization for complex values
    }

    template <typename T>
    std::vector<T> NLM_denoise(const std::vector<T>& y, int search_window, int patch_size, double h) {
        int N = y.size();
        std::vector<T> y_denoised(N, 0.0);

        for (int i = 0; i < N; ++i) {
            int start_search = std::max(0, i - search_window);
            int end_search = std::min(N - 1, i + search_window);
            std::vector<double> weights(end_search - start_search + 1, 0.0);

            for (int j = start_search, w_idx = 0; j <= end_search; ++j, ++w_idx) {
                int start_patch_i = std::max(0, i - patch_size);
                int end_patch_i = std::min(N - 1, i + patch_size);

                int start_patch_j = std::max(0, j - patch_size);
                int end_patch_j = std::min(N - 1, j + patch_size);

                int d = std::min(end_patch_i - start_patch_i + 1, end_patch_j - start_patch_j + 1);

                double squared_diff_sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    squared_diff_sum += squaredDifference(y[start_patch_i + k], y[start_patch_j + k]);
                }

                weights[w_idx] = std::exp(-squared_diff_sum / (h * h));
            }

            double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
            for (double& w : weights) {
                w /= weight_sum;
            }

            for (int j = start_search, w_idx = 0; j <= end_search; ++j, ++w_idx) {
                y_denoised[i] += y[j] * weights[w_idx];
            }
        }

        return y_denoised;
    }
    /////////////////////////////////////////////////////////////////////////////
    // NLM denoising for real-valued signals
    ////////////////////////////////////////////////////////////////////////////////
    public:
        std::vector<double> fNLM_denoise(const std::vector<double>& y, int search_window, int patch_size, double h) {
            return NLM_denoise(y, search_window, patch_size, h);
        }
        // NLM denoising for real-valued signals
        std::vector<std::complex<double>> fNLM_denoise(const std::vector<std::complex<double>>& y, int search_window, int patch_size, double h) {
            return NLM_denoise(y, search_window, patch_size, h);
        }
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // Windowing functions
     ///////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Bessel Function of the first kind
        // Reference: https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/bessel/bessel.html
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Approximation to the Modified Bessel function of the first kind, order zero.
        double bessel_i0(double x) {
            double sum = 1.0;
            double term = 1.0;
            double squaredXDividedBy4 = x * x / 4.0;
            double i = 1.0;

            while (term > 1e-6 * sum) {
                term *= squaredXDividedBy4 / (i * i);
                sum += term;
                i++;
            }

            return sum;
        }
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 1. Hanning Window
        // Description: A Hanning window is a type of window function that reduces side lobes ringing in the frequency domain and helps mitigate segmentation noise. 
        // It's defined using a raised cosine
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        std::vector<double> hanning(int win_size) {
            std::vector<double> window(win_size);

            for (int n = 0; n < win_size; n++) {
                window[n] = 0.5 * (1 - cos((2 * M_PI * n) / (win_size - 1)));
            }

            return window;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 2. Tukey Window
        //Description: A Tukey window is a taper formed by joining a cosine taper to the constant section.
        // It can be seen as a blend between a rectangular window and a Hanning window, with a top that allows for flat-fading.
        // This occurs at the expense of tolerable ripple sidelobes.
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<double> tukey(int win_size) {
            std::vector<double> window(win_size);

            for (int n = 0; n < win_size; n++) {
                if (n < win_size / 2) {
                    window[n] = 0.5 * (1 + cos((2 * M_PI * n) / (win_size - 1)));
                }
                else if (n > win_size / 2) {
                    window[n] = 0.5 * (1 + cos((2 * M_PI * (n - win_size)) / (win_size - 1)));
                }
                else {
                    window[n] = 1;
                }
            }

            return window;


        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 3. Blackman Window
        // Description: A Blackman window provides greater reduction in sidelobe levels at the expense of
        //a wider mainlobe compared to the Hanning window. It's defined using a polynomial function
        //  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        std::vector<double> blackman(int win_size) {
            std::vector<double> window(win_size);

            for (int n = 0; n < win_size; n++) {
                window[n] = 0.42 - 0.5 * cos((2 * M_PI * n) / (win_size - 1)) + 0.08 * cos((4 * M_PI * n) / (win_size - 1));
            }

            return window;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 4. Kaiser Window
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<double> kaiser(int win_size, double beta) {
            std::vector<double> window(win_size);

            for (int n = 0; n < win_size; n++) {
                window[n] = bessel_i0(beta * sqrt(1 - pow((2 * n) / (win_size - 1) - 1, 2))) / bessel_i0(beta);
            }

            return window;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 5. Gaussian Window
        // Description: A Gaussian window is defined by a Gaussian function. It provides a better trade-off between mainlobe width and sidelobe levels.
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<double> gaussian(int win_size, double sigma) {
            std::vector<double> window(win_size);

            for (int n = 0; n < win_size; n++) {
                window[n] = exp(-0.5 * pow((n - (win_size - 1) / 2) / (sigma * (win_size - 1) / 2), 2));
            }

            return window;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Statistical Filters
        // 1. Mean Filter
        // 2. Median Filter
        // 3. Mode Filter
        // 4. Midpoint Filter
        // 5. Alpha-trimmed Mean Filter
        // 6. Weiner Filter
        // 7. Non-Local Means Filter
        //////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Wiener filter for real-valued signals
    //////////////////////////////////////////////////////////////////////////////////////////////
    
        
    std::vector<std::complex<double>> weiner(const std::vector<std::complex<double>>& Y, int n_fft, double Fs, double alpha, double beta) {
        // Compute the magnitude squared (power) of Y
        std::vector<double> Y_power(Y.size());

        // Magnitude squared of complex number is the norm
        std::transform(Y.begin(), Y.end(), Y_power.begin(), [](const std::complex<double>& y) { return std::norm(y); });

        // Compute the noise power (beta)'s spectral estimate
        double N_val = beta * std::accumulate(Y_power.begin(), Y_power.end(), 0.0) / Y_power.size();

        // Compute the signal power spectral estimate
        std::vector<double> S_temp(Y_power.size());
        std::transform(Y_power.begin(), Y_power.end(), S_temp.begin(), [N_val](double y_power) { return std::max(y_power - N_val, 0.0); });

        // Add regularization parameter
        std::vector<double> S_reg(Y_power.size());
        std::transform(S_temp.begin(), S_temp.end(), S_reg.begin(), [N_val, alpha](double s) { return s + alpha * N_val; });

        // Compute Wiener filter coefficients
        std::vector<double> W(Y_power.size());
        std::transform(S_reg.begin(), S_reg.end(), Y_power.begin(), W.begin(), [](double s_reg, double y_power) { return s_reg / y_power; });

        // Apply Wiener filter to input power spectrum
        std::vector<std::complex<double>> Y_filt(Y.size());
        std::transform(W.begin(), W.end(), Y.begin(), Y_filt.begin(), [](double w, const std::complex<double>& y) { return w * y; });

        return Y_filt;
    }

};
