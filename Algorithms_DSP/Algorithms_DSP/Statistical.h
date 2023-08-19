#ifndef STATISTICAL_H
#define STATISTICAL_H
#endif // STATISTICAL_H

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>  // Add this
#include <iostream>
#include <valarray>
#include <numeric>    // For std::accumulate
#include <fftw3.h>
#include <cstdint>

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
    std::vector<T> fNLM_denoise_impl(const std::vector<T>& y, int search_window, int patch_size, double h) {
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
    //////////////////////////////////////////////////////////
    // NLM denoising for real-valued signals
    ??
    public:
        std::vector<double> fNLM_denoise(const std::vector<double>& y, int search_window, int patch_size, double h) {
            return fNLM_denoise_impl(y, search_window, patch_size, h);
        }
        // NLM denoising for real-valued signals
        std::vector<std::complex<double>> fNLM_denoise(const std::vector<std::complex<double>>& y, int search_window, int patch_size, double h) {
            return fNLM_denoise_impl(y, search_window, patch_size, h);
        }

        // Weiner Filter for real-valued signals
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
