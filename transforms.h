#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include<valarray>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class transforms {

public: std::vector<std::complex<double>> hermitian(const std::vector<std::complex<double>>& vec) {
    std::vector<std::complex<double>> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = std::conj(vec[i]);
    }

    return result;
}
public:
    // Cooley-Tukey FFT (in-place, divide-and-conquer)
    void fft(std::vector<std::complex<double>>& data) {
        int n = data.size();

        // Base case: if input size is 1, nothing to do
        if (n <= 1) return;

        // Divide step: Split the input into even and odd parts
        std::vector<std::complex<double>> data_even(n / 2), data_odd(n / 2);
        for (int i = 0; i < n / 2; i++) {
            data_even[i] = data[i * 2]; // Even indexed elements
            data_odd[i] = data[i * 2 + 1]; // Odd indexed elements
        }

        // Conquer step: Recursively apply FFT to both halves
        fft(data_even);
        fft(data_odd);


        // Combine step: Merge the results of the two halves
        for (int k = 0; k < n / 2; k++) {

            // Compute the complex exponential factor
            std::complex<double> t = std::polar(1.0, -2 * M_PI * k / n) * data_odd[k];

            // Combine even and odd parts with the complex exponential factor
            data[k] = data_even[k] + t;
            data[k + n / 2] = data_even[k] - t;
        }

    }

    // Inverse FFT
    void ifft(std::vector<std::complex<double>>& data) {

        int n = data.size();
        if (n <= 1) return;

        // Divide
        std::vector<std::complex<double>> data_even(n / 2), data_odd(n / 2);
        for (int i = 0; i < n / 2; ++i) {
            data_even[i] = data[i * 2];
            data_odd[i] = data[i * 2 + 1];
        }

        // Recurse
        ifft(data_even);
        ifft(data_odd);


        // Combine
        for (int k = 0; k < n / 2; ++k) {
            std::complex<double> t = std::polar(1.0, -2 * M_PI * k / n) * data_odd[k];
            data[k] = data_even[k] + t;
            data[k + n / 2] = data_even[k] - t;

        }

    }

public:
    // STFT
    std::vector<std::vector<std::complex<double>>> stft(const std::vector<std::complex<double>>& signal, int window_size, int overlap) {
        // Step 1: Split the signal into overlapping windows
        int hop_size = window_size - overlap;
        int num_windows = (signal.size() - overlap) / hop_size;

        std::vector<std::vector<std::complex<double>>> result; // 2D vector to hold the STFT result

        // Hanning window function
        auto hanning = [window_size](int n) {
            return 0.5 * (1 - cos((2 * M_PI * n) / (window_size - 1)));
            };

        for (int w = 0; w < num_windows; w++) {
            // Extract a window from the signal
            std::vector<std::complex<double>> window(window_size);

            for (int n = 0; n < window_size; n++) {
                // Step 2: Apply the window function
                window[n] = signal[w * hop_size + n] * hanning(n);
            }

            // Step 3: Compute the FFT of the window
            fft(window);

            result.push_back(window);
        }

        return result;
    }

    public:
        std::vector<std::complex<double>> istft(const std::vector<std::vector<std::complex<double>>>& stft_output, int window_size, int overlap) {
            int hop_size = window_size - overlap;
            int signal_length = stft_output.size() * hop_size + overlap;

            // Hanning window function
            auto hanning = [window_size](int n) {
                return 0.5 * (1 - cos((2 * M_PI * n) / (window_size - 1)));
                };

            std::vector<std::complex<double>> reconstructed_signal(signal_length, 0);

            for (size_t w = 0; w < stft_output.size(); w++) {
                // Take the inverse FFT of each segment
                std::vector<std::complex<double>> segment = stft_output[w];
                ifft(segment);

                // Multiply by the window function and overlap-add into the result
                for (int n = 0; n < window_size; n++) {
                    reconstructed_signal[w * hop_size + n] += segment[n] * hanning(n);
                }
            }

            return reconstructed_signal;
        }



    // Placeholder for Discrete Cosine Transform
    //void dct(std::vector<double>& a) {
        // Implementation here
    //}

    // Placeholder for Inverse Discrete Cosine Transform
    //void idct(std::vector<double>& a) {
        // Implementation here
    //}

    // ... you can add more transforms as required
};


//////////////////////////////////////////////////////////////
//Example on how to instantiate the object
//
//
//
//////////////////////////////////////////////////////////////
//int main() {
//    transforms dsp;
//                      Test your FFT or other transforms here
//    return 0;
//}
/////////////////////////////////////////////////////////////////////





