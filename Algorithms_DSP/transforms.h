#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <valarray>
#include <fftw3.h>
#include <cstdint>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum WindowType {
    HANNING,
    HAMMING,
    BLACKMAN,
    // Add more window types as needed...
    CUSTOM
};

template <typename T>
std::vector<T> zeroPadToNextPowerOf2(const std::vector<T>& input) {
    int currentSize = input.size();

    // If the current size is already a power of 2, just return the input.
    if ((currentSize & (currentSize - 1)) == 0); {
        return input;
    }


    // Calculate the next power of 2 greater than currentSize.
    int nextPowerOf2 = std::pow(2, std::ceil(std::log2(currentSize)));

    // Create a new vector with the desired size, initializing it with zeros.
    std::vector<T> padded(nextPowerOf2, 0);

    // Copy the input data to the beginning of the padded vector.
    std::copy(input.begin(), input.end(), padded.begin());

    return padded;
}

class Transforms {
    // Firstly, we need to declare the methods that we will be using in the class
   



public: std::vector<std::complex<double>> hermitian(const std::vector<std::complex<double>>& vec) {
    std::vector<std::complex<double>> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = std::conj(vec[i]);
    }

    return result;
}
      // Bit-reversal permutation

      int reverseBits(int num, int log2n) {
          int reversed = 0;
          for (int i = 0; i < log2n; i++) {
              reversed = (reversed << 1) | (num & 1);
              num >>= 1;
          }
          return reversed;
      }
      // Stride permutation

      std::vector<std::complex<double>> stridepermutation(const std::vector<std::complex<double>>& data) {
          int n = data.size();
          int log2n = std::log2(n);
          std::vector<std::complex<double>> result(n);

          for (int i = 0; i < n; ++i) {
              int j = reverseBits(i, log2n);
              result[j] = data[i];
          }

          return result;
      }

public: std::vector<std::complex<double>> convolution(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& y) {
    int n = x.size();
    int m = y.size();
    int N = 1;
    while (N < n + m - 1) N <<= 1; // Find the smallest power of 2 >= n+m-1

    std::vector<std::complex<double>> padded_x(N, 0.0);
    std::vector<std::complex<double>> padded_y(N, 0.0);

    std::copy(x.begin(), x.end(), padded_x.begin());
    std::copy(y.begin(), y.end(), padded_y.begin());

    // Perform the FFT on x and y
    std::vector<std::complex<double>> X = fft(padded_x);
    std::vector<std::complex<double>> Y = fft(padded_y);

    // Perform the element-wise multiplication
    for (int i = 0; i < N; i++) {
        X[i] *= Y[i];
    }

    // Perform the IFFT on the result
    ifft(X);
    X.resize(n + m - 1);  // Truncate or resize to get the correct convolution size

    return X;
}
private:
    std::vector<std::complex<double>> twiddle_factors;

    void twiddle(int N) {
        if (twiddle_factors.size() != N) {
            twiddle_factors.resize(N / 2);
            for (int i = 0; i < N / 2; ++i) {
                twiddle_factors[i] = std::polar(1.0, -2 * M_PI * i / N);
            }
        }
    }

    std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& data) {
        int n = data.size();


        // Base case: if input size is 1, return the input
        if (n <= 1) return data;

        // Divide step: Split the input into even and odd parts
        std::vector<std::complex<double>> data_even(n / 2), data_odd(n / 2);
        for (int i = 0; i < n / 2; i++) {
            data_even[i] = data[i * 2]; // Even indexed elements
            data_odd[i] = data[i * 2 + 1]; // Odd indexed elements
        }

        // Conquer step: Recursively apply FFT to both halves
        std::vector<std::complex<double>> Y_even = fft(data_even);
        std::vector<std::complex<double>> Y_odd = fft(data_odd);

        // Get the twiddle factors
        twiddle(n);
        

        // Combine step: Merge the results of the two halves
        std::vector<std::complex<double>> Y(n);  // The result
        for (int k = 0; k < n / 2; k++) {
            std::complex<double> t = twiddle_factors[k] * Y_odd[k];
            Y[k] = Y_even[k] + t;
            Y[k + n / 2] = Y_even[k] - t;
        }

        return Y;
    }

      // Inverse FFT
      void ifft(std::vector<std::complex<double>>& data) {
          int n = data.size();

          // Step 1: Conjugate the input data 
          for (auto& d : data) {
              d = std::conj(d);
          }

          if (n <= 1) return;

          // Divide
          std::vector<std::complex<double>> data_even(n / 2), data_odd(n / 2);
          for (int i = 0; i < n / 2; ++i) {
              data_even[i] = data[i * 2];
              data_odd[i] = data[i * 2 + 1];
          }

          // Recurse
          ifft(data_even);  // Recursion on even-indexed elements possibly complex valued
          ifft(data_odd);   // Recursion on odd-indexed elements possibly complex valued

          twiddle(n);

          // Combine
          for (int k = 0; k < n / 2; ++k) {
              std::complex<double> t = std::conj(twiddle_factors[k]) * data_odd[k];  // Remember to conjugate the twiddle factor

              data[k] = data_even[k] + t;
              data[k + n / 2] = data_even[k] - t;
          }

          // Step 2: Conjugate the data again and scale
          for (auto& d : data) {
              d = std::conj(d) / static_cast<double>(n);
          }
      }

    
    // Begin Stride-Permutation FFT Algorithm
    std::vector<std::complex<double>> fft_stride(const std::vector<std::complex<double>>& data) {
        int n = data.size();
        if (n <= 1) return data;

        // Re-order the input data by stride permutation (bit-reversal)
        std::vector<std::complex<double>> output_array = stridepermutation(data);

        std::vector<std::complex<double>> twiddle(n);
        for (int k = 0; k < n; ++k) {
            twiddle[k] = std::polar(1.0, -2 * M_PI * k / n);
        }

        for (int s = 1; s <= std::log2(n); ++s) {
            int m = std::pow(2, s);
            int half_m = m / 2;

            for (int k = 0; k < n; k += m) {
                for (int j = 0; j < half_m; ++j) {
                    std::complex<double> t = twiddle[j * n / m] * output_array[k + half_m + j];
                    std::complex<double> u = output_array[k + j];
                    output_array[k + j] = u + t;
                    output_array[k + half_m + j] = u - t;
                }
            }
        }
        return output_array;
    }
    // IFFT By Stride Permutation
    std::vector<std::complex<double>> ifft_stride(const std::vector<std::complex<double>>& data) {
        int n = data.size();
        if (n <= 1) return data;

        // Re-order the input data by stride permutation (bit-reversal)
        std::vector<std::complex<double>> output_array = stridepermutation(data);

        std::vector<std::complex<double>> twiddle(n);
        for (int k = 0; k < n; ++k) {
            twiddle[k] = std::polar(1.0, 2 * M_PI * k / n);
        }

        for (int s = 1; s <= std::log2(n); ++s) {
            int m = std::pow(2, s);
            int half_m = m / 2;

            for (int k = 0; k < n; k += m) {
                for (int j = 0; j < half_m; ++j) {
                    std::complex<double> t = twiddle[j * n / m] * output_array[k + half_m + j];
                    std::complex<double> u = output_array[k + j];
                    output_array[k + j] = u + t;
                    output_array[k + half_m + j] = u - t;
                }
            }
        }

        // Normalization after IFFT
        for (int i = 0; i < n; ++i) {
            output_array[i] /= n;
        }

        return output_array;
    }

public:
    // STFT
    std::vector<std::vector<std::complex<double>>> stft(const std::vector<std::complex<double>>& signal, int window_size, int overlap, const std::vector<double>& custom_window = {}) {
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
    std::vector<std::complex<double>> istft(const std::vector<std::vector<std::complex<double>>>&stft_output, int window_size, int overlap) {
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

};
#endif // TRANSFORMS_H

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



        




    





