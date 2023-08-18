#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
// #include "transforms.h" // Include the header file for the Transform class
// #include <fftw3.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



using namespace std;

std::vector<std::complex<double>> hermitian(const std::vector<std::complex<double>>& vec) {
    std::vector<std::complex<double>> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = std::conj(vec[i]);
    }

    return result;
}

// bit reversal permutation for real numbers
int reverseBits(int num, int log2n) {
    int reversed = 0;
    for (int i = 0; i < log2n; i++) {
        reversed = (reversed << 1) | (num & 1);
        num >>= 1;
    }
    return reversed;
}

// bit reversal permutation for complex numbers
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

// Begin Stride-Permutation FFT Algorithm
std::vector<std::complex<double>> fft_stride(const std::vector<std::complex<double>>& data) {
    int n = data.size();
    if (n <= 1) return data;


    // Re-order the input data by stride permutation (bit-reversal)
    std::vector<std::complex<double>> output_array = stridepermutation(data);

    // Pre-compute the twiddle factors
    std::vector<std::complex<double>> twiddle(n);
    for (int k = 0; k < n; ++k) {
        twiddle[k] = std::polar(1.0, -2 * M_PI * k / n);
    }

    // Perform the FFT
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
    // Return output spectra
    return output_array;
}

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
// Define the convolution function
std::vector<std::complex<double>> convolution(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& y) {
    int n = x.size();
    int m = y.size();
    int N = 1;
    while (N < n + m - 1) N <<= 1; // Find the smallest power of 2 >= n+m-1

    std::vector<std::complex<double>> padded_x(N, 0.0);
    std::vector<std::complex<double>> padded_y(N, 0.0);

    std::copy(x.begin(), x.end(), padded_x.begin());
    std::copy(y.begin(), y.end(), padded_y.begin());

    // Perform the FFT on x and y
    std::vector<std::complex<double>> X = fft_stride(padded_x);
    std::vector<std::complex<double>> Y = fft_stride(padded_y);

    // Perform the element-wise multiplication
    for (int i = 0; i < N; i++) {
        X[i] *= Y[i];
    }

    // Perform the IFFT on the result
    ifft_stride(X);
    X.resize(n + m - 1);  // Truncate or resize to get the correct convolution size

    return X;
}
int main() {
    // Sample input sequence of 8 numbers
    vector<complex<double>> dataIn = { 1, 2, 3, 4, 5, 6, 7, 8 }, Y, Y_conv, Y_T;
    Y_T = hermitian(dataIn);
    Y_conv = convolution(Y_T, dataIn);


    // Call the fft method to obtain output signal from the input data
    Y = fft_stride(dataIn);

    
    cout << "8-Point Fast-Fourier Transform (FFT) using Stride Permutation" << "\n" << endl;

    // Print the results of the FFT
    for (auto& val : Y) { // C++11 range-based for loop
        cout << val <<  endl;
        
        

    }
    cout <<"\n " << "Convolution using FFT: " << "\n" << endl;

    // Print the results of the FFT
    for (auto& val : Y_conv) { // C++11 range-based for loop
        cout << val << endl;



    }
       
      return 0;


}
