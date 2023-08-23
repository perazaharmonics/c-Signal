#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include "matplotlibcpp.h"

// #include "transforms.h" // Include the header file for the Transform class
// #include <fftw3.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



using namespace std;
namespace plt = matplotlibcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hermitian Adjoint to handle complex conjugates
// Reference: https://en.wikipedia.org/wiki/Hermitian_adjoint
/////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::complex<double>> hermitian(const std::vector<std::complex<double>>& vec) {
    std::vector<std::complex<double>> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i)
    {
        result[i] = std::conj(vec[i]);
    }

    return result;
}

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
// Twiddle Factor Function
// Reference: https://en.wikipedia.org/wiki/Twiddle_factor
/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::complex<double>> twiddle_factors;

void twiddle(int N) {
    if (twiddle_factors.size() != N) {
        twiddle_factors.resize(N / 2);
        for (int i = 0; i < N / 2; ++i) {
            twiddle_factors[i] = std::polar(1.0, -2 * M_PI * i / N);
        }
    }
}

// Fast Fourier Transform
// Description: A fast Fourier transform (FFT) is an algorithm that computes the discrete Fourier transform (DFT) of a sequence, using divide and conquer approach.
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Windowing Functions
// 1. Hanning Window
// 2. Tukey Window
// 3. Blackman Window
// 4. Kaiser Window
// 5. Gaussian Window
// Reference: https://en.wikipedia.org/wiki/Window_function
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 1. Hanning Window
// Description: A Hanning window is a type of window function that reduces side lobes ringing in the frequency domain and helps mitigate segmentation noise. 
// It's defined using a raised cosine.
std::vector<double> hanning(int win_size) {
	std::vector<double> window(win_size);

    for (int n = 0; n < win_size; n++) {
		window[n] = 0.5 * (1 - cos((2 * M_PI * n) / (win_size - 1)));
	}

	return window;
}

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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> gaussian(int win_size, double sigma) {
	std::vector<double> window(win_size);

    for (int n = 0; n < win_size; n++) {
		window[n] = exp(-0.5 * pow((n - (win_size - 1) / 2) / (sigma * (win_size - 1) / 2), 2));
	}

	return window;
}



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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// STFT
// Reference:https://en.wikipedia.org/wiki/Short-time_Fourier_transform
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<std::complex<double>>> STFT(std::vector<std::complex<double>>& signal, const std::vector<double>& window, int overlap) {
    int stepSize = window.size() - overlap;
    int numFrames = (signal.size() - overlap) / stepSize;

    std::vector<std::vector<std::complex<double>>> result(numFrames);

    for (int frame = 0; frame < numFrames; frame++) {
        std::vector<std::complex<double>> segment(window.size());

        // Windowing
        for (int i = 0; i < window.size(); i++) {
            segment[i] = signal[frame * stepSize + i] * window[i];
        }

        // Fourier Transform
        result[frame] = fft(segment); // Assuming your FFT can handle complex numbers and you've handled zero-padding
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plot the Spectrogram
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plotSpectrum(const std::vector<std::vector<double>>& energy) {
    // Convert 2D energy spectrum to 1D for plotting
    std::vector<double> energy_1d;
    for (const auto& row : energy) {
        for (const auto& val : row) {
            energy_1d.push_back(val);
        }
    }

    // Plot the 1D energy spectrum
    plt::plot(energy_1d);
    plt::title("Spectrum Energy");
    plt::xlabel("Frequency Bin");
    plt::ylabel("Energy");
    plt::show();
}


int main() {
    // Generate a random signal
    std::default_random_engine generator; // default random engine
    std::normal_distribution<double> distribution(0.0, 1.0); // Gaussian distribution with mean = 0 and standard deviation = 1

    int fs = 44100; // Sampling frequency
    int Noise_size = fs / 2 + 1; // Size of the noise signal
    int window_size = 256; // Size of the window
    int hop_size = window_size / 2; // Size of the hop

    std::vector<std::complex<double>> X_noise(Noise_size);
    for (int i = 0; i < Noise_size; i++) {
        X_noise[i] = distribution(generator);
    }
    std::vector<std::complex<double>> Y_noisy;
    Y_noisy = convolution(X_noise, hermitian(X_noise));



    for (const auto& val : Y_noisy) {
        cout << val.real() << " + " << val.imag() << "i" << endl;
    }
    cout << "\nChoose a Window for STFT:" << endl;
    cout << "1. Hanning Window" << endl;
    cout << "2. Tukey Window" << endl;
    cout << "3. Blackman Window" << endl;
    cout << "4. Kaiser Window" << endl;
    cout << "5. Gaussian Window" << endl;

    int choice;
    cout << "Enter your choice (1-5): ";
    cin >> choice;

    std::vector<double> selected_window;

    int window_size = 8; // Sample value, modify accordingly
    switch (choice) {
    case 1:
        selected_window = hanning(window_size);
        break;
    case 2:
        selected_window = tukey(window_size);
        break;
    case 3:
        selected_window = blackman(window_size);
        break;
    case 4:
        double beta; // Kaiser window needs a beta value
        cout << "Enter beta value for Kaiser Window: ";
        cin >> beta;
        selected_window = kaiser(window_size, beta);
        break;
    case 5:
        double sigma; // Gaussian window needs a sigma value
        cout << "Enter sigma value for Gaussian Window: ";
        cin >> sigma;
        selected_window = gaussian(window_size, sigma);
        break;
    default:
        cout << "Invalid choice. Exiting." << endl;
        return 0;
    }
    std::vector<double> win = selected_window;
    std::vector<std::vector<std::complex<double>>> result = STFT(Y_noisy, win, 0);
    // Convert complex STFT result to energy
    std::vector<std::vector<double>> energy(result.size());
    for (int i = 0; i < result.size(); i++) {
        energy[i].resize(result[i].size());
        for (int j = 0; j < result[i].size(); j++) {
            energy[i][j] = std::norm(result[i][j]);
        }
    }

    plotSpectrum(energy);

    return 0;



}

