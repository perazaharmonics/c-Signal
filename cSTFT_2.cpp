#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
//#include "imgui.h"


// #include "transforms.h" // Include the header file for the Transform class
// #include <fftw3.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parseval's Theorem
// Reference: https://en.wikipedia.org/wiki/Parseval%27s_theorem
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double myParseval(const std::vector<std::complex<double>>& X) {
    double energy = 0.0;
    cout << "Executing myParseval()" << endl; 
    int N = X.size();
    for (const auto& val : X) {
        energy += std::norm(val); // norm returns magnitude squared
    }
    cout << "Calculated energy: " << energy << endl; // Debug: Check the calculated energy
    return energy / N;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Zero Padding
// Reference: https://en.wikipedia.org/wiki/Zero-padding
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Zero Padding
// Reference: https://en.wikipedia.org/wiki/Zero-padding
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::complex<double>> zeroPadToNextPowerOf2(const std::vector<std::complex<double>>& input) {
    std::cout << "Executing zeroPadToNextPowerOf2()" << std::endl;

    int currentSize = input.size();


    // If the current size is already a power of 2, just return the input.
    if ((currentSize & (currentSize - 1)) == 0) { // Removed extraneous semicolon after if condition
        std::cout << "Size is already a power of 2: " << currentSize << std::endl;  // Debug: Inform that the size is already a power of 2
        return input;
    }

    // Compute the next power of 2
    int nextPowerOf2 = 1;
    while (nextPowerOf2 < currentSize) {
        nextPowerOf2 <<= 1;
    }

    std::cout << "Padding to size: " << nextPowerOf2 << std::endl; // Debug: Inform the new size after padding

    std::vector<std::complex<double>> result = input;
    result.resize(nextPowerOf2, std::complex<double>(0, 0)); // Zero-padding with complex zeros.

    return result;
}


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
// Description: A Kaiser window is a taper formed by using a Bessel function of the first kind, J0, as the shape function.
// Reference: https://en.wikipedia.org/wiki/Kaiser_window
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Window Chooser
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Get user's window choice
std::vector<double> getWindow(int choice, int windowSize) {
    std::vector<double> selectedWindow;
    switch (choice) {
    case 1:
        selectedWindow = hanning(windowSize);
        break;
    case 2:
        selectedWindow = tukey(windowSize);
        break;
    case 3:
        selectedWindow = blackman(windowSize);
        break;
    case 4:
        float beta;
        std::cout << "\nEnter the value of beta: ";
        std::cin >> beta;
        selectedWindow = kaiser(windowSize, beta);
        break;
    case 5:
        float sigma;
        std::cout << "\nEnter the value of sigma: ";
        std::cin >> sigma;
        selectedWindow = gaussian(windowSize, sigma);
        break;
    default:
        std::cout << "Invalid choice. Using Hanning window." << std::endl;
        selectedWindow = hanning(windowSize);
        break;
    }
    return selectedWindow;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convolution
// Description: Convolution is a mathematical operation on two functions that produces a third function expressing how the shape of one is modified by the other.
// The term convolution refers to both the result function and to the process of computing it.
//  Convolution is similar to cross-correlation. In cross-correlation, the filter is flipped before sliding over the signal.
//  In convolution, the filter is not flipped. Convolution is commutative while cross-correlation is not.
// Convolution is used in the mathematics of many fields, such as probability and statistics.
// It is used to calculate the probability of a sum of two independent random variables.
// It is the key to defining the Fourier transform and other Fourier-related transforms.
// In signal processing, it describes the process of filtering a signal by another signal.
//  The convolution of one function (the input) with a second function (the impulse response) produces a third function (the output).
// The Fourier transform of the output function is the product of the Fourier transform of the input function with the Fourier transform of the impulse response function.
// Reference: https://en.wikipedia.org/wiki/Convolution
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<std::complex<double>> convolution(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& y) {
    int n = x.size();
    int m = y.size();
    int N = 1;
    while (N < n + m - 1) {
        N <<= 1;  // Doubling N until it's large enough
    }

    std::cout << "Input sizes: n = " << n << ", m = " << m << ", N = " << N << std::endl; // Debug: Print sizes

    std::vector<std::complex<double>> padded_x(N, 0.0);
    std::vector<std::complex<double>> padded_y(N, 0.0);

    std::copy(x.begin(), x.end(), padded_x.begin());
    std::copy(y.begin(), y.end(), padded_y.begin());

    std::cout << "Padded vector sizes: padded_x = " << padded_x.size() << ", padded_y = " << padded_y.size() << std::endl; // Debug: Print sizes

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


std::vector<std::vector<std::complex<double>>> STFT(const std::vector<std::complex<double>>& signal, const std::vector<double>& window, int overlap) {
    int Nframe = window.size();
    int Nnoise = signal.size();

    // Percentage of overlap
    double poverl = static_cast<double>(overlap) / Nframe;

    // Calculate the number of overlapped frames
    int koverlap = (Nnoise - overlap) / (Nframe - overlap);

    if (koverlap > Nnoise) {
        std::cout << "Error, overlap of less than one sample" << std::endl;
        return {}; // Return an empty vector
    }

    std::vector<std::vector<std::complex<double>>> result(koverlap);

    for (int i = 0; i < koverlap; ++i) {
        int nstart = i * (Nframe - overlap);
        int nend = nstart + Nframe - 1;

        // Create the windowed segment
        std::vector<std::complex<double>> segment(Nframe);
        for (int j = nstart; j <= nend; ++j) {
            segment[j - nstart] = signal[j] * window[j - nstart];
        }

        result[i] = fft(segment);

        // No need for segshifted and sout in this context, as they are used for plotting and accumulation.
        // If you need them for further processing, you can create similar vectors as needed.

        // You can also plot or process the 'segment' vector here if needed.
    }

    return result;
}



/* std::vector<std::complex<double>> generateNoisySignal(int signalSize) {
    std::default_random_engine generator; // Default random engine
    std::normal_distribution<double> AWGN(0.0, 1.0);

    std::vector<std::complex<double>> signal(signalSize);

    // Add AWGN to the signal
    for (auto& val : signal) {
        double realNoise = AWGN(generator);
        double imagNoise = AWGN(generator);
        val = std::complex<double>(realNoise, imagNoise);
    }

    return signal;
}

*/

int main() {
    // Create the input signal (using the complex type as required)
    std::vector<double> signalInRealPart = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    std::vector<std::complex<double>> signalIn(signalInRealPart.size());
    for (size_t i = 0; i < signalInRealPart.size(); ++i) {
        signalIn[i] = std::complex<double>(signalInRealPart[i], 0.0);
    }

    int window_size = 8; // Window size, change as needed
    int overlap = 4;     // Overlap, change as needed
    int hop_size = window_size - overlap;


    std::cout << "\n The input signal sequence: \n";
    for (const double& value : signalInRealPart) {
        std::cout << value << " ";
    }
    std::cout << std::endl;


    int choice;
    std::cout << "Choose a window type:" << std::endl;
    std::cout << "1. Hanning\n2. Tukey\n3. Blackman\n4. Kaiser\n5. Gaussian" << std::endl;
    std::cin >> choice;

    std::vector<double> selected_window = getWindow(choice, window_size);


    // Calculate the Short-Time Fourier Transform (STFT)
    auto stft_result = STFT(signalIn, selected_window, hop_size);


    // Display STFT results
    for (size_t frame = 0; frame < stft_result.size(); ++frame) {
        std::cout << "Frame " << frame << ": ";
        for (size_t bin = 0; bin < stft_result[frame].size(); ++bin) {
            std::cout << std::abs(stft_result[frame][bin]) << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

