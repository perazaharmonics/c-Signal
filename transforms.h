#include <vector>
#include <complex>
#include <cmath>

class Transform {
public:
    // Cooley-Tukey FFT (in-place, divide-and-conquer)
    // Fast Fourier Transform (FFT) is an algorithm that computes the discrete Fourier transform (DFT)
    // and its inverse. This implementation uses the Cooley-Tukey algorithm.
    void fft(std::vector<std::complex<double>>& a) {
        int n = a.size();
        
        // Base case: if input size is 1, nothing to do
        if (n <= 1) return;

        // Divide step: Split the input into even and odd parts
        std::vector<std::complex<double>> a_even(n / 2), a_odd(n / 2);
        for (int i = 0; i < n / 2; i++) {
            a_even[i] = a[i * 2]; // Even indexed elements
            a_odd[i] = a[i * 2 + 1]; // Odd indexed elements
        }

        // Conquer step: Recursively apply FFT to both halves
        fft(a_even);
        fft(a_odd);

        // Define PI constant
        const double PI = 3.14159265358979323846;

        // Combine step: Merge the results of the two halves
        for (int k = 0; k < n / 2; k++) {
            // Compute the complex exponential factor
            std::complex<double> t = std::polar(1.0, -2 * PI * k / n) * a_odd[k];
            
            // Combine even and odd parts with the complex exponential factor
            a[k] = a_even[k] + t;
            a[k + n / 2] = a_even[k] - t;
        }
    }
    
    // Other transformation functions could go here, such as inverse FFT, DCT, etc.
};
