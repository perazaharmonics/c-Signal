#ifndef _SignalOps_
#define _SignalOps_

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <valarray>
#include <cstdint>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// Define alias templates for complex vectors and matrices.
template <typename T>
using complexV_t = vector<complex<T>>;
template <typename R>
using realV_t = vector<R>; // Template for real vectors
template <typename M>
using matrix_t = vector<vector<complex<M>>>; // Template for matrix vectors.

class SignalOps
{
    // Your SignalOps class implementation
public:
    enum class WindowType
    {
        Hanning,
        Hamming,
        BlackmanHarris,
        ExactBlackman,
        Blackman,
        FlatTop,
        FourTermBHarris,
        SevenTermBHarris,
        LowSideLobe
    };

    SignalOps(const complexV_t<double> &s);

    SignalOps(const complexV_t<double> &s, const WindowType &w);

    // Parametrized constructor
    SignalOps(const complexV_t<double> &s, WindowType w, int windowSize, int overlap);

    // Destructor
    ~SignalOps();

    // Public Methods
    complexV_t<double> Convolution(const complexV_t<double> &s, const complexV_t<double> &h, const WindowType &w);
    complexV_t<double> Convolution(const complexV_t<double> &s, const complexV_t<double> &h);
    complexV_t<double> Hermitian(const complexV_t<double> &vec);
    complexV_t<double> FFT(const complexV_t<double> &s, const WindowType &w);
    complexV_t<double> FFT(const complexV_t<double> &s);
    complexV_t<double> IFFT(const complexV_t<double> &s);
    matrix_t<double> STFT(const complexV_t<double> &s, const WindowType &w, const int wSiz, const int overlap);
    complexV_t<double> ISTFT(const matrix_t<double> &sMat, const WindowType &w, const int wSiz, const int overlap);
    complexV_t<double> OLAProcessor(const complexV_t<double> &s, const complexV_t<double> &h, const WindowType &w, const int wSiz, const int overlap);

protected:
    realV_t<double> GenerateWindow(const WindowType &w, const int N);
    complexV_t<double> ZeroPad(const complexV_t<double> &s);

private:
    void TwiddleFactor(int N);   // Precompute the Twiddle for the FFT.
    complexV_t<double> signal;   // The signal to process.
    realV_t<double> window;      // The window to apply to the signal.
    WindowType windowType;       // The Window function.
    const int windowSize;        // The size of the window.
    const int overlap;           // The overlap factor.
    complexV_t<double> twiddles; // Precomputed twiddle factors.
};
// ------------------------------------ //
// Constructors and Destructors
// ------------------------------------ //
SignalOps::SignalOps(const complexV_t<double> &s)
    : signal{s}, windowType{WindowType::Hanning}, windowSize{0}, overlap{0}
{
    // No window generation is needed for this constructor
}

SignalOps::SignalOps(const complexV_t<double> &s, const WindowType &w)
    : signal{s}, windowType{w}, windowSize{0}, overlap{0}
{
    // Generate the appropriate window for the given window size
    window = GenerateWindow(windowType, signal.size());
}

SignalOps::SignalOps(const complexV_t<double> &s, WindowType w, int winSize, int ovlp)
    : signal{s}, windowType{w}, windowSize{winSize}, overlap{ovlp}
{
    // Generate the appropriate window for the given window size
    window = GenerateWindow(windowType, windowSize);

    // Zero-pad the signal to the next power of 2
    signal = ZeroPad(signal);
}

// Object destructor uses the default vector garbage collection.
SignalOps::~SignalOps(void) {}

// ======================== Utility Methods ================================= //
// Twiddle Factor Precomputation method.
// ==========================================================================
void SignalOps::TwiddleFactor(int N)
{
    if (twiddles.size() != N / 2)                        // Did we precompute N/2 twiddle before?
    {                                                    // No, so we..
        twiddles.resize(N / 2);                          // Resize the twiddle factor vector.
        for (int i = 0; i < N / 2; ++i)                  //  loop for the N/2 points and
            twiddles[i] = polar(1.0, -2 * M_PI * i / N); //  compute the twiddle factors.
    }
}

complexV_t<double> SignalOps::ZeroPad(const complexV_t<double> &s)
{
    int N = s.size();                       // Get the length of the input signal.
    if (N && ((N & (N - 1)) == 0))          // Is N already a power of 2?
        return s;                           // Yes, so just return the signal.
    int powerOfTwo = pow(2, ceil(log2(N))); // Get next power of two.
    // ---------------------------------- //
    // Generate a vector with the next power of two.
    // ---------------------------------- //
    complexV_t<double> sPad(powerOfTwo, complex<double>(0, 0));
    // ---------------------------------- //
    // Copy the input data to the vector of quadratic length.
    // ---------------------------------- //
    copy(s.begin(), s.end(), sPad.begin());
    return sPad; // Return quadratic length vector.
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hermitian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Method to compute the Hermitian transform of a vector. The hermitian
// transform is important in multidimensional signal processing and FFTs as it
// reduces the overhead of inverting matrices by just taking a complex conjugate.
// This is legal because the Vandermonde matrix (DFT matrix) is symmetric.
// Reference: https://www.seas.upenn.edu/~ese2240/wiki/Lecture%20Notes/sip_PCA.pdf
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
complexV_t<double> SignalOps::Hermitian(const complexV_t<double> &vec)
{
    complexV_t<double> Ah(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        Ah[i] = conj(vec[i]);
    return Ah;
}
// ======================== Computational Methods =========================== //
// Spectral methods utilized to obtain Time-Frequency Transforms.
// ========================================================================== //
// A function for generating Spectral Windows for windowed Fast Fourier
// Transforms. The need of windowing the FFT comes from minimizing the effects
// of the Gibbs phenomenon or "ringing" that arises from truncating the infinite
// series produced by the Fourier Transform.
// Reference: https://www.mathworks.com/help/signal/ug/windows.html
//
// ======================== Computational Methods =========================== //
// Spectral methods utilized to obtain Time-Frequency Transforms.
// ==========================================================================
realV_t<double> SignalOps::GenerateWindow(const WindowType &w, int N)
{
    realV_t<double> window(N); // Create a window vector for N points.
    switch (w)
    {
    case WindowType::Hanning:
        for (int n = 0; n < N; ++n)
            window[n] = 0.5 * (1 - cos(2 * M_PI * n / N));
        break;
    case WindowType::Hamming:
        for (int n = 0; n < N; ++n)
            window[n] = 0.54 - 0.46 * cos(2 * M_PI * n / N);
        break;
    case WindowType::BlackmanHarris:
        for (int n = 0; n < N; ++n)
            window[n] = 0.35875 - 0.48829 * cos(2 * M_PI * n / (N - 1)) + 0.14128 * cos(4 * M_PI * n / (N - 1)) - 0.01168 * cos(6 * M_PI * n / (N - 1));
        break;
    case WindowType::ExactBlackman:
        for (int n = 0; n < N; ++n)
            window[n] = 0.42659071 - 0.49656062 * cos(2 * M_PI * n / (N - 1)) + 0.07684867 * cos(4 * M_PI * n / (N - 1));
        break;
    case WindowType::Blackman:
        for (int n = 0; n < N; ++n)
            window[n] = 0.42 - 0.5 * cos(2 * M_PI * n / (N - 1)) + 0.08 * cos(4 * M_PI * n / (N - 1));
        break;
    case WindowType::FlatTop:
        for (int n = 0; n < N; ++n)
            window[n] = 0.21557895 - 0.41663158 * cos(2 * M_PI * n / (N - 1)) + 0.277263158 * cos(4 * M_PI * n / (N - 1)) - 0.083578947 * cos(6 * M_PI * n / (N - 1)) + 0.006947368 * cos(8 * M_PI * n / (N - 1));
        break;
    case WindowType::FourTermBHarris:
        for (int n = 0; n < N; ++n)
            window[n] = 0.35875 - 0.48829 * cos(2 * M_PI * n / (N - 1)) + 0.14128 * cos(4 * M_PI * n / (N - 1)) - 0.01168 * cos(6 * M_PI * n / (N - 1));
        break;
    case WindowType::SevenTermBHarris:
        for (int n = 0; n < N; ++n)
            window[n] = 0.27105140069342 - 0.43329793923448 * cos(2 * M_PI * n / (N - 1)) + 0.21812299954311 * cos(4 * M_PI * n / (N - 1)) - 0.06592544638803 * cos(6 * M_PI * n / (N - 1)) + 0.01081276974919 * cos(8 * M_PI * n / (N - 1));
        break;
    case WindowType::LowSideLobe:
        for (int n = 0; n < N; ++n)
            window[n] = 0.35875 - 0.48829 * cos(2 * M_PI * n / (N - 1)) + 0.14128 * cos(4 * M_PI * n / (N - 1)) - 0.01168 * cos(6 * M_PI * n / (N - 1));
        break;
    default:
        for (int n = 0; n < N; ++n)
            window[n] = 0.5 * (1 - cos(2 * M_PI * n / N));
        break;
    }
    return window;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// So to get the FFT of a signal x(n) of length N we have to divide and conquer.
// We do this by using the Cooley-Tukey Algorithm, described here as:
// 1. Divide the signal into even and odd samples, so we have:
//      for (i.begin(); i.end() )
//        evenT[i]=x[i*2]
//        oddT[i]=x[i*2+1]
// 2. Conquer the signal; we recursively apply the FFT on both halves:
//      X_even(k) = FFT(evenT)
//      X_odd(k) = FFT(oddT)
// 3. Now we precompute the twiddle factors, this saves A LOT of time:
//      TwiddleFactor(k) = exp( -j * (2*pi/n)*k)
// 4. Having the Twiddle Factors we compute the FFT butterfly to obtain the full
//    frequency spectrum - the amplitude and phase of sines and cosines that
//    composit it.
//      for (k.begin; k.end/2)
//        t=TwiddleFactor(k)*X_odd(k)
//        X(k)=X_even(k)+t
//        X(k+N/2)=X_even(k)-t
// 5. Return the spectrum of the signal
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
complexV_t<double> SignalOps::FFT( // -------------- FFT --------------- //
    const complexV_t<double> &s)   // Our input time-domain signal.
{                                  // ---------------------------------- //
    int N = s.size();              // The length of the input signal.
    // ---------------------------------- //
    // 1. Base Case: if the input size is 1, return the signal.
    // ---------------------------------- //
    if (N <= 1)   // Is the signal smaller than one sample?
        return s; // Yes. Return the signal.
    // ---------------------------------- //
    // 2. Divide Step: Divide the signal into even and odd samples.
    // ---------------------------------- //
    complexV_t<double> evenT(N / 2), oddT(N / 2);
    for (int i = 0; i < N / 2; i++) // For every sample up to N/2
    {                               // Perform decimation in time (DIT)
        evenT[i] = s[i * 2];        // Collect the even-indexed samples.
        oddT[i] = s[i * 2 + 1];     // Collect the odd-indexed samples.
    } // Done decimating in time.
    // ---------------------------------- //
    // 3. Conquer Step: Recursively apply the FFT to both halves.
    // ---------------------------------- //
    complexV_t<double> evenF = FFT(evenT); // Get the FFT of the even samples.
    complexV_t<double> oddF = FFT(oddT);   // Get the FFT of the odd samples.
    // ---------------------------------- //
    // Precompute the twiddle factors for a length N signal
    // ---------------------------------- //
    TwiddleFactor(N); // Precomputed twiddle-factors.
    // ---------------------------------- //
    // Compute the FFT butterfly from the even and odd parts.
    // ---------------------------------- //
    complexV_t<double> sigF(N);     // Where to store the entire signal.
    complex<double> t{0.0, 0.0};    // Twiddle factor coefficients.
    for (int k = 0; k < N / 2; k++) // For discrete-time samples k
    {                               // Perform the Cooley-Tukey DFT.
        // -------------------------------- //
        // Get amplitude and phase contribution for the current index.
        // -------------------------------- //
        t = twiddles[k] * oddF[k]; // Scale and change phase of odd input.
        // -------------------------------- //
        // Combine the results for the next stage in the butterfly.
        // -------------------------------- //
        sigF[k] = evenF[k] + t;         // Produce the even result.
        sigF[k + N / 2] = evenF[k] - t; // Produce the odd result.
    } // Done computing butterfly.
    // ---------------------------------- //
    // Step 6. Return the computed FFT.
    // ---------------------------------- //
    return sigF; // Return the spectrum of the signal.
} // ------------- FFT ---------------- //
// Windowed FFT
// Similar to vanilla FFT, but uses a spectral window to minimize "ringing" or
// spectral leakage that occurs outside due to the Gibbs phenomenon.
complexV_t<double> SignalOps::FFT( // ------------- FFT ---------------- //
    const complexV_t<double> &s,   // Our input signal.
    const WindowType &w)           // Our chosen window.
{                                  // ---------------------------------- //
    int N = s.size();              // The length of the input signal.
    // ---------------------------------- //
    // Generate the desired window to use on the signal.
    // ---------------------------------- //
    realV_t<double> window = GenerateWindow(w, N); // Get a window of this length.
    complexV_t<double> wsig(N);                    // Store windowed signal here.
    for (int i = 0; i < N; ++i)                    // For the duration of the signal, get
        wsig[i] = s[i] * window[i];                // the product of the windowed signal.
    // ---------------------------------- //
    // 1. Base Case: if the input size is 1, return the windowed signal.
    // ---------------------------------- //
    if (N <= 1)      // Is the signal smaller than one sample?
        return wsig; // Yes. Return windowed signal.
    // ---------------------------------- //
    // 2. Divide Step: Divide the signal into even and odd samples.
    // ---------------------------------- //
    complexV_t<double> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; i++) // For every sample up to N/2
    {                               // Perform decimation in time (DIT)
        even[i] = wsig[i * 2];      // Collect the even-indexed samples.
        odd[i] = wsig[i * 2 + 1];   // Collect the odd-indexed samples.
    } // Done decimating in time.
    // ---------------------------------- //
    // 3. Conquer Step: Recursively apply the FFT to both halves.
    // ---------------------------------- //
    complexV_t<double> evenS = FFT(even, w); // Get the FFT of the even samples.
    complexV_t<double> oddS = FFT(odd, w);   // Get the FFT of the odd samples.
    // ---------------------------------- //
    // Precompute the twiddle factors for a length N signal
    // ---------------------------------- //
    TwiddleFactor(N); // Precomputed twiddle factors.
    // ---------------------------------- //
    // Compute the FFT butterfly from the even and odd parts.
    // ---------------------------------- //
    complexV_t<double> sigF(N);     // Where to store the entire signal.
    complex<double> t{0.0, 0.0};    // Twiddle factor coefficients.
    for (int k = 0; k < N / 2; k++) // For discrete-time samples k
    {                               // Perform the Cooley-Tukey DFT.
        // -------------------------------- //
        // Get amplitude and phase contribution for the current index.
        // -------------------------------- //
        t = twiddles[k] * oddS[k]; // Scale and change phase of odd input.
        // -------------------------------- //
        // Combine the results for the next stage in the butterfly.
        // -------------------------------- //
        sigF[k] = evenS[k] + t;         // Produce the even result.
        sigF[k + N / 2] = evenS[k] - t; // Produce the odd result.
    } // Done computing butterfly.
    // ---------------------------------- //
    // Step 6. Return the computed FFT.
    // ---------------------------------- //
    return sigF; // Return the spectrum of the signal.
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ IFFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// So how do we calculate the Inverse FFT? To get the inverse of a signal X(k)
// of length k an elegant trick is performed.
// The IFFT is nothing more than the FFT multiplied by 1/N, and with a twiddle
// factor that rotates clockwise, instead of counter-clockwise.
// It is nothing more than the conjugate of the FFT multiplied by (1/N).
//
// The operation goes as follows:
// 1. Conjugate the discrete frequency input signal X(k):
//    X_conjugate(k) = Re(X(k)) - j*Im(X(k))
// 2. Next we perform the FFT on the conjugated signal, this performs the IDFT
//    but does not scale - that comes afterwards:
//    X(k) = SUM[n=0 to N-1] {x(n) * exp(-j*(2*pi/N)*k*n) }
// 3. Now we conjugate again and multiply by (1/N), this returns our signal to
//    the time domain: (wizardry at hand!)
//    x(n) = (1/N) * SUM[k=0 to N-1] * {X(k) * exp(-j*(2*pi/N)*k*n)
// 4. Return the time-domain samples of the signal.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
complexV_t<double> SignalOps::IFFT( // -------------- IFFT -------------- //
    const complexV_t<double> &s)    // Input frequency-domain signal.
{                                   // ---------------------------------- //
    int N = s.size();               // Get the length of the signal.
    // ---------------------------------- //
    // 1. Conjugate the input signal, just prepares the twiddle factors.
    // ---------------------------------- //
    complexV_t<double> sConj(N); // Conjugated signal buffer.
    for (int i = 0; i < N; ++i)  // For all discrete-freq samples
        sConj[i] = conj(s[i]);   // Conjugate the input signal.
    // ---------------------------------- //
    // 2. Perform FFT on conjugated signal.
    // ---------------------------------- //
    complexV_t<double> sigF = FFT(sConj); // The signal spectrum.
    // ---------------------------------- //
    // 3. Conjugate again essentially flipping the FFT's twiddle factor to positive.
    // ---------------------------------- //
    complexV_t<double> sigT(N);                           // Where to store the time domain signal.
    for (int i = 0; i < N; ++i)                           // For all discrete-freq samples
        sigT[i] = conj(sigF[i]) / static_cast<double>(N); // Conjugate and normalize.
    return sigT;                                          // Return time-domain signal.
} // ------------ IFFT ---------------- //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Convolution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Method to perform the convolution of two signals. Typically in our context
// the convolution is done between an input signal s(n) of length N  and
// filter of length N h(n).
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
complexV_t<double> SignalOps::Convolution( // Obtain the product of two signals.
    const complexV_t<double> &s,           // Our input signal.
    const complexV_t<double> &h)           // Our filter.
{                                          // ---------- Convolution ----------- //
    int n = s.size();                      // Length of the input signal.
    int m = h.size();                      // Length of the filter.
    int N = 1;                             // Size of the N-point FFT.
    while (N < n + m - 1)                  // While N is less than the convolution
        N <<= 1;                           // Find smallest power of 2 >= n+m-1
    // ---------------------------------- //
    // Zero-Pad the signal and filter to the length of the N-Point FFT.
    // ---------------------------------- //
    complexV_t<double> sPad = ZeroPad(s); // Zero pad the input signal.
    complexV_t<double> hPad = ZeroPad(h); // Zero pad the filter.
    // ---------------------------------- //
    // Apply the FFT on the Zero-Padded signal.
    // ---------------------------------- //
    complexV_t<double> S = FFT(sPad); // The FFT of the input signal.
    complexV_t<double> H = FFT(hPad); // The FFT of the filter.
    // ---------------------------------- //
    // Now the filtered signal is just the product of their spectrum.
    // ---------------------------------- //
    complexV_t<double> Y(N);    // Place to store resulting spectrum.
    for (int i = 0; i < N; ++i) // For the N-Points of the FFT.
        Y[i] = S[i] * H[i];     // Get the filtered spectrum.
    // ---------------------------------- //
    // Obtain the time-domain resulting signal using the IFFT.
    // ---------------------------------- //
    complexV_t<double> y(N); // Place to store resulting signal.
    y = IFFT(Y);             // Get the resulting signal.
    y.resize(n + m - 1);     // Truncate to the original size.
    return y;                // Return the filtered signal.
} // ---------- Convolution ----------- //
// Short-Time Fourier Transform Overlap Method
matrix_t<double> SignalOps::STFT(const complexV_t<double> &s, const WindowType &w, int wSiz, int overlap)
{
    // Calculate the step size based on the overlap percentage
    int step = wSiz * (1 - overlap / 100.0);

    // Calculate the number of segs needed to cover the entire signal
    int nSegs = (s.size() - wSiz + step) / step;

    // Initialize the resulting STFT matrix (each row is the FFT of a seg)
    matrix_t<double> sMat(nSegs, complexV_t<double>(wSiz));

    // Generate the window to be applied to each seg
    realV_t<double> window = GenerateWindow(windowType, wSiz);

    // Process each seg of the signal
    for (int i = 0; i < nSegs; ++i)
    {
        int start = i * step;              // Calculate the starting index of the seg
        complexV_t<double> seg(wSiz, 0.0); // Initialize the seg with zeros

        // Apply the window to the seg and copy the windowed signal
        // For the size of the window and the remaining signal seg
        for (int j = 0; j < wSiz && start + j < s.size(); ++j)
        {
            seg[j] = s[start + j] * window[j];
        }

        // Compute the FFT of the windowed seg and store it in the STFT matrix
        sMat[i] = FFT(seg, windowType);
    }

    // Return the STFT result (matrix of FFTs of the segs)
    return sMat;
}

// Inverse Short-Time-Fourier-Transform Method.
complexV_t<double> SignalOps::ISTFT(
    const matrix_t<double> &sMat, // Input STFT matrix
    const WindowType &w,          // Window type
    const int wSiz,               // Window size (should be power of 2)
    const int overlap             // Overlap percentage (e.g., 50%)
)
{
    // Calculate the step size based on the overlap percentage
    int step = wSiz * (1 - overlap / 100.0);

    // Calculate the length of the original signal from the STFT segs
    int len = step * (sMat.size() - 1) + wSiz;

    // Initialize the result signal and the overlap count for normalization
    complexV_t<double> sig(len, 0.0);    // Initialize the result signal
    realV_t<double> nOverlaps(len, 0.0); // Initialize the overlap count

    // Generate the window to be applied during the inverse process
    realV_t<double> window = GenerateWindow(w, wSiz);

    // Process each seg of the STFT matrix
    for (int i = 0; i < sMat.size(); ++i)
    {
        int start = i * step; // Calculate the starting index of the seg

        // Compute the IFFT of the current seg
        complexV_t<double> seg = IFFT(sMat[i]);

        // Overlap-add the IFFT result to the output signal
        for (int j = 0; j < wSiz && start + j < sig.size(); ++j)
        {
            sig[start + j] += seg[j] * window[j];
            nOverlaps[start + j] += window[j]; // Keep track of the windowing
        }
    }

    // Normalize the result by dividing by the overlap count
    for (int i = 0; i < sig.size(); ++i)
    {
        if (nOverlaps[i] != 0.0)
        {
            sig[i] /= nOverlaps[i];
        }
    }

    // Return the reconstructed time-domain signal
    return sig;
}
complexV_t<double> SignalOps::OLAProcessor(
    const complexV_t<double> &s,        // The input signal.
    const complexV_t<double> &h,        // The desired FIR filter.
    const WindowType &w,                // The window used.
    const int wSiz,                     // The size of the window.
    const int overlap)                  // The percentage of overlap
{                                       // ----------OLAProcessor----------- //
    matrix_t<double> sMat = STFT(s, w, wSiz, overlap); // STFT of input signal.
    complexV_t<double> H = FFT(h);      // FFT of the FIR filter.
    const int frames = sMat.size();     // Number of frames in the STFT matrix.
    matrix_t<double> sig(frames, complexV_t<double>(wSiz)); // The filtered signal.
    // -------------------------------- //
    // Perform element-wise multiplication of the STFTs
    // -------------------------------- //
    for (int i = 0; i < frames; ++i)    // For each frame in the STFT matrix.
    {                                   //
        for (int j = 0; j < wSiz; ++j)  // And for each phase-bin in the frame.
        {                               //
            sig[i][j] = sMat[i][j] * H[j];// Get the filtered signal.
        }                               // Done with each phase-bin
    }                                   // Done with each frame.
    // -------------------------------- //
    // Perform the inverse STFT to get the filtered time-domain signal.
    // -------------------------------- //
    return ISTFT(sig, w, wSiz, overlap); // Return the filtered signal.
} // ----------OLAProcessor----------- //

#endif // _FFT_
