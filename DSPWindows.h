#ifndef DSP_WINDOWS_H
#define DSP_WINDOWS_H

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <valarray>
#include <cstdint>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif
// Missing Theory
using namespace std;

template <typename T>
class Window
{
 friend class SpectralOps;
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
        LowSideLobe,
        Rectangular
    };

    Window(const int N);
    Window(const WindowType &w, const int N);
    ~Window(void);
    void SetWindowType(const WindowType &w, const int N);
    vector<T> Hanning(const int N);
    vector<T> Hamming(const int N);
    vector<T> BlackmanHarris(const int N);
    vector<T> ExactBlackman(const int N);
    vector<T> Blackman(const int N);
    vector<T> FlatTop(const int N);
    vector<T> FourTermBHarris(const int N);
    vector<T> SevenTermBHarris(const int N);
    vector<T> LowSideLobe(const int N);
    vector<T> Rectangular(const int N);
    ~Window(void);
    vector<T> GenerateWindow(const WindowType &w, const int N);
    
    // Accessors
   const inline Window<T> GetWindow (void) const { return window; }
   const inline int GetWindowSize (void) const { return windowSize; }
   const inline vector<T> GetDefaultWindow (void) { return Rectangular(windowSize);}  
   inline void SetWindowSize (const int wSiz) {windowSize=wsiz;}
   

private:
    void SetWindowType (const WindowType& w, const int N);
    int windowSize;
    WindowType window;
};

// ------------------------------------ //
// Constructors and Destructors
// ------------------------------------ //

template <typename T>
Window<T>::Window(const int N) : windowSize{N}, window{WindowType::Rectangular} {}

template <typename T>
Window<T>::Window(const WindowType &w, const int N) 
{
  SetWindowType(w,N);
}


template <typename T>
Window<T>::~Window(void) {}

// ------------------------------------ //
// Accessor Methods
// ------------------------------------ //

template <typename T>
const int Window<T>::GetWindowSize(void) const
{
    return windowSize;
}

template <typename T>
void Window<T>::SetWindowSize(const int N)
{
    windowSize = N;
}

template <typename T>
vector<T> Window<T>::GetWindowType(void) const
{
    return window;
}

template <typename T>
void Window<T>::SetWindowType(const WindowType &w, const int N)
{
    switch (w):
        {
        case:
            WindowType::Hanning : window = Hanning(N);
            break;
        case:
            WindowType::Hamming : window = Hamming(N);
            break;
        case:
            WindowType::BlackmanHarris : window = BlackmanHarris(N);
        case:
            WindowType::ExactBlackman : window = ExactBlackman(N);
            break;
        case:
            WindowType::Blackman : window = Blackman(N);
            break;
        case:
            WindowType::FlatTop : window = FlatTop(N);
            break;
        case:
            WindowType::FourTermBHarris : window = FourTermBHarris(N);
            break;
        case:
            WindowType::SevenTermBHarris : window = SevenTermBHarris(N);
            break;
        case:
            WindowType::LowSideLobe : window = LowSideLobe(N);
            break;
        case:
            WindowType::Rectangular : window = Rectangular(N);
            break;
        default:
            window = Rectangular(N);
            break;
        }
}

template <typename T>
vector<T> Window<T>::GetDefaultWindow(void)
{
    return Rectangular(windowSize);
}

// ------------------------------------ //
// Window Definition Methods
// Reference: https://en.wikipedia.org/wiki/Window_function
// https://web.archive.org/web/20050113013738id_/http://congres.cran.uhp-nancy.fr/ICASSP_2001/MAIN/papers/pap45.pdf
// https://www.ni.com/docs/en-US/bundle/labwindows-cvi/page/advancedanalysisconcepts/lvac_low_sidelobe.html?srsltid=AfmBOoq24bE811jsNCA5Frywall7E4fABxA6kj3FgSxqYY_808W37dA1

// ------------------------------------ //
template <typename T>
vector<T> Window<T>::Hanning(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
        w[n] = 0.5 * (1 - cos(2 * M_PI * n / (N - 1)));
    return w;
}

template <typename T>
vector<T> Window<T>::Hamming(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
        w[n] = 0.5383553946707251 - 0.4616446053292749 * cos(2 * M_PI * n / (N - 1));
    return w;
}

template <typename T>
vector<T> Window<T>::BlackmanHarris(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
        w[n] = 0.35875 - 0.48829 * cos(2 * M_PI * n / (N - 1)) + 0.14128 * cos(4 * M_PI * n / (N - 1)) - 0.01168 * cos(6 * M_PI * n / (N - 1));
    return w;
}
template <typename T>
vector<T> Window<T>::ExactBlackman(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
        w[i] = 0.4243800934609435 - 0.4973406350967378 * cos(2 * M_PI * n / (N - 1)) + 0.7827927144231873 * cos(4 * M_PI * n / (N - 1));
    return w;
}

template <typename T>
vector<T> Window<T>::Blackman(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
        w[n] = 0.42 - 0.5 * cos(2 * M_PI * n / (N - 1)) + 0.08 * cos(4 * M_PI * n / (N - 1));
    return w;
}

template <typename T>
vector<T> Window<T>::FlatTop(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
        w[n] = 0.21557895 - 0.41663158 * cos(2 * M_PI * n / (N - 1)) + 0.277263158 * cos(4 * M_PI * n / (N - 1)) - 0.083578947 * cos(6 * M_PI * n / (N - 1)) + 0.006947368 * cos(8 * M_PI * n / (N - 1));
    return w;
}

template <typename T>
vector<T> Window<T>::FourTermBHarris(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
    {
        w[n] = 0.3635819267707608 - 0.4891774371450171 * cos(2 * M_PI * n / (N - 1)) + 0.1365995139786921 * cos(4 * M_PI * n / (N - 1)) - 0.01064112210553003 * cos(6 * M_PI * n / (N - 1));
    }
    return w;
}

template <typename T>
vector<T> Window<T>::SevenTermBHarris(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
    {
        w[n] = 0.27105140069342 - 0.43329793923448 * cos(2 * M_PI * n / (N - 1)) + 0.21812299954311 * cos(4 * M_PI * n / (N - 1)) - 0.06592544638803 * cos(6 * M_PI * n / (N - 1)) + 0.01081174209837 * cos(8 * M_PI * n / (N - 1)) - 0.00077658482522 * cos(10 * M_PI * n / (N - 1)) + 0.00001388721735 * cos(12 * M_PI * n / (N - 1));
    }
    return w;
}

template <typename T>
vector<T> Window<T>::LowSideLobe(const int N)
{
    vector<T> w(N, T(0));
    for (int n = 0; n < N; ++n)
    {
        w[n] = 0.471492057 - 0.17553428 * cos(2 * M_PI * n / (N - 1)) + 0.028497078 * cos(4 * M_PI * n / (N - 1)) - 0.001261367 * cos(6 * M_PI * n / (N - 1));
    }
    return w;
}

template <typename T>
vector<T> Window<T>::Rectangular(const int N)
{
    vector<T> w(N, T(1));
    return w;
}

template <typename T>
vector<T> Window<T>::GenerateWindow(const WindowType& w, const int N)
{
  SetWindowType(w,N);
  return window;
}
#endif // DSP_WINDOWS_H
