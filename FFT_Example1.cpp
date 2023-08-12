#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <sciplot/sciplot.hpp>
#include "transforms.h" // Include the header file for the Transform class

using namespace std;

int main() {
    // Sample input sequence of 8 numbers
    vector<complex<double>> dataIn = { 1, 2, 3, 4, 5, 6, 7, 8 };

    // Create an instance of the Transform class
    transforms dsp;

    // Call the fft method from the Transform class on the input data
    dsp.fft(dataIn);

    

    // Print the results of the FFT
    for (auto& val : dataIn) { // C++11 range-based for loop
        cout << val << endl;
        
        

    }
    // IFFT portion 
    cout << "\n" << endl;
    // Call the hermitian to perform the firrst step of the IFFT()
    vector<complex<double>> hermit = dsp.hermitian(dataIn);
    dsp.ifft(hermit);

    for (auto& val : hermit) {
        cout << val << endl;
        
    }
        
      return 0;


}
