#include "AudioHandler.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <valarray>
#include <numeric>

#include <cstdint>
#include "transforms.h"
#include "Statistical.h"
#include "Audiohandler.h"
using namespace std; // Added for convenience. You can also prefix with std:: wherever required.

void AudioHandler::readAudio(const std::string& filepath) {
    ifstream infile(filepath, ios::binary);

    if (!infile) {
        cerr << "Error opening file." << endl;
        return;
    }

    // Read header
    infile.read(reinterpret_cast<char*>(&meta), sizeof(Header));

    // Read audio data
    short int value;
    while (infile.read(reinterpret_cast<char*>(&value), sizeof(short int))) {
        double normalized = value / 32768.0;
        audioData.push_back(normalized);
    }
}

int AudioHandler::getSampleRate() const {
    return meta.sample_rate;
}

void AudioHandler::writeToFile(const std::string& outFile) {  // Added the 'AudioHandler::' qualifier.
    ofstream outfile(outFile, ios::binary);

    // Write the header
    outfile.write(reinterpret_cast<char*>(&meta), sizeof(Header));

    // Write audio data
    for (const auto& val : audioData) {
        short int sample = static_cast<short int>(val * 32768.0);
        outfile.write(reinterpret_cast<char*>(&sample), sizeof(short int));
    }
}

void AudioHandler::setAudioData(const std::vector<double>& data) {
    audioData = data;
}

std::vector<double>& AudioHandler::getAudioData() {
    return audioData;
}

int AudioHandler::padToNextPowerOfTwo(int sample_rate) {
    // Calculate the next power of 2 greater than the current size
    int nextPowerOfTwo = 1;
    while (nextPowerOfTwo < audioData.size()) {
        nextPowerOfTwo *= 2;
    }

    return nextPowerOfTwo;

    // Resize the vector. New elements are initialized to 0.
    audioData.resize(nextPowerOfTwo);
}

// bit reversal permutation
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


//  std::vector<std::complex<double>> weiner(const std::vector<std::complex<double>>& Y, int n_fft, double Fs, double alpha, double beta)
int main() {
    Statistical stats;
	AudioHandler audio;
    
    int f_sample = audio.getSampleRate();
    int n_fft = audio.padToNextPowerOfTwo(f_sample);
    const std::vector<std::complex<double>> x; // = audio.getAudioData();
	audio.readAudio("C:\\Users\\User\\Desktop\\Baby.wav");
    audio.getSampleRate();
	// int f_upsample = 2 * audio.padToNextPowerOfTwo(f_sample);

    // Get audio data
    std::vector<double> audioData = audio.getAudioData();
   
    // Convert real audio data to complex
    std::vector<std::complex<double>> x(audioData.begin(), audioData.end());
    fft_stride(x);
    std::vector<std::complex<double>> Y = stats.weiner(x, n_fft, f_sample, 0.5, 0.5);

	audio.writeToFile("C:\\Users\\User\\Desktop\\Baby_Anal2_C_Weiner.wav");
	return 0;
}