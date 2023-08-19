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
#include <fftw3.h>
#include <cstdint>
#include "transforms.h"

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
