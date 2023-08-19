// AudioHandler.h

#ifndef AUDIO_HANDLER_H
#define AUDIO_HANDLER_H

#include <vector>
#include <string>

// WAVE PCM soundfile format
struct Header {
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short int audio_format;
    short int num_channels;
    int sample_rate;
    int byte_rate;
    short int block_align;
    short int bits_per_sample;
    char subchunk2_id[4];
    int subchunk2_size;
};

class AudioHandler {
private:
    Header meta;
    std::vector<double> audioData;

public:
    void readAudio(const std::string& filepath);
    void writeToFile(const std::string& outFile);

    std::vector<double>& getAudioData();
    void setAudioData(const std::vector<double>& data);
    int getSampleRate() const;
    int padToNextPowerOfTwo(int sample_rate);
};

#endif // AUDIO_HANDLER_