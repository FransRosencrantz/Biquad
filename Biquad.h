/*
  ==============================================================================

    Biquad.h
    Created: 26 Jul 2018 7:50:20pm
    Author:  Frans Rosencrantz

  ==============================================================================
*/

#pragma once

struct biquad_context {
    float b_0 = 0.0f;
    float b_1 = 0.0f;
    float b_2 = 0.0f;
    float alpha = 0.0f;
    float a_0 = 0.0f;
    float a_1 = 0.0f;
    float a_2 = 0.0f;
    float Fs = 44100;
    float Fc = 500;
    float dBgain = 0.0f;
    float A = 0.0f;
    float Q = 1.0f;
    float w_0 = 0.0f;
    int filterType = 0;
    int bufferSize = 0;
    float *x;
    float *y;
    float factor = 1;
    float N0, N1, N2 = 1;
};

void set_filter_type(biquad_context&, int filterType, float Fs, int buffersize);

void set_coefficient(biquad_context&, float Q, float dBgain, float Fc, int filterType);

void get_magnitude_response(biquad_context&, float* magnitude_response, int length);

void process_audio_channels(biquad_context&, float *buffer);
