/*
  ==============================================================================

    Biquad.cpp
    Created: 26 Jul 2018 7:50:20pm
    Author:  Frans Rosencrantz

  ==============================================================================
*/

#include "Biquad.h"
#include <stdlib.h>
#include <cmath>

void set_filter_type(biquad_context& biquad, int filterType, float Fs, int buffersize){
    biquad.Fs = Fs;
    biquad.filterType = filterType;
    biquad.bufferSize = buffersize;
    biquad.x = (float*)malloc(3*sizeof(float));
    biquad.y = (float*)malloc(3*sizeof(float));
    
    for (int i = 0; i < 3; i++) {
        biquad.x[i] = 0;
        biquad.y[i] = 0;
    }
}

void set_coefficient(biquad_context& biquad, float Q, float dBgain, float Fc, int filterType){
    biquad.Q = Q;
    biquad.dBgain = dBgain;
    biquad.Fc = Fc;
    biquad.filterType = filterType;
    
    if (biquad.filterType == 0){
        biquad.A = std::sqrt(std::pow(10.0f, biquad.dBgain/20.0f));
        biquad.w_0 = 2.0f * M_PI * (biquad.Fc/biquad.Fs);
        biquad.alpha = std::sin(biquad.w_0)/(2.0f*biquad.Q);
        
        biquad.b_0 = (1.0f - std::cos(biquad.w_0))/2.0f;
        biquad.b_1 = 1.0f - std::cos(biquad.w_0);
        biquad.b_2 = (1.0f - std::cos(biquad.w_0))/2.0f;
        
        biquad.a_0 = 1.0f + biquad.alpha;
        biquad.a_1 = -2.0f * std::cos(biquad.w_0);
        biquad.a_2 = 1.0f - biquad.alpha;
        
    }
    
    if (biquad.filterType == 1){
        biquad.A = std::sqrt(std::pow(10.0f, biquad.dBgain/20.0f));
        biquad.w_0 = 2.0f * M_PI * (biquad.Fc/biquad.Fs);
        biquad.alpha = std::sin(biquad.w_0)/(2.0f*biquad.Q);
        
        biquad.b_0 = (1.0f + std::cos(biquad.w_0))/2.0f;
        biquad.b_1 = -(1.0f + std::cos(biquad.w_0));
        biquad.b_2 = (1.0f + std::cos(biquad.w_0))/2.0f;
        
        biquad.a_0 = 1.0f + biquad.alpha;
        biquad.a_1 = -2.0f * std::cos(biquad.w_0);
        biquad.a_2 = 1.0f - biquad.alpha;
    }
    
}

void process_biquad(biquad_context& biquad, float &sample){
        biquad.x[0] = sample;
        
        biquad.y[0] = ((biquad.b_0/biquad.a_0) * biquad.x[0])
                    + ((biquad.b_1/biquad.a_0) * biquad.x[1])
                    + ((biquad.b_2/biquad.a_0) * biquad.x[2])
                    - ((biquad.a_1/biquad.a_0) * biquad.y[1])
                    - ((biquad.a_2/biquad.a_0) * biquad.y[2]);
        
        sample = biquad.y[0];
        biquad.y[2] = biquad.y[1];
        biquad.y[1] = biquad.y[0];
        biquad.x[2] = biquad.x[1];
        biquad.x[1] = biquad.x[0];
}

void process_audio_channels(biquad_context& biquad, float *buffer){
    for (int sample = 0; sample < biquad.bufferSize; sample++) {
        process_biquad(biquad, buffer[sample]);
    }
}

void get_magnitude_response(biquad_context& biquad, float* magnitude_response, int length){
    static float* sigma = (float*)malloc(length*sizeof(float));
    
    for (int i = 0; i < length; ++i) {
        sigma[i] = std::pow(std::sin(((float)i*(M_PI/(float)length))/2),2);
    }
    for (int i = 0; i < length; ++i) {
        magnitude_response[i] = (std::pow((biquad.b_0+biquad.b_1+biquad.b_2)/2.0f,2) - sigma[i]*(4*biquad.b_0*biquad.b_2*(1-sigma[i]) + biquad.b_1*(biquad.b_0+biquad.b_2))) / (std::pow((biquad.a_0+biquad.a_1+biquad.a_2)/2.0f,2) - sigma[i]*(4*biquad.a_0*biquad.a_2*(1-sigma[i]) + biquad.a_1*(biquad.a_0+biquad.a_2)));
        if (magnitude_response[i] > 0){
            magnitude_response[i] = std::log(magnitude_response[i]);
        }
        else{
            magnitude_response[i] = magnitude_response[i-1];
        }
    }
}

