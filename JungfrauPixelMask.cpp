#include <algorithm> // std::count
#include <stdint.h>
#include <fstream>
#include <math.h>
#include "jungfrau.h"

JungfrauPixelMask::JungfrauPixelMask(uint32_t *in_mask) : mask(in_mask) {
    std::fill_n(mask, CHANNELS_PER_WORKER, 0); // initialise with zeros

}

void JungfrauPixelMask::MaskIfGainNot(int expected, uint16_t *imagedata) {
    // if gain of a pixel is not as expected, mask pixel
    for (int i = 0; i < CHANNELS_PER_WORKER; i++) {
        if (((imagedata[i]&0xc000) >> 14) != expected) {
            if (expected == 0) mask[i] |= 2;
            else if (expected == 1) mask[i] |= 4;
            else if (expected == 3) mask[i] |= 8;
        }
    }
}

void JungfrauPixelMask::MaskIfPedeRMSG0GreaterThan(float* pedestalG0RMSs, int val) {
    for (int i = 0; i < CHANNELS_PER_WORKER; i++) {
        if (pedestalG0RMSs[i] > val) {
            mask[i] |= 16;
        }
    }
}

void JungfrauPixelMask::MaskIfPedestalStep(uint16_t *imagedata, float* pedestalG0, int step) {
    for (int i = 0; i < CHANNELS_PER_WORKER; i++) {
        if ( fabs(pedestalG0[i] - (imagedata[i]&0x3fff)) > step) {
                      mask[i] |= 8;
        }
    }    
}


//  int jungfrauPixelMask::getNMasked(int* mask) {
//    return count(mask, mask + NCH, 1);
//  }
