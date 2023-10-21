#include "fir_filter.h"

// FIR3 Filter Function Implementation
void fir3(fir3Vars* a) {
    word input;
    for (size_t i = 0; i < a->Ain.size(); i++) {
        input = a->Ain[i];
        a->Aout[i] = a->b0 * input + a->b1 * a->s1 + a->b2 * a->s2;
        a->s2 = a->s1;
        a->s1 = input;
    }
}
