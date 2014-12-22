#ifndef INTERNAL___SF_SIN_H
#define INTERNAL___SF_SIN_H

#include <sharemind/3rdparty/libsoftfloat/softfloat.h>
#include <stdint.h>


sf_result32f __sf_float32_sin(sf_float64 x, sf_fpu_state fpu);
sf_result64f __sf_float64_sin(sf_float64 x, sf_float64 y, int64_t iy, sf_fpu_state fpu);

#endif /* INTERNAL___SF_SIN_H */
