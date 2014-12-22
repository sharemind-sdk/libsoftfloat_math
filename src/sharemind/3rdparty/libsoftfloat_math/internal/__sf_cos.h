#ifndef INTERNAL___SF_COS_H
#define INTERNAL___SF_COS_H

#include <sharemind/3rdparty/libsoftfloat/softfloat.h>


sf_result32f __sf_float32_cos(sf_float64 x, sf_fpu_state fpu);
sf_result64f __sf_float64_cos(sf_float64 x, sf_float64 y, sf_fpu_state fpu);

#endif /* INTERNAL___SF_COS_H */
