#ifndef SF_EXP_H
#define SF_EXP_H

#include <sharemind/3rdparty/libsoftfloat/softfloat.h>
#include <sharemind/extern_c.h>


SHAREMIND_EXTERN_C_BEGIN

sf_result32f sf_float32_exp(sf_float32 x, sf_fpu_state fpu);
sf_result64f sf_float64_exp(sf_float64 x, sf_fpu_state fpu);

SHAREMIND_EXTERN_C_END

#endif /* SF_EXP_H */
