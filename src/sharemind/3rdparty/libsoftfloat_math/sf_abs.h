#ifndef SF_ABS_H
#define SF_ABS_H

#include <sharemind/3rdparty/libsoftfloat/softfloat.h>
#include <sharemind/extern_c.h>


SHAREMIND_EXTERN_C_BEGIN

sf_float32 sf_float32_abs(sf_float32 x);
sf_float64 sf_float64_abs(sf_float64 x);

SHAREMIND_EXTERN_C_END

#endif /* SF_ABS_H */
