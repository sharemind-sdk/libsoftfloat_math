#ifndef INTERNAL___SF_REM_PIO2_LARGE_H
#define INTERNAL___SF_REM_PIO2_LARGE_H

#include <sharemind/3rdparty/libsoftfloat/softfloat.h>
#include <stdint.h>


sf_result32i __sf_float64_rem_pio2_large(sf_float64 * x, sf_float64 * y,
        int64_t e0, int64_t nx, int64_t prec, sf_fpu_state fpu);

#endif /* INTERNAL___SF_REM_PIO2_LARGE_H */
