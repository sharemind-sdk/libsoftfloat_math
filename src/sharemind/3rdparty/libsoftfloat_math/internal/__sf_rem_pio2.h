#ifndef INTERNAL___SF_REM_PIO2_H
#define INTERNAL___SF_REM_PIO2_H

#include <sharemind/3rdparty/libsoftfloat/softfloat.h>
#include <stdint.h>


sf_result32i __sf_float32_rem_pio2(sf_float32 x, sf_float64 * y, sf_fpu_state fpu);
sf_result32i __sf_float64_rem_pio2(sf_float64 x, sf_float64 * y, sf_fpu_state fpu);

#endif /* INTERNAL___SF_REM_PIO2_H */
