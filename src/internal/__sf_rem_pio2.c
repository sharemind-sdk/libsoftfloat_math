/*
 * Copyright (c) 2015, Cybernetica AS
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *   * Neither the name of Cybernetica AS nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL CYBERNETICA AS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* origin: musl-libc src/math/__rem_pio2.c */
/* origin: musl-libc src/math/__rem_pio2f.c */
/*
 * ----------------------------------------------------------------------
 * Copyright © 2005-2014 Rich Felker, et al.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ----------------------------------------------------------------------
 */
/* origin: FreeBSD /usr/src/lib/msun/src/e_rem_pio2f.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 * Debugged and optimized by Bruce D. Evans.
 */
/* origin: FreeBSD /usr/src/lib/msun/src/e_rem_pio2.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 *
 * Optimized by Bruce D. Evans.
 */
/* __rem_pio2(x,y)
 *
 * return the remainder of x rem pi/2 in y[0]+y[1]
 * use __rem_pio2_large() for large x
 */

#include "__sf_rem_pio2.h"
#include "__sf_rem_pio2_large.h"
#include "macros.h"


/*
 * invpio2:  53 bits of 2/pi
 */
static const sf_float64
toint   = 0x4338000000000000, /* 1.5 / 2.22044604925031308085e-16 */
invpio2 = 0x3fe45f306dc9c883; /* 6.36619772367581382433e-01       */

/*
 * pio2_1:   first 25 bits of pi/2
 * pio2_1t:  pi/2 - pio2_1
 */
static const sf_float64
f32pio2_1  = 0x3ff921fb50000000, /* 1.57079631090164184570e+00 */
f32pio2_1t = 0x3e5110b4611a6263; /* 1.58932547735281966916e-08 */

/*
 * pio2_1:   first  33 bit of pi/2
 * pio2_1t:  pi/2 - pio2_1
 * pio2_2:   second 33 bit of pi/2
 * pio2_2t:  pi/2 - (pio2_1+pio2_2)
 * pio2_3:   third  33 bit of pi/2
 * pio2_3t:  pi/2 - (pio2_1+pio2_2+pio2_3)
 */
static const sf_float64
f64pio2_1  = 0x3ff921fb54400000, /* 1.57079632673412561417e+00       */
f64pio2_1t = 0x3dd0b4611a626331, /* 6.07710050650619224932e-11       */
f64pio2_2  = 0x3dd0b4611a600000, /* 6.07710050630396597660e-11       */
f64pio2_2t = 0x3ba3198a2e037073, /* 2.02226624879595063154e-21       */
f64pio2_3  = 0x3ba3198a2e000000, /* 2.02226624871116645580e-21       */
f64pio2_3t = 0x397b839a252049c1; /* 8.47842766036889956997e-32       */

sf_result32i __sf_float32_rem_pio2(sf_float32 x, sf_float64 * y, sf_fpu_state fpu) {
    sf_float64 fn, tx[1], ty[1], tmp1, tmp2;
    sf_int32 n;
    uint32_t ix, sign;
    int64_t e0;

    ix = x & 0x7fffffff;
    /* 25+53 bit pi is good enough for medium size */
    if (ix < 0x4dc90fdb) {  /* |x| ~< 2^28*(pi/2), medium size */
        /* Use a specialized rint() to get fn.  Assume round-to-nearest. */
        { /* fn = x*invpio2 + toint - toint; */
            SF_F64_OP(fpu, fn, sf_float32_to_float64, x);
            SF_F64_OP(fpu, fn, sf_float64_mul, fn, invpio2);
            SF_F64_OP(fpu, fn, sf_float64_add, fn, toint);
            SF_F64_OP(fpu, fn, sf_float64_sub, fn, toint);
        }
        SF_I32_OP(fpu, n, sf_float64_to_int32_round_to_zero, fn); /* n  = (int32_t)fn; */
        { /* *y = x - fn*pio2_1 - fn*pio2_1t; */
            SF_F64_OP(fpu, tmp1, sf_float32_to_float64, x);

            SF_F64_OP(fpu, tmp2, sf_float64_mul, fn, f32pio2_1);
            SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, tmp2);

            SF_F64_OP(fpu, tmp2, sf_float64_mul, fn, f32pio2_1t);
            SF_F64_OP(fpu, *y, sf_float64_sub, tmp1, tmp2);
        }
        return (sf_result32i) { n, fpu };
    }
    if (ix >= 0x7f800000) {  /* x is inf or NaN */
        { /* *y = x-x; */
            SF_F64_OP(fpu, tmp1, sf_float32_to_float64, x);
            SF_F64_OP(fpu, *y, sf_float64_sub, tmp1, tmp1);
        }
        return (sf_result32i) { 0, fpu };
    }
    /* scale x into [2^23, 2^24-1] */
    sign = x >> 31;
    e0 = (ix >> 23) - (0x7f + 23);  /* e0 = ilogb(|x|)-23, positive */
    SF_F64_OP(fpu, tx[0], sf_float32_to_float64, (sf_bits32)ix - (sf_bits32)(e0 << 23));
    SF_I32_OP(fpu, n, __sf_float64_rem_pio2_large, tx, ty, e0, 1, 0); /* n  =  __sf_float64_rem_pio2_large(tx, ty, e0, 1, 0); */
    if (sign) {
        *y = sf_float64_neg(ty[0]); /* *y = -ty[0]; */
        return (sf_result32i) { -n, fpu };
    }
    *y = ty[0];
    return (sf_result32i) { n, fpu };
}

/* caller must handle the case when reduction is not needed: |x| ~<= pi/4 */
sf_result32i __sf_float64_rem_pio2(sf_float64 x, sf_float64 * y, sf_fpu_state fpu) {
    sf_float64 z, w, t, r, fn, tmp1, tmp2;
    sf_float64 tx[3], ty[2];
    uint64_t ix, sign;
    int64_t ex, ey, i;
    sf_int32 i32tmp;
    sf_int32 n;

    sign = x >> 63;
    ix = x & 0x7fffffff00000000;
    if (ix <= 0x400f6a7a00000000) {  /* |x| ~<= 5pi/4 */
        if ((ix & 0xfffff00000000) == 0x921fb00000000)  /* |x| ~= pi/2 or 2pi/2 */
            goto medium;  /* cancellation -- use medium case */
        if (ix <= 0x4002d97c00000000) {  /* |x| ~<= 3pi/4 */
            if (!sign) {
                SF_F64_OP(fpu, z, sf_float64_sub, x, f64pio2_1); /* z = x - pio2_1; */ /* one round good to 85 bits */
                SF_F64_OP(fpu, y[0], sf_float64_sub, z, f64pio2_1t); /* y[0] = z - pio2_1t; */
                { /* y[1] = (z-y[0]) - pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_sub, y[1], f64pio2_1t);
                }
                return (sf_result32i) { 1, fpu };
            } else {
                SF_F64_OP(fpu, z, sf_float64_add, x, f64pio2_1); /* z = x + pio2_1; */
                SF_F64_OP(fpu, y[0], sf_float64_add, z, f64pio2_1t); /* y[0] = z + pio2_1t; */
                { /* y[1] = (z-y[0]) + pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_add, y[1], f64pio2_1t);
                }
                return (sf_result32i) { -1, fpu };
            }
        } else {
            if (!sign) {
                { /* z = x - 2*pio2_1; */
                    SF_F64_OP(fpu, z, sf_float64_mul, 0x4000000000000000, f64pio2_1);
                    SF_F64_OP(fpu, z, sf_float64_sub, x, z);
                }
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4000000000000000, f64pio2_1t);
                SF_F64_OP(fpu, y[0], sf_float64_sub, z, tmp1); /* y[0] = z - 2*pio2_1t; */
                { /* y[1] = (z-y[0]) - 2*pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_sub, y[1], tmp1);
                }
                return (sf_result32i) { 2, fpu };
            } else {
                { /* z = x + 2*pio2_1; */
                    SF_F64_OP(fpu, z, sf_float64_mul, 0x4000000000000000, f64pio2_1);
                    SF_F64_OP(fpu, z, sf_float64_add, x, z);
                }
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4000000000000000, f64pio2_1t);
                SF_F64_OP(fpu, y[0], sf_float64_add, z, tmp1); /* y[0] = z + 2*pio2_1t; */
                { /* y[1] = (z-y[0]) + 2*pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_add, y[1], tmp1);
                }
                return (sf_result32i) { -2, fpu };
            }
        }
    }
    if (ix <= 0x401c463b00000000) {  /* |x| ~<= 9pi/4 */
        if (ix <= 0x4015fdbc00000000) {  /* |x| ~<= 7pi/4 */
            if (ix == 0x4012d97c00000000)  /* |x| ~= 3pi/2 */
                goto medium;
            if (!sign) {
                { /* z = x - 3*pio2_1; */
                    SF_F64_OP(fpu, z, sf_float64_mul, 0x4008000000000000, f64pio2_1);
                    SF_F64_OP(fpu, z, sf_float64_sub, x, z);
                }
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4008000000000000, f64pio2_1t);
                SF_F64_OP(fpu, y[0], sf_float64_sub, z, tmp1); /* y[0] = z - 3*pio2_1t; */
                { /* y[1] = (z-y[0]) - 3*pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_sub, y[1], tmp1);
                }
                return (sf_result32i) { 3, fpu };
            } else {
                { /* z = x + 3*pio2_1; */
                    SF_F64_OP(fpu, z, sf_float64_mul, 0x4008000000000000, f64pio2_1);
                    SF_F64_OP(fpu, z, sf_float64_add, x, z);
                }
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4008000000000000, f64pio2_1t);
                SF_F64_OP(fpu, y[0], sf_float64_add, z, tmp1); /* y[0] = z + 3*pio2_1t; */
                { /* y[1] = (z-y[0]) + 3*pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_add, y[1], tmp1);
                }
                return (sf_result32i) { -3, fpu };
            }
        } else {
            if (ix == 0x401921fb00000000)  /* |x| ~= 4pi/2 */
                goto medium;
            if (!sign) {
                { /* z = x - 4*pio2_1; */
                    SF_F64_OP(fpu, z, sf_float64_mul, 0x4010000000000000, f64pio2_1);
                    SF_F64_OP(fpu, z, sf_float64_sub, x, z);
                }
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4010000000000000, f64pio2_1t);
                SF_F64_OP(fpu, y[0], sf_float64_sub, z, tmp1); /* y[0] = z - 4*pio2_1t; */
                { /* y[1] = (z-y[0]) - 4*pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_sub, y[1], tmp1);
                }
                return (sf_result32i) { 4, fpu };
            } else {
                { /* z = x + 4*pio2_1; */
                    SF_F64_OP(fpu, z, sf_float64_mul, 0x4010000000000000, f64pio2_1);
                    SF_F64_OP(fpu, z, sf_float64_add, x, z);
                }
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4010000000000000, f64pio2_1t);
                SF_F64_OP(fpu, y[0], sf_float64_add, z, tmp1); /* y[0] = z + 4*pio2_1t; */
                { /* y[1] = (z-y[0]) + 4*pio2_1t; */
                    SF_F64_OP(fpu, y[1], sf_float64_sub, z, y[0]);
                    SF_F64_OP(fpu, y[1], sf_float64_add, y[1], tmp1);
                }
                return (sf_result32i) { -4, fpu };
            }
        }
    }
    if (ix < 0x413921fb00000000) {  /* |x| ~< 2^20*(pi/2), medium size */
medium:
        /* rint(x/(pi/2)), Assume round-to-nearest. */
        { /* fn = x*invpio2 + toint - toint; */
            SF_F64_OP(fpu, fn, sf_float64_mul, x, invpio2);
            SF_F64_OP(fpu, fn, sf_float64_add, fn, toint);
            SF_F64_OP(fpu, fn, sf_float64_sub, fn, toint);
        }
        SF_I32_OP(fpu, n, sf_float64_to_int32_round_to_zero, fn); /* n = (int32_t)fn; */
        { /* r = x - fn*pio2_1; */
            SF_F64_OP(fpu, r, sf_float64_mul, fn, f64pio2_1);
            SF_F64_OP(fpu, r, sf_float64_sub, x, r);
        }
        SF_F64_OP(fpu, w, sf_float64_mul, fn, f64pio2_1t); /* w = fn*pio2_1t; */ /* 1st round, good to 85 bits */
        SF_F64_OP(fpu, y[0], sf_float64_sub, r, w); /* y[0] = r - w; */
        ey = y[0] >> 52 & 0x7ff;
        ex = (int64_t)ix >> 52;
        if (ex - ey > 16) { /* 2nd round, good to 118 bits */
            t = r;
            SF_F64_OP(fpu, w, sf_float64_mul, fn, f64pio2_2); /* w = fn*pio2_2; */
            SF_F64_OP(fpu, r, sf_float64_sub, t, w); /* r = t - w; */
            { /* w = fn*pio2_2t - ((t-r)-w); */
                SF_F64_OP(fpu, tmp1, sf_float64_mul, fn, f64pio2_2t);
                SF_F64_OP(fpu, tmp2, sf_float64_sub, t, r);
                SF_F64_OP(fpu, tmp2, sf_float64_sub, tmp2, w);
                SF_F64_OP(fpu, w, sf_float64_sub, tmp1, tmp2);
            }
            SF_F64_OP(fpu, y[0], sf_float64_sub, r, w); /* y[0] = r - w; */
            ey = y[0] >> 52 & 0x7ff;
            if (ex - ey > 49) {  /* 3rd round, good to 151 bits, covers all cases */
                t = r;
                SF_F64_OP(fpu, w, sf_float64_mul, fn, f64pio2_3); /* w = fn*pio2_3; */
                SF_F64_OP(fpu, r, sf_float64_sub, t, w); /* r = t - w; */
                { /* w = fn*pio2_3t - ((t-r)-w); */
                    SF_F64_OP(fpu, tmp1, sf_float64_mul, fn, f64pio2_3t);
                    SF_F64_OP(fpu, tmp2, sf_float64_sub, t, r);
                    SF_F64_OP(fpu, tmp2, sf_float64_sub, tmp2, w);
                    SF_F64_OP(fpu, w, sf_float64_sub, tmp1, tmp2);
                }
                SF_F64_OP(fpu, y[0], sf_float64_sub, r, w); /* y[0] = r - w; */
            }
        }
        { /* y[1] = (r - y[0]) - w; */
            SF_F64_OP(fpu, y[1], sf_float64_sub, r, y[0]);
            SF_F64_OP(fpu, y[1], sf_float64_sub, y[1], w);
        }
        return (sf_result32i) { n, fpu };
    }
    /*
     * all other (large) arguments
     */
    if (ix >= 0x7ff0000000000000) {  /* x is inf or NaN */
        { /* y[0] = y[1] = x - x; */
            SF_F64_OP(fpu, y[0], sf_float64_sub, x, x);
            SF_F64_OP(fpu, y[1], sf_float64_sub, x, x);
        }
        return (sf_result32i) { 0, fpu };
    }
    /* set z = scalbn(|x|,-ilogb(x)+23) */
    z = x;
    z &= (uint64_t)-1 >> 12;
    z |= (uint64_t)(0x3ff + 23) << 52;
    for (i = 0; i < 2; i++) {
        { /* tx[i] = (double)(int32_t)z; */
            SF_I32_OP(fpu, i32tmp, sf_float64_to_int32_round_to_zero, z);
            tx[i] = sf_int32_to_float64(i32tmp);
        }
        { /* z = (z-tx[i])*0x1p24; */
            SF_F64_OP(fpu, z, sf_float64_sub, z, tx[i]);
            SF_F64_OP(fpu, z, sf_float64_mul, z, 0x4170000000000000);
        }
    }
    tx[i] = z;
    /* skip zero terms, first term is non-zero */
    while (tx[i] == 0x0)
        i--;
    { /* n = __rem_pio2_large(tx,ty,(int)(ix>>52)-(0x3ff+23),i+1,1); */
        SF_I32_OP(fpu, n, __sf_float64_rem_pio2_large, tx, ty,
                  (int64_t)(ix >> 52) - (0x3ff + 23), i + 1, 1);
    }
    if (sign) {
        y[0] = sf_float64_neg(ty[0]); /* y[0] = -ty[0]; */
        y[1] = sf_float64_neg(ty[1]); /* y[1] = -ty[1]; */
        return (sf_result32i) { -n, fpu };
    }
    y[0] = ty[0];
    y[1] = ty[1];
    return (sf_result32i) { n, fpu };
}
