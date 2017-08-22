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

/* origin: musl-libc src/math/sin.c */
/* origin: musl-libc src/math/sinf.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/s_sinf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 * Optimized by Bruce D. Evans.
 */
/* origin: FreeBSD /usr/src/lib/msun/src/s_sin.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/* sin(x)
 * Return sine function of x.
 *
 * kernel function:
 *      __sin            ... sine function on [-pi/4,pi/4]
 *      __cos            ... cose function on [-pi/4,pi/4]
 *      __rem_pio2       ... argument reduction routine
 *
 * Method.
 *      Let S,C and T denote the sin, cos and tan respectively on
 *      [-PI/4, +PI/4]. Reduce the argument x to y1+y2 = x-k*pi/2
 *      in [-pi/4 , +pi/4], and let n = k mod 4.
 *      We have
 *
 *          n        sin(x)      cos(x)        tan(x)
 *     ----------------------------------------------------------
 *          0          S           C             T
 *          1          C          -S            -1/T
 *          2         -S          -C             T
 *          3         -C           S            -1/T
 *     ----------------------------------------------------------
 *
 * Special cases:
 *      Let trig be any of sin, cos, or tan.
 *      trig(+-INF)  is NaN, with signals;
 *      trig(NaN)    is that NaN;
 *
 * Accuracy:
 *      TRIG(x) returns trig(x) nearly rounded
 */

#include "internal/__sf_cos.h"
#include "internal/__sf_rem_pio2.h"
#include "internal/__sf_sin.h"
#include "internal/macros.h"
#include "sf_sin.h"


/* Small multiples of pi/2 rounded to double precision. */
static const sf_float64
s1pio2 = 0x3ff921fb54442d18, /* 1*M_PI_2 */
s2pio2 = 0x400921fb54442d18, /* 2*M_PI_2 */
s3pio2 = 0x4012d97c7f3321d2, /* 3*M_PI_2 */
s4pio2 = 0x401921fb54442d18; /* 4*M_PI_2 */

sf_result32f sf_float32_sin(sf_float32 x, sf_fpu_state fpu) {
    sf_float32 f32tmp;
    sf_float64 y, f64tmp;
    uint32_t ix, sign;
    sf_int32 n;

    sign = x >> 31;
    ix = x & 0x7fffffff;

    if (ix <= 0x3f490fda) {  /* |x| ~<= pi/4 */
        if (ix < 0x39800000) {  /* |x| < 2**-12 */
            /* raise inexact if x!=0 and underflow if subnormal */
            { /* FORCE_EVAL(ix < 0x00800000 ? x/0x1p120f : x+0x1p120f); */
                if (ix < 0x00800000) {
                    SF_F32_OP(fpu, f32tmp, sf_float32_div, x, 0x7b800000);
                } else {
                    SF_F32_OP(fpu, f32tmp, sf_float32_add, x, 0x7b800000);
                }
            }
            return (sf_result32f) { x, fpu };
        }
        { /* return __sindf(x); */
            SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
            return __sf_float32_sin(f64tmp, fpu);
        }
    }
    if (ix <= 0x407b53d1) {  /* |x| ~<= 5*pi/4 */
        if (ix <= 0x4016cbe3) {  /* |x| ~<= 3pi/4 */
            if (sign) {
                { /* return -__cosdf(x + s1pio2); */
                    SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
                    SF_F64_OP(fpu, f64tmp, sf_float64_add, f64tmp, s1pio2);
                    SF_F32_OP(fpu, f32tmp, __sf_float32_cos, f64tmp);
                    return (sf_result32f) { sf_float32_neg(f32tmp), fpu };
                }
            } else {
                { /* return __cosdf(x - s1pio2); */
                    SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
                    SF_F64_OP(fpu, f64tmp, sf_float64_sub, f64tmp, s1pio2);
                    return __sf_float32_cos(f64tmp, fpu);
                }
            }
        }
        { /* return __sindf(sign ? -(x + s2pio2) : -(x - s2pio2)); */
            SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
            if (sign) {
                SF_F64_OP(fpu, f64tmp, sf_float64_add, f64tmp, s2pio2);
            } else {
                SF_F64_OP(fpu, f64tmp, sf_float64_sub, f64tmp, s2pio2);
            }
            return __sf_float32_sin(sf_float64_neg(f64tmp), fpu);
        }
    }
    if (ix <= 0x40e231d5) {  /* |x| ~<= 9*pi/4 */
        if (ix <= 0x40afeddf) {  /* |x| ~<= 7*pi/4 */
            if (sign) {
                /* return __cosdf(x + s3pio2); */
                SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
                SF_F64_OP(fpu, f64tmp, sf_float64_add, f64tmp, s3pio2);
                return __sf_float32_cos(f64tmp, fpu);
            } else {
                /* return -__cosdf(x - s3pio2); */
                SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
                SF_F64_OP(fpu, f64tmp, sf_float64_sub, f64tmp, s3pio2);
                SF_F32_OP(fpu, f32tmp, __sf_float32_cos, f64tmp);
                return (sf_result32f) { sf_float32_neg(f32tmp), fpu };
            }
        }
        { /* return __sindf(sign ? x + s4pio2 : x - s4pio2); */
            SF_F64_OP(fpu, f64tmp, sf_float32_to_float64, x);
            if (sign) {
                SF_F64_OP(fpu, f64tmp, sf_float64_add, f64tmp, s4pio2);
            } else {
                SF_F64_OP(fpu, f64tmp, sf_float64_sub, f64tmp, s4pio2);
            }
            return __sf_float32_sin(f64tmp, fpu);
        }
    }

    /* sin(Inf or NaN) is NaN */
    if (ix >= 0x7f800000)
        return sf_float32_sub(x, x, fpu); /* return x - x; */

    /* general argument reduction needed */
    SF_I32_OP(fpu, n, __sf_float32_rem_pio2, x, &y); /* n = __sf_float32_rem_pio2(x, &y); */
    switch (n & 0x3) {
    case 0:
        return __sf_float32_sin(y, fpu); /* return __sindf(y); */
    case 1:
        return __sf_float32_cos(y, fpu); /* return __cosdf(y); */
    case 2:
        return __sf_float32_sin(sf_float64_neg(y), fpu); /* return __sindf(-y); */
    default:
        { /* return -__cosdf(y); */
            SF_F32_OP(fpu, f32tmp, __sf_float32_cos, y);
            return (sf_result32f) { sf_float32_neg(f32tmp), fpu };
        }
    }
}

sf_result64f sf_float64_sin(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 y[2], tmp;
    uint64_t ix;
    sf_int32 n;

    /* High word of x. */
    ix = x & 0x7fffffff00000000;

    /* |x| ~< pi/4 */
    if (ix <= 0x3fe921fb00000000) {
        if (ix < 0x3e50000000000000) {  /* |x| < 2**-26 */
            /* raise inexact if x != 0 and underflow if subnormal*/
            { /* FORCE_EVAL(ix < 0x00100000 ? x/0x1p120f : x+0x1p120f); */
                if (ix < 0x0010000000000000) {
                    SF_F64_OP(fpu, tmp, sf_float64_div, x, 0x4770000000000000);
                } else {
                    SF_F64_OP(fpu, tmp, sf_float64_add, x, 0x4770000000000000);
                }
            }
            return (sf_result64f) { x, fpu };
        }
        return __sf_float64_sin(x, 0x0, 0, fpu);
    }

    /* sin(Inf or NaN) is NaN */
    if (ix >= 0x7ff0000000000000)
        return sf_float64_sub(x, x, fpu); /* x - x */

    /* argument reduction needed */
    SF_I32_OP(fpu, n, __sf_float64_rem_pio2, x, y); /* n = __sf_float64_rem_pio2(x, y); */
    switch (n & 0x3) {
    case 0:
        return __sf_float64_sin(y[0], y[1], 1, fpu);
    case 1:
        return __sf_float64_cos(y[0], y[1], fpu);
    case 2:
        SF_F64_OP(fpu, tmp, __sf_float64_sin, y[0], y[1], 1);
        return (sf_result64f) { sf_float64_neg(tmp), fpu };
    default:
        SF_F64_OP(fpu, tmp, __sf_float64_cos, y[0], y[1]);
        return (sf_result64f) { sf_float64_neg(tmp), fpu };
    }
}
