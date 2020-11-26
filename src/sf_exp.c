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

/* origin: musl-libc src/math/exp.c */
/* origin: musl-libc src/math/expf.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/e_expf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */
/* origin: FreeBSD /usr/src/lib/msun/src/e_exp.c */
/*
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/* exp(x)
 * Returns the exponential of x.
 *
 * Method
 *   1. Argument reduction:
 *      Reduce x to an r so that |r| <= 0.5*ln2 ~ 0.34658.
 *      Given x, find r and integer k such that
 *
 *               x = k*ln2 + r,  |r| <= 0.5*ln2.
 *
 *      Here r will be represented as r = hi-lo for better
 *      accuracy.
 *
 *   2. Approximation of exp(r) by a special rational function on
 *      the interval [0,0.34658]:
 *      Write
 *          R(r**2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r**4/360 + ...
 *      We use a special Remez algorithm on [0,0.34658] to generate
 *      a polynomial of degree 5 to approximate R. The maximum error
 *      of this polynomial approximation is bounded by 2**-59. In
 *      other words,
 *          R(z) ~ 2.0 + f64P1*z + f64P2*z**2 + f64P3*z**3 + f64P4*z**4 + f64P5*z**5
 *      (where z=r*r, and the values of f64P1 to f64P5 are listed below)
 *      and
 *          |                  5          |     -59
 *          | 2.0+f64P1*z+...+f64P5*z   -  R(z) | <= 2
 *          |                             |
 *      The computation of exp(r) thus becomes
 *                              2*r
 *              exp(r) = 1 + ----------
 *                            R(r) - r
 *                                 r*c(r)
 *                     = 1 + r + ----------- (for better accuracy)
 *                                2 - c(r)
 *      where
 *                              2       4             10
 *              c(r) = r - (f64P1*r  + f64P2*r  + ... + f64P5*r   ).
 *
 *   3. Scale back to obtain exp(x):
 *      From step 1, we have
 *         exp(x) = 2^k * exp(r)
 *
 * Special cases:
 *      exp(INF) is INF, exp(NaN) is NaN;
 *      exp(-INF) is 0, and
 *      for finite argument, only exp(0)=1 is exact.
 *
 * Accuracy:
 *      according to an error analysis, the error is always less than
 *      1 ulp (unit in the last place).
 *
 * Misc. info.
 *      For IEEE double
 *          if x >  709.782712893383973096 then exp(x) overflows
 *          if x < -745.133219101941108420 then exp(x) underflows
 */

#include "sf_exp.h"

#include <stdint.h>
#include "internal/macros.h"
#include "sf_scalbn.h"


static const sf_float32
f32half[2] = { 0x3f000000,
            0xbf000000 },
f32ln2hi   =   0x3f317200, /* 6.9314575195e-1f */
f32ln2lo   =   0x35bfbe8e, /* 1.4286067653e-6f */
f32invln2  =   0x3fb8aa3b, /* 1.4426950216e+0f */
/*
 * Domain [-0.34568, 0.34568], range ~[-4.278e-9, 4.447e-9]:
 * |x*(exp(x)+1)/(exp(x)-1) - p(x)| < 2**-27.74
 */
f32P1      =   0x3e2aaa8f, /*  1.6666625440e-1f */
f32P2      =   0xbb355215; /* -2.7667332906e-3f */

static const sf_float64
f64half[2] = { 0x3fe0000000000000,   /*  0.5                        */
               0xbfe0000000000000 }, /* -0.5                        */
f64ln2hi   =   0x3fe62e42fee00000,   /*  6.93147180369123816490e-01 */
f64ln2lo   =   0x3dea39ef35793c76,   /*  1.90821492927058770002e-10 */
f64invln2  =   0x3ff71547652b82fe,   /*  1.44269504088896338700e+00 */
f64P1      =   0x3fc555555555553e,   /*  1.66666666666666019037e-01 */
f64P2      =   0xbf66c16c16bebd93,   /* -2.77777777770155933842e-03 */
f64P3      =   0x3f11566aaf25de2c,   /*  6.61375632143793436117e-05 */
f64P4      =   0xbebbbd41c5d26bf1,   /* -1.65339022054652515390e-06 */
f64P5      =   0x3e66376972bea4d0;   /*  4.13813679705723846039e-08 */

sf_result32f sf_float32_exp(sf_float32 x, sf_fpu_state fpu) {
    sf_float32 hi, lo, c, xx, y, tmp;
    sf_int32 k;
    uint32_t hx, sign;

    sign = x >> 31;   /* sign bit of x */
    hx = x & 0x7fffffff;  /* high word of |x| */

    /* special cases */
    if (hx >= 0x42aeac50) {  /* if |x| >= -87.33655f or NaN */
        if (hx >= 0x42b17218 && !sign) {  /* x >= 88.722839f */
            /* overflow */
            SF_F32_OP(fpu, x, sf_float32_mul, x, 0x7f000000); /* x *= 0x1p127f; */
            return (sf_result32f) { x, fpu };
        }
        if (sign) {
            /* underflow */
            SF_F32_OP(fpu, tmp, sf_float32_div, 0x80000001, x); /* FORCE_EVAL(-0x1p-149f/x); */
            if (hx >= 0x42cff1b5)  /* x <= -103.972084f */
                return (sf_result32f) { 0x0, fpu };
        }
    }

    /* argument reduction */
    if (hx > 0x3eb17218) {  /* if |x| > 0.5 ln2 */
        if (hx > 0x3f851592) {  /* if |x| > 1.5 ln2 */
            /* k = f32invln2*x + f32half[sign]; */
            SF_F32_OP(fpu, tmp, sf_float32_mul, f32invln2, x);
            SF_F32_OP(fpu, tmp, sf_float32_add, tmp, f32half[sign]);
            SF_I32_OP(fpu, k, sf_float32_to_int32_round_to_zero, tmp);
        } else {
            k = sign ? -1 : 1; /* k = 1 - sign - sign; */
        }
        { /* hi = x - k * f32ln2hi; */ /* k*f32ln2hi is exact here */
            SF_F32_OP(fpu, hi, sf_int32_to_float32, k);
            SF_F32_OP(fpu, hi, sf_float32_mul, hi, f32ln2hi);
            SF_F32_OP(fpu, hi, sf_float32_sub, x, hi);
        }
        { /* lo = k * f32ln2lo; */
            SF_F32_OP(fpu, lo, sf_int32_to_float32, k);
            SF_F32_OP(fpu, lo, sf_float32_mul, lo, f32ln2lo);
        }
        SF_F32_OP(fpu, x, sf_float32_sub, hi, lo); /* x = hi - lo; */
    } else if (hx > 0x39000000) {  /* |x| > 2**-14 */
        k = 0;
        hi = x;
        lo = 0x0;
    } else {
        /* raise inexact */
        SF_F32_OP(fpu, tmp, sf_float32_add, 0x7f000000, x); /* FORCE_EVAL(0x1p127f + x); */
        return sf_float32_add(0x3f800000, x, fpu);
    }

    /* x is now in primary range */
    SF_F32_OP(fpu, xx, sf_float32_mul, x, x); /* xx = x*x; */
    { /* c = x - xx*(P1+xx*P2); */
        SF_F32_OP(fpu, c, sf_float32_mul, xx, f32P2);
        SF_F32_OP(fpu, c, sf_float32_add, f32P1, c);
        SF_F32_OP(fpu, c, sf_float32_mul, xx, c);
        SF_F32_OP(fpu, c, sf_float32_sub, x, c);
    }
    { /* y = 1 + (x*c/(2-c) - lo + hi); */
        SF_F32_OP(fpu, y, sf_float32_mul, x, c);
        SF_F32_OP(fpu, tmp, sf_float32_sub, 0x40000000, c);
        SF_F32_OP(fpu, y, sf_float32_div, y, tmp);
        SF_F32_OP(fpu, y, sf_float32_sub, y, lo);
        SF_F32_OP(fpu, y, sf_float32_add, y, hi);
        SF_F32_OP(fpu, y, sf_float32_add, 0x3f800000, y);
    }
    if (k == 0)
        return (sf_result32f) { y, fpu };

    return sf_float32_scalbn(y, k, fpu);
}

sf_result64f sf_float64_exp(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 hi, lo, c, xx, y, tmp;
    sf_flag flag;
    sf_int64 k;
    uint64_t hx, sign;

    sign = x >> 63;
    hx = x & 0x7fffffff00000000;  /* high word of |x| */

    /* special cases */
    if (hx >= 0x4086232b00000000) {  /* if |x| >= 708.39... */
        if (sf_float64_is_signaling_nan(x))
            return (sf_result64f) { x, fpu };
        SF_FLAG_OP(fpu, flag, sf_float64_lt, 0x40862e42fefa39ef, x);
        if (flag) { /* x > 709.782712893383973096 */
            /* overflow if x!=inf */
            SF_F64_OP(fpu, x, sf_float64_mul, x, 0x7fe0000000000000); /* x *= 0x1p1023; */
            return (sf_result64f) { x, fpu };
        }
        SF_FLAG_OP(fpu, flag, sf_float64_lt, x, 0xc086232bdd7abcd2);
        if (flag) { /* x < -708.39641853226410622 */
            /* underflow if x!=-inf */
            SF_F64_OP(fpu, tmp, sf_float64_div, 0xb6a0000000000000, x); /* FORCE_EVAL((float)(-0x1p-149/x)); */
            SF_FLAG_OP(fpu, flag, sf_float64_lt, x, 0xc0874910d52d3051);
            if (flag) /* x < -745.13321910194110842 */
                return (sf_result64f) { 0x0, fpu };
        }
    }

    /* argument reduction */
    if (hx > 0x3fd62e4200000000) {  /* if |x| > 0.5 ln2 */
        if (hx >= 0x3ff0a2b200000000) {  /* if |x| >= 1.5 ln2 */
            /* k = (int)(f64invln2*x + f64half[sign]); */
            SF_F64_OP(fpu, tmp, sf_float64_mul, f64invln2, x);
            SF_F64_OP(fpu, tmp, sf_float64_add, tmp, f64half[sign]);
            SF_I64_OP(fpu, k, sf_float64_to_int64_round_to_zero, tmp);
        } else {
            k = sign ? -1 : 1; /* 1 - sign - sign; */
        }

        { /* hi = x - k*f64ln2hi; */  /* k*f64ln2hi is exact here */
            SF_F64_OP(fpu, hi, sf_int64_to_float64, k);
            SF_F64_OP(fpu, hi, sf_float64_mul, hi, f64ln2hi);
            SF_F64_OP(fpu, hi, sf_float64_sub, x, hi);
        }
        { /* lo = k*f64ln2lo; */
            SF_F64_OP(fpu, lo, sf_int64_to_float64, k);
            SF_F64_OP(fpu, lo, sf_float64_mul, lo, f64ln2lo);
        }
        SF_F64_OP(fpu, x, sf_float64_sub, hi, lo); /* x = hi - lo; */
    } else if (hx > 0x3e30000000000000)  {  /* if |x| > 2**-28 */
        k = 0;
        hi = x;
        lo = 0;
    } else {
        /* inexact if x!=0 */
        SF_F64_OP(fpu, tmp, sf_float64_add, 0x7fe0000000000000, x); /* FORCE_EVAL(0x1p1023 + x); */
        return sf_float64_add(0x3ff0000000000000, x, fpu); /* return 1 + x; */
    }

    /* x is now in primary range */
    SF_F64_OP(fpu, xx, sf_float64_mul, x, x); /* xx = x*x; */

    { /* c = x - xx*(f64P1+xx*(f64P2+xx*(f64P3+xx*(f64P4+xx*f64P5)))); */
        SF_F64_OP(fpu, c, sf_float64_mul, xx, f64P5);
        SF_F64_OP(fpu, c, sf_float64_add, f64P4, c);
        SF_F64_OP(fpu, c, sf_float64_mul, xx, c);
        SF_F64_OP(fpu, c, sf_float64_add, f64P3, c);
        SF_F64_OP(fpu, c, sf_float64_mul, xx, c);
        SF_F64_OP(fpu, c, sf_float64_add, f64P2, c);
        SF_F64_OP(fpu, c, sf_float64_mul, xx, c);
        SF_F64_OP(fpu, c, sf_float64_add, f64P1, c);
        SF_F64_OP(fpu, c, sf_float64_mul, xx, c);
        SF_F64_OP(fpu, c, sf_float64_sub, x, c);
    }
    { /* y = 1 + (x*c/(2-c) - lo + hi); */
        SF_F64_OP(fpu, y, sf_float64_mul, x, c);
        SF_F64_OP(fpu, tmp, sf_float64_sub, 0x4000000000000000, c);
        SF_F64_OP(fpu, y, sf_float64_div, y, tmp);
        SF_F64_OP(fpu, y, sf_float64_sub, y, lo);
        SF_F64_OP(fpu, y, sf_float64_add, y, hi);
        SF_F64_OP(fpu, y, sf_float64_add, 0x3ff0000000000000, y);
    }
    if (k == 0)
        return (sf_result64f) { y, fpu };

    return sf_float64_scalbn(y, k, fpu);
}
