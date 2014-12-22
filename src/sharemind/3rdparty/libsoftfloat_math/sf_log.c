/* origin: musl-libc src/math/log.c */
/* origin: musl-libc src/math/logf.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/e_logf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */
/* origin: FreeBSD /usr/src/lib/msun/src/e_log.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/* log(x)
 * Return the logarithm of x
 *
 * Method :
 *   1. Argument Reduction: find k and f such that
 *                      x = 2^k * (1+f),
 *         where  sqrt(2)/2 < 1+f < sqrt(2) .
 *
 *   2. Approximation of log(1+f).
 *      Let s = f/(2+f) ; based on log(1+f) = log(1+s) - log(1-s)
 *               = 2s + 2/3 s**3 + 2/5 s**5 + .....,
 *               = 2s + s*R
 *      We use a special Remez algorithm on [0,0.1716] to generate
 *      a polynomial of degree 14 to approximate R The maximum error
 *      of this polynomial approximation is bounded by 2**-58.45. In
 *      other words,
 *                      2      4      6      8      10      12      14
 *          R(z) ~ Lg1*s +Lg2*s +Lg3*s +Lg4*s +Lg5*s  +Lg6*s  +Lg7*s
 *      (the values of Lg1 to Lg7 are listed in the program)
 *      and
 *          |      2          14          |     -58.45
 *          | Lg1*s +...+Lg7*s    -  R(z) | <= 2
 *          |                             |
 *      Note that 2s = f - s*f = f - hfsq + s*hfsq, where hfsq = f*f/2.
 *      In order to guarantee error in log below 1ulp, we compute log
 *      by
 *              log(1+f) = f - s*(f - R)        (if f is not too large)
 *              log(1+f) = f - (hfsq - s*(hfsq+R)).     (better accuracy)
 *
 *      3. Finally,  log(x) = k*ln2 + log(1+f).
 *                          = k*ln2_hi+(f-(hfsq-(s*(hfsq+R)+k*ln2_lo)))
 *         Here ln2 is split into two floating point number:
 *                      ln2_hi + ln2_lo,
 *         where n*ln2_hi is always exact for |n| < 2000.
 *
 * Special cases:
 *      log(x) is NaN with signal if x < 0 (including -INF) ;
 *      log(+INF) is +INF; log(0) is -INF with signal;
 *      log(NaN) is that NaN with no signal.
 *
 * Accuracy:
 *      according to an error analysis, the error is always less than
 *      1 ulp (unit in the last place).
 *
 * Constants:
 * The hexadecimal values are the intended ones for the following
 * constants. The decimal values may be used, provided that the
 * compiler will convert from decimal to binary accurately enough
 * to produce the hexadecimal values shown.
 */

#include <math.h>
#include <stdint.h>
#include "internal/macros.h"
#include "sf_log.h"


static const sf_float32
f32ln2_hi = 0x3f317180, /* 6.9313812256e-01 */
f32ln2_lo = 0x3717f7d1, /* 9.0580006145e-06 */
/* |(log(1+s)-log(1-s))/s - Lg(s)| < 2**-34.24 (~[-4.95e-11, 4.97e-11]). */
f32Lg1    = 0x3f2aaaaa, /* 0xaaaaaa.0p-24 */
f32Lg2    = 0x3eccce13, /* 0xccce13.0p-25 */
f32Lg3    = 0x3e91e9ee, /* 0x91e9ee.0p-25 */
f32Lg4    = 0x3e789e26; /* 0xf89e26.0p-26 */

static const sf_float64
f64ln2_hi = 0x3fe62e42fee00000, /* 6.93147180369123816490e-01 */
f64ln2_lo = 0x3dea39ef35793c76, /* 1.90821492927058770002e-10 */
f64Lg1    = 0x3fe5555555555593, /* 6.666666666666735130e-01   */
f64Lg2    = 0x3fd999999997fa04, /* 3.999999999940941908e-01   */
f64Lg3    = 0x3fd2492494229359, /* 2.857142874366239149e-01   */
f64Lg4    = 0x3fcc71c51d8e78af, /* 2.222219843214978396e-01   */
f64Lg5    = 0x3fc7466496cb03de, /* 1.818357216161805012e-01   */
f64Lg6    = 0x3fc39a09d078c69f, /* 1.531383769920937332e-01   */
f64Lg7    = 0x3fc2f112df3e5244; /* 1.479819860511658591e-01   */

sf_result32f sf_float32_log(sf_float32 x, sf_fpu_state fpu) {
    sf_float32 hfsq, f, s, z, R, w, t1, t2, dk, tmp1, tmp2;
    uint32_t hx;
    int32_t k;

    hx = x;
    k = 0;
    if (hx < 0x00800000 || hx >> 31) {  /* x < 2**-126  */
        if (hx << 1 == 0) {
            /* return -1/(x*x); */ /* log(+-0)=-inf */
            SF_F32_OP(fpu, tmp1, sf_float32_mul, x, x);
            return sf_float32_div(0xbf800000, tmp1, fpu);
        }
        if (hx >> 31) {
            /* return (x-x)/0.0f; */ /* log(-#) = NaN */
            SF_F32_OP(fpu, tmp1, sf_float32_sub, x, x);
            return sf_float32_div(tmp1, 0x0, fpu);
        }
        /* subnormal number, scale up x */
        k -= 25;
        SF_F32_OP(fpu, x, sf_float32_mul, x, 0x4c000000); /* x *= 0x1p25f; */
        hx = x;
    } else if (hx >= 0x7f800000) {
        return (sf_result32f) { x, fpu };
    } else if (hx == 0x3f800000) {
        return (sf_result32f) { 0x0, fpu };
    }

    /* reduce x into [sqrt(2)/2, sqrt(2)] */
    hx += 0x3f800000 - 0x3f3504f3;
    k += (int32_t)(hx >> 23) - 0x7f;
    hx = (hx & 0x007fffff) + 0x3f3504f3;
    x = hx;

    SF_F32_OP(fpu, f, sf_float32_sub, x, 0x3f800000); /* f = x - 1.0f; */
    { /* s = f/(2.0f + f); */
        SF_F32_OP(fpu, s, sf_float32_add, 0x40000000, f);
        SF_F32_OP(fpu, s, sf_float32_div, f, s);
    }
    SF_F32_OP(fpu, z, sf_float32_mul, s, s); /* z = s*s; */
    SF_F32_OP(fpu, w, sf_float32_mul, z, z); /* w = z*z; */
    { /* t1= w*(Lg2+w*Lg4); */
        SF_F32_OP(fpu, t1, sf_float32_mul, w, f32Lg4);
        SF_F32_OP(fpu, t1, sf_float32_add, f32Lg2, t1);
        SF_F32_OP(fpu, t1, sf_float32_mul, w, t1);
    }
    { /* t2= z*(Lg1+w*Lg3); */
        SF_F32_OP(fpu, t2, sf_float32_mul, w, f32Lg3);
        SF_F32_OP(fpu, t2, sf_float32_add, f32Lg1, t2);
        SF_F32_OP(fpu, t2, sf_float32_mul, z, t2);
    }
    SF_F32_OP(fpu, R, sf_float32_add, t2, t1); /* R = t2 + t1; */
    { /* hfsq = 0.5f*f*f; */
        SF_F32_OP(fpu, hfsq, sf_float32_mul, 0x3f000000, f);
        SF_F32_OP(fpu, hfsq, sf_float32_mul, hfsq, f);
    }
    SF_F32_OP(fpu, dk, sf_int32_to_float32, k); /* dk = k; */
    { /* return s*(hfsq+R) + dk*ln2_lo - hfsq + f + dk*ln2_hi; */
        SF_F32_OP(fpu, tmp1, sf_float32_add, hfsq, R);
        SF_F32_OP(fpu, tmp1, sf_float32_mul, s, tmp1);

        SF_F32_OP(fpu, tmp2, sf_float32_mul, dk, f32ln2_lo);
        SF_F32_OP(fpu, tmp1, sf_float32_add, tmp1, tmp2);

        SF_F32_OP(fpu, tmp1, sf_float32_sub, tmp1, hfsq);

        SF_F32_OP(fpu, tmp1, sf_float32_add, tmp1, f);

        SF_F32_OP(fpu, tmp2, sf_float32_mul, dk, f32ln2_hi);
        return sf_float32_add(tmp1, tmp2, fpu);
    }
}

sf_result64f sf_float64_log(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 hfsq, f, s, z, R, w, t1, t2, dk, tmp1, tmp2;
    uint64_t hx;
    int64_t k;

    hx = x & 0xffffffff00000000;
    k = 0;
    if (hx < 0x0010000000000000 || hx >> 63) {
        if (x << 1 == 0) {
            /* return -1/(x*x); */  /* log(+-0)=-inf */
            SF_F64_OP(fpu, tmp1, sf_float64_mul, x, x);
            return sf_float64_div(0xbff0000000000000, tmp1, fpu);
        }
        if (hx >> 63) {
            /* return (x-x)/0.0; */ /* log(-#) = NaN */
            SF_F64_OP(fpu, tmp1, sf_float64_sub, x, x);
            return sf_float64_div(tmp1, 0x0, fpu);
        }
        /* subnormal number, scale x up */
        k -= 54;
        SF_F64_OP(fpu, x, sf_float64_mul, x, 0x4350000000000000); /* x *= 0x1p54; */
        hx = x & 0xffffffff00000000;
    } else if (hx >= 0x7ff0000000000000) {
        return (sf_result64f) { x, fpu };
    } else if (hx == 0x3ff0000000000000 && x << 32 == 0) {
        return (sf_result64f) { 0x0, fpu };
    }

    /* reduce x into [sqrt(2)/2, sqrt(2)] */
    hx += 0x3ff0000000000000 - 0x3fe6a09e00000000;
    k += (int64_t)(hx >> 52) - 0x3ff;
    hx = (hx & 0x000fffff00000000) + 0x3fe6a09e00000000;
    x = (sf_float64)hx | (x & 0xffffffff);

    SF_F64_OP(fpu, f, sf_float64_sub, x, 0x3ff0000000000000); /* f = x - 1.0; */
    { /* hfsq = 0.5*f*f; */
        SF_F64_OP(fpu, hfsq, sf_float64_mul, 0x3fe0000000000000, f);
        SF_F64_OP(fpu, hfsq, sf_float64_mul, hfsq, f);
    }
    { /* s = f/(2.0+f); */
        SF_F64_OP(fpu, s, sf_float64_add, 0x4000000000000000, f);
        SF_F64_OP(fpu, s, sf_float64_div, f, s);
    }
    SF_F64_OP(fpu, z, sf_float64_mul, s, s); /* z = s*s; */
    SF_F64_OP(fpu, w, sf_float64_mul, z, z); /* w = z*z; */
    { /* t1 = w*(f64Lg2+w*(f64Lg4+w*f64Lg6)); */
        SF_F64_OP(fpu, t1, sf_float64_mul, w, f64Lg6);
        SF_F64_OP(fpu, t1, sf_float64_add, f64Lg4, t1);
        SF_F64_OP(fpu, t1, sf_float64_mul, w, t1);
        SF_F64_OP(fpu, t1, sf_float64_add, f64Lg2, t1);
        SF_F64_OP(fpu, t1, sf_float64_mul, w, t1);
    }
    { /* t2 = z*(f64Lg1+w*(f64Lg3+w*(f64Lg5+w*f64Lg7))); */
        SF_F64_OP(fpu, t2, sf_float64_mul, w, f64Lg7);
        SF_F64_OP(fpu, t2, sf_float64_add, f64Lg5, t2);
        SF_F64_OP(fpu, t2, sf_float64_mul, w, t2);
        SF_F64_OP(fpu, t2, sf_float64_add, f64Lg3, t2);
        SF_F64_OP(fpu, t2, sf_float64_mul, w, t2);
        SF_F64_OP(fpu, t2, sf_float64_add, f64Lg1, t2);
        SF_F64_OP(fpu, t2, sf_float64_mul, z, t2);
    }
    SF_F64_OP(fpu, R, sf_float64_add, t2, t1); /* R = t2 + t1; */
    SF_F64_OP(fpu, dk, sf_int64_to_float64, k); /* dk = k; */
    { /* return s*(hfsq+R) + dk*f64ln2_lo - hfsq + f + dk*f64ln2_hi; */
        SF_F64_OP(fpu, tmp1, sf_float64_add, hfsq, R);
        SF_F64_OP(fpu, tmp1, sf_float64_mul, s, tmp1);

        SF_F64_OP(fpu, tmp2, sf_float64_mul, dk, f64ln2_lo);
        SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, tmp2);

        SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, hfsq);

        SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, f);

        SF_F64_OP(fpu, tmp2, sf_float64_mul, dk, f64ln2_hi);
        return sf_float64_add(tmp1, tmp2, fpu);
    }
}
