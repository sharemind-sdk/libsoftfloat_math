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

/* origin: musl-libc src/math/__sin.c */
/* origin: musl-libc src/math/__sindf.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/k_sin.c */
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
/* __sin( x, y, iy)
 * kernel sin function on ~[-pi/4, pi/4] (except on -0), pi/4 ~ 0.7854
 * Input x is assumed to be bounded by ~pi/4 in magnitude.
 * Input y is the tail of x.
 * Input iy indicates whether y is 0. (if iy=0, y assume to be 0).
 *
 * Algorithm
 *      1. Since sin(-x) = -sin(x), we need only to consider positive x.
 *      2. Callers must return sin(-0) = -0 without calling here since our
 *         odd polynomial is not evaluated in a way that preserves -0.
 *         Callers may do the optimization sin(x) ~ x for tiny x.
 *      3. sin(x) is approximated by a polynomial of degree 13 on
 *         [0,pi/4]
 *                               3            13
 *              sin(x) ~ x + S1*x + ... + S6*x
 *         where
 *
 *      |sin(x)         2     4     6     8     10     12  |     -58
 *      |----- - (1+S1*x +S2*x +S3*x +S4*x +S5*x  +S6*x   )| <= 2
 *      |  x                                               |
 *
 *      4. sin(x+y) = sin(x) + sin'(x')*y
 *                  ~ sin(x) + (1-x*x/2)*y
 *         For better accuracy, let
 *                   3      2      2      2      2
 *              r = x *(S2+x *(S3+x *(S4+x *(S5+x *S6))))
 *         then                   3    2
 *              sin(x) = x + (S1*x + (x *(r-y/2)+y))
 */

#include "__sf_sin.h"
#include "macros.h"


/* |sin(x)/x - s(x)| < 2**-37.5 (~[-4.89e-12, 4.824e-12]). */
static const sf_float64
f32S1 = 0xbfc5555554cbac77, /* -0.166666666416265235595 */
f32S2 = 0x3f811110896efbb2, /*  0.0083333293858894631756 */
f32S3 = 0xbf2a00f9e2cae774, /* -0.000198393348360966317347 */
f32S4 = 0x3ec6cd878c3b46a7; /*  0.0000027183114939898219064 */

static const sf_float64
f64S1 = 0xbfc5555555555549, /* -1.66666666666666324348e-01 */
f64S2 = 0x3f8111111110f8a6, /*  8.33333333332248946124e-03 */
f64S3 = 0xbf2a01a019c161d5, /* -1.98412698298579493134e-04 */
f64S4 = 0x3ec71de357b1fe7d, /*  2.75573137070700676789e-06 */
f64S5 = 0xbe5ae5e68a2b9ceb, /* -2.50507602534068634195e-08 */
f64S6 = 0x3de5d93a5acfd57c; /*  1.58969099521155010221e-10 */

sf_result32f __sf_float32_sin(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 r, s, w, z, tmp1, tmp2;

    /* Try to optimize for parallel evaluation as in __tandf.c. */
    SF_F64_OP(fpu, z, sf_float64_mul, x, x); /* z = x*x; */
    SF_F64_OP(fpu, w, sf_float64_mul, z, z); /* w = z*z; */
    { /* r = S3 + z*S4; */
        SF_F64_OP(fpu, r, sf_float64_mul, z, f32S4);
        SF_F64_OP(fpu, r, sf_float64_add, f32S3, r);
    }
    SF_F64_OP(fpu, s, sf_float64_mul, z, x); /* s = z*x; */
    { /* return (x + s*(S1 + z*S2)) + s*w*r; */
        SF_F64_OP(fpu, tmp1, sf_float64_mul, z, f32S2);
        SF_F64_OP(fpu, tmp1, sf_float64_add, f32S1, tmp1);
        SF_F64_OP(fpu, tmp1, sf_float64_mul, s, tmp1);
        SF_F64_OP(fpu, tmp1, sf_float64_add, x, tmp1);

        SF_F64_OP(fpu, tmp2, sf_float64_mul, s, w);
        SF_F64_OP(fpu, tmp2, sf_float64_mul, tmp2, r);

        SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, tmp2);
        return sf_float64_to_float32(tmp1, fpu);
    }
}

sf_result64f __sf_float64_sin(sf_float64 x, sf_float64 y, int64_t iy, sf_fpu_state fpu) {
    sf_float64 z, r, v, w, tmp1, tmp2;

    SF_F64_OP(fpu, z, sf_float64_mul, x, x); /* z = x*x; */
    SF_F64_OP(fpu, w, sf_float64_mul, z, z); /* w = z*z; */
    { /* r = S2 + z*(S3 + z*S4) + z*w*(S5 + z*S6); */
        SF_F64_OP(fpu, r, sf_float64_mul, z, f64S4);
        SF_F64_OP(fpu, r, sf_float64_add, f64S3, r);
        SF_F64_OP(fpu, r, sf_float64_mul, z, r);

        SF_F64_OP(fpu, tmp1, sf_float64_mul, z, f64S6);
        SF_F64_OP(fpu, tmp1, sf_float64_add, f64S5, tmp1);
        SF_F64_OP(fpu, tmp2, sf_float64_mul, z, w);
        SF_F64_OP(fpu, tmp1, sf_float64_mul, tmp2, tmp1);

        SF_F64_OP(fpu, r, sf_float64_add, f64S2, r);
        SF_F64_OP(fpu, r, sf_float64_add, r, tmp1);
    }
    SF_F64_OP(fpu, v, sf_float64_mul, z, x); /* v = z*x; */

    if (iy == 0) {
        { /* return x + v*(f64S1 + z*r); */
            SF_F64_OP(fpu, tmp1, sf_float64_mul, z, r);
            SF_F64_OP(fpu, tmp1, sf_float64_add, f64S1, tmp1);
            SF_F64_OP(fpu, tmp1, sf_float64_mul, v, tmp1);
            return sf_float64_add(x, tmp1, fpu);
        }
    } else {
        { /* return x - ((z*(0.5*y - v*r) - y) - v*S1); */
            SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x3fe0000000000000, y);
            SF_F64_OP(fpu, tmp2, sf_float64_mul, v, r);
            SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, tmp2);
            SF_F64_OP(fpu, tmp1, sf_float64_mul, z, tmp1);
            SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, y);
            SF_F64_OP(fpu, tmp2, sf_float64_mul, v, f64S1);
            SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, tmp2);
            return sf_float64_sub(x, tmp1, fpu);
        }
    }
}
