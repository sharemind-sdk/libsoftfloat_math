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

/* origin: musl-libc src/math/__cos.c */
/* origin: musl-libc src/math/__cosdf.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/k_cosf.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 * Debugged and optimized by Bruce D. Evans.
 */
/* origin: FreeBSD /usr/src/lib/msun/src/k_cos.c */
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
/*
 * __cos( x,  y )
 * kernel cos function on [-pi/4, pi/4], pi/4 ~ 0.785398164
 * Input x is assumed to be bounded by ~pi/4 in magnitude.
 * Input y is the tail of x.
 *
 * Algorithm
 *      1. Since cos(-x) = cos(x), we need only to consider positive x.
 *      2. if x < 2^-27 (hx<0x3e400000 0), return 1 with inexact if x!=0.
 *      3. cos(x) is approximated by a polynomial of degree 14 on
 *         [0,pi/4]
 *                                       4            14
 *              cos(x) ~ 1 - x*x/2 + C1*x + ... + C6*x
 *         where the remez error is
 *
 *      |              2     4     6     8     10    12     14 |     -58
 *      |cos(x)-(1-.5*x +C1*x +C2*x +C3*x +C4*x +C5*x  +C6*x  )| <= 2
 *      |                                                      |
 *
 *                     4     6     8     10    12     14
 *      4. let r = C1*x +C2*x +C3*x +C4*x +C5*x  +C6*x  , then
 *             cos(x) ~ 1 - x*x/2 + r
 *         since cos(x+y) ~ cos(x) - sin(x)*y
 *                        ~ cos(x) - x*y,
 *         a correction term is necessary in cos(x) and hence
 *              cos(x+y) = 1 - (x*x/2 - (r - x*y))
 *         For better accuracy, rearrange to
 *              cos(x+y) ~ w + (tmp + (r-x*y))
 *         where w = 1 - x*x/2 and tmp is a tiny correction term
 *         (1 - x*x/2 == w + tmp exactly in infinite precision).
 *         The exactness of w + tmp in infinite precision depends on w
 *         and tmp having the same precision as x.  If they have extra
 *         precision due to compiler bugs, then the extra precision is
 *         only good provided it is retained in all terms of the final
 *         expression for cos().  Retention happens in all cases tested
 *         under FreeBSD, so don't pessimize things by forcibly clipping
 *         any extra precision in w.
 */

#include "__sf_cos.h"
#include "macros.h"


/* |cos(x) - c(x)| < 2**-34.1 (~[-5.37e-11, 5.295e-11]). */
static const sf_float64
f32C0 = 0xbfdffffffd0c5e81, /* -0.499999997251031003120 */
f32C1 = 0x3fa55553e1053a42, /*  0.0416666233237390631894 */
f32C2 = 0xbf56c087e80f1e27, /* -0.00138867637746099294692 */
f32C3 = 0x3ef99342e0ee5069; /*  0.0000243904487962774090654 */

static const sf_float64
f64C1 = 0x3fa555555555554c, /*  4.16666666666666019037e-02 */
f64C2 = 0xbf56c16c16c15177, /* -1.38888888888741095749e-03 */
f64C3 = 0x3efa01a019cb1590, /*  2.48015872894767294178e-05 */
f64C4 = 0xbe927e4f809c52ad, /* -2.75573143513906633035e-07 */
f64C5 = 0x3e21ee9ebdb4b1c4, /*  2.08757232129817482790e-09 */
f64C6 = 0xbda8fae9be8838d4; /* -1.13596475577881948265e-11 */

sf_result32f __sf_float32_cos(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 r, w, z, tmp1, tmp2;

    /* Try to optimize for parallel evaluation as in __tandf.c. */
    SF_F64_OP(fpu, z, sf_float64_mul, x, x); /* z = x*x; */
    SF_F64_OP(fpu, w, sf_float64_mul, z, z); /* w = z*z; */
    { /* r = C2+z*C3; */
        SF_F64_OP(fpu, r, sf_float64_mul, z, f32C3);
        SF_F64_OP(fpu, r, sf_float64_add, f32C2, r);
    }
    { /* return ((1.0+z*C0) + w*C1) + (w*z)*r; */
        SF_F64_OP(fpu, tmp1, sf_float64_mul, z, f32C0);
        SF_F64_OP(fpu, tmp1, sf_float64_add, 0x3ff0000000000000, tmp1);

        SF_F64_OP(fpu, tmp2, sf_float64_mul, w, f32C1);
        SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, tmp2);

        SF_F64_OP(fpu, tmp2, sf_float64_mul, w, z);
        SF_F64_OP(fpu, tmp2, sf_float64_mul, tmp2, r);
        SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, tmp2);

        return sf_float64_to_float32(tmp1, fpu);
    }
}

sf_result64f __sf_float64_cos(sf_float64 x, sf_float64 y, sf_fpu_state fpu) {
    sf_float64 hz, z, r, w, tmp1, tmp2;

    SF_F64_OP(fpu, z, sf_float64_mul, x, x); /* z  = x*x; */
    SF_F64_OP(fpu, w, sf_float64_mul, z, z); /* w  = z*z; */
    { /* r  = z*(C1+z*(C2+z*C3)) + w*w*(C4+z*(C5+z*C6)); */
        SF_F64_OP(fpu, r, sf_float64_mul, z, f64C3);
        SF_F64_OP(fpu, r, sf_float64_add, f64C2, r);
        SF_F64_OP(fpu, r, sf_float64_mul, z, r);
        SF_F64_OP(fpu, r, sf_float64_add, f64C1, r);
        SF_F64_OP(fpu, r, sf_float64_mul, z, r);

        SF_F64_OP(fpu, tmp1, sf_float64_mul, z, f64C6);
        SF_F64_OP(fpu, tmp1, sf_float64_add, f64C5, tmp1);
        SF_F64_OP(fpu, tmp1, sf_float64_mul, z, tmp1);
        SF_F64_OP(fpu, tmp1, sf_float64_add, f64C4, tmp1);
        SF_F64_OP(fpu, tmp2, sf_float64_mul, w, w);
        SF_F64_OP(fpu, tmp1, sf_float64_mul, tmp2, tmp1);

        SF_F64_OP(fpu, r, sf_float64_add, r, tmp1);
    }
    SF_F64_OP(fpu, hz, sf_float64_mul, 0x3fe0000000000000, z); /* hz = 0.5*z; */
    SF_F64_OP(fpu, w, sf_float64_sub, 0x3ff0000000000000, hz); /* w  = 1.0-hz; */
    { /* return w + (((1.0-w)-hz) + (z*r-x*y)); */
        SF_F64_OP(fpu, tmp1, sf_float64_mul, z, r);
        SF_F64_OP(fpu, tmp2, sf_float64_mul, x, y);
        SF_F64_OP(fpu, tmp2, sf_float64_sub, tmp1, tmp2);

        SF_F64_OP(fpu, tmp1, sf_float64_sub, 0x3ff0000000000000, w);
        SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, hz);

        SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, tmp2);
        return sf_float64_add(w, tmp1, fpu);
    }
}
