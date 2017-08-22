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

/* origin: musl-libc src/math/erf.c */
/* origin: musl-libc src/math/erff.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/s_erff.c */
/*
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */
/* origin: FreeBSD /usr/src/lib/msun/src/s_erf.c */
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
/* double erf(double x)
 * double erfc(double x)
 *                           x
 *                    2      |\
 *     erf(x)  =  ---------  | exp(-t*t)dt
 *                 sqrt(pi) \|
 *                           0
 *
 *     erfc(x) =  1-erf(x)
 *  Note that
 *              erf(-x) = -erf(x)
 *              erfc(-x) = 2 - erfc(x)
 *
 * Method:
 *      1. For |x| in [0, 0.84375]
 *          erf(x)  = x + x*R(x^2)
 *          erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
 *                  = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
 *         where R = P/Q where P is an odd poly of degree 8 and
 *         Q is an odd poly of degree 10.
 *                                               -57.90
 *                      | R - (erf(x)-x)/x | <= 2
 *
 *
 *         Remark. The formula is derived by noting
 *          erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
 *         and that
 *          2/sqrt(pi) = 1.128379167095512573896158903121545171688
 *         is close to one. The interval is chosen because the fix
 *         point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
 *         near 0.6174), and by some experiment, 0.84375 is chosen to
 *         guarantee the error is less than one ulp for erf.
 *
 *      2. For |x| in [0.84375,1.25], let s = |x| - 1, and
 *         c = 0.84506291151 rounded to single (24 bits)
 *              erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
 *              erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
 *                        1+(c+P1(s)/Q1(s))    if x < 0
 *              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
 *         Remark: here we use the taylor series expansion at x=1.
 *              erf(1+s) = erf(1) + s*Poly(s)
 *                       = 0.845.. + P1(s)/Q1(s)
 *         That is, we use rational approximation to approximate
 *                      erf(1+s) - (c = (single)0.84506291151)
 *         Note that |P1/Q1|< 0.078 for x in [0.84375,1.25]
 *         where
 *              P1(s) = degree 6 poly in s
 *              Q1(s) = degree 6 poly in s
 *
 *      3. For x in [1.25,1/0.35(~2.857143)],
 *              erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
 *              erf(x)  = 1 - erfc(x)
 *         where
 *              R1(z) = degree 7 poly in z, (z=1/x^2)
 *              S1(z) = degree 8 poly in z
 *
 *      4. For x in [1/0.35,28]
 *              erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
 *                      = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
 *                      = 2.0 - tiny            (if x <= -6)
 *              erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
 *              erf(x)  = sign(x)*(1.0 - tiny)
 *         where
 *              R2(z) = degree 6 poly in z, (z=1/x^2)
 *              S2(z) = degree 7 poly in z
 *
 *      Note1:
 *         To compute exp(-x*x-0.5625+R/S), let s be a single
 *         precision number and s := x; then
 *              -x*x = -s*s + (s-x)*(s+x)
 *              exp(-x*x-0.5626+R/S) =
 *                      exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
 *      Note2:
 *         Here 4 and 5 make use of the asymptotic series
 *                        exp(-x*x)
 *              erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
 *                        x*sqrt(pi)
 *         We use rational approximation to approximate
 *              g(s)=f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
 *         Here is the error bound for R1/S1 and R2/S2
 *              |R1/S1 - f(x)|  < 2**(-62.57)
 *              |R2/S2 - f(x)|  < 2**(-61.52)
 *
 *      5. For inf > x >= 28
 *              erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
 *              erfc(x) = tiny*tiny (raise underflow) if x > 0
 *                      = 2 - tiny if x<0
 *
 *      7. Special case:
 *              erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
 *              erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2,
 *              erfc/erf(NaN) is NaN
 */

#include <stdint.h>
#include "internal/macros.h"
#include "sf_abs.h"
#include "sf_erf.h"
#include "sf_exp.h"


static const sf_float32
f32erx  = 0x3f58560b, /*  8.4506291151e-01 */
/*
 * Coefficients for approximation to  erf on [0,0.84375]
 */
f32efx8 = 0x3f8375d4, /*  1.0270333290e+00 */
f32pp0  = 0x3e0375d4, /*  1.2837916613e-01 */
f32pp1  = 0xbea66beb, /* -3.2504209876e-01 */
f32pp2  = 0xbce9528f, /* -2.8481749818e-02 */
f32pp3  = 0xbbbd1489, /* -5.7702702470e-03 */
f32pp4  = 0xb7c756b1, /* -2.3763017452e-05 */
f32qq1  = 0x3ecbbbce, /*  3.9791721106e-01 */
f32qq2  = 0x3d852a63, /*  6.5022252500e-02 */
f32qq3  = 0x3ba68116, /*  5.0813062117e-03 */
f32qq4  = 0x390aee49, /*  1.3249473704e-04 */
f32qq5  = 0xb684e21a, /* -3.9602282413e-06 */
/*
 * Coefficients for approximation to  erf  in [0.84375,1.25]
 */
f32pa0  = 0xbb1acdc6, /* -2.3621185683e-03 */
f32pa1  = 0x3ed46805, /*  4.1485610604e-01 */
f32pa2  = 0xbebe9208, /* -3.7220788002e-01 */
f32pa3  = 0x3ea2fe54, /*  3.1834661961e-01 */
f32pa4  = 0xbde31cc2, /* -1.1089469492e-01 */
f32pa5  = 0x3d1151b3, /*  3.5478305072e-02 */
f32pa6  = 0xbb0df9c0, /* -2.1663755178e-03 */
f32qa1  = 0x3dd9f331, /*  1.0642088205e-01 */
f32qa2  = 0x3f0a5785, /*  5.4039794207e-01 */
f32qa3  = 0x3d931ae7, /*  7.1828655899e-02 */
f32qa4  = 0x3e013307, /*  1.2617121637e-01 */
f32qa5  = 0x3c5f6e13, /*  1.3637083583e-02 */
f32qa6  = 0x3c445aa3, /*  1.1984500103e-02 */
/*
 * Coefficients for approximation to  erfc in [1.25,1/0.35]
 */
f32ra0  = 0xbc21a093, /* -9.8649440333e-03 */
f32ra1  = 0xbf31a0b7, /* -6.9385856390e-01 */
f32ra2  = 0xc128f022, /* -1.0558626175e+01 */
f32ra3  = 0xc2798057, /* -6.2375331879e+01 */
f32ra4  = 0xc322658c, /* -1.6239666748e+02 */
f32ra5  = 0xc3389ae7, /* -1.8460508728e+02 */
f32ra6  = 0xc2a2932b, /* -8.1287437439e+01 */
f32ra7  = 0xc11d077e, /* -9.8143291473e+00 */
f32sa1  = 0x419d35ce, /*  1.9651271820e+01 */
f32sa2  = 0x4309a863, /*  1.3765776062e+02 */
f32sa3  = 0x43d9486f, /*  4.3456588745e+02 */
f32sa4  = 0x442158c9, /*  6.4538726807e+02 */
f32sa5  = 0x43d6810b, /*  4.2900814819e+02 */
f32sa6  = 0x42d9451f, /*  1.0863500214e+02 */
f32sa7  = 0x40d23f7c, /*  6.5702495575e+00 */
f32sa8  = 0xbd777f97, /* -6.0424413532e-02 */
/*
 * Coefficients for approximation to  erfc in [1/.35,28]
 */
f32rb0  = 0xbc21a092, /* -9.8649431020e-03 */
f32rb1  = 0xbf4c9dd4, /* -7.9928326607e-01 */
f32rb2  = 0xc18e104b, /* -1.7757955551e+01 */
f32rb3  = 0xc320a2ea, /* -1.6063638306e+02 */
f32rb4  = 0xc41f6441, /* -6.3756646729e+02 */
f32rb5  = 0xc480230b, /* -1.0250950928e+03 */
f32rb6  = 0xc3f1c275, /* -4.8351919556e+02 */
f32sb1  = 0x41f2b459, /*  3.0338060379e+01 */
f32sb2  = 0x43a2e571, /*  3.2579251099e+02 */
f32sb3  = 0x44c01759, /*  1.5367296143e+03 */
f32sb4  = 0x4547fdbb, /*  3.1998581543e+03 */
f32sb5  = 0x451f90ce, /*  2.5530502930e+03 */
f32sb6  = 0x43ed43a7, /*  4.7452853394e+02 */
f32sb7  = 0xc1b38712; /* -2.2440952301e+01 */


static const sf_float64
f64erx  = 0x3feb0ac160000000, /* 8.45062911510467529297e-01 */
/*
 * Coefficients for approximation to  erf on [0,0.84375]
 */
f64efx8 = 0x3ff06eba8214db69, /*  1.02703333676410069053e+00 */
f64pp0  = 0x3fc06eba8214db68, /*  1.28379167095512558561e-01 */
f64pp1  = 0xbfd4cd7d691cb913, /* -3.25042107247001499370e-01 */
f64pp2  = 0xbf9d2a51dbd7194f, /* -2.84817495755985104766e-02 */
f64pp3  = 0xbf77a291236668e4, /* -5.77027029648944159157e-03 */
f64pp4  = 0xbef8ead6120016ac, /* -2.37630166566501626084e-05 */
f64qq1  = 0x3fd97779cddadc09, /*  3.97917223959155352819e-01 */
f64qq2  = 0x3fb0a54c5536ceba, /*  6.50222499887672944485e-02 */
f64qq3  = 0x3f74d022c4d36b0f, /*  5.08130628187576562776e-03 */
f64qq4  = 0x3f215dc9221c1a10, /*  1.32494738004321644526e-04 */
f64qq5  = 0xbed09c4342a26120, /* -3.96022827877536812320e-06 */
/*
 * Coefficients for approximation to  erf  in [0.84375,1.25]
 */
f64pa0  = 0xbf6359b8bef77538, /* -2.36211856075265944077e-03 */
f64pa1  = 0x3fda8d00ad92b34d, /*  4.14856118683748331666e-01 */
f64pa2  = 0xbfd7d240fbb8c3f1, /* -3.72207876035701323847e-01 */
f64pa3  = 0x3fd45fca805120e4, /*  3.18346619901161753674e-01 */
f64pa4  = 0xbfbc63983d3e28ec, /* -1.10894694282396677476e-01 */
f64pa5  = 0x3fa22a36599795eb, /*  3.54783043256182359371e-02 */
f64pa6  = 0xbf61bf380a96073f, /* -2.16637559486879084300e-03 */
f64qa1  = 0x3fbb3e6618eee323, /*  1.06420880400844228286e-01 */
f64qa2  = 0x3fe14af092eb6f33, /*  5.40397917702171048937e-01 */
f64qa3  = 0x3fb2635cd99fe9a7, /*  7.18286544141962662868e-02 */
f64qa4  = 0x3fc02660e763351f, /*  1.26171219808761642112e-01 */
f64qa5  = 0x3f8bedc26b51dd1c, /*  1.36370839120290507362e-02 */
f64qa6  = 0x3f888b545735151d, /*  1.19844998467991074170e-02 */
/*
 * Coefficients for approximation to  erfc in [1.25,1/0.35]
 */
f64ra0  = 0xbf843412600d6435, /* -9.86494403484714822705e-03 */
f64ra1  = 0xbfe63416e4ba7360, /* -6.93858572707181764372e-01 */
f64ra2  = 0xc0251e0441b0e726, /* -1.05586262253232909814e+01 */
f64ra3  = 0xc04f300ae4cba38d, /* -6.23753324503260060396e+01 */
f64ra4  = 0xc0644cb184282266, /* -1.62396669462573470355e+02 */
f64ra5  = 0xc067135cebccabb2, /* -1.84605092906711035994e+02 */
f64ra6  = 0xc054526557e4d2f2, /* -8.12874355063065934246e+01 */
f64ra7  = 0xc023a0efc69ac25c, /* -9.81432934416914548592e+00 */
f64sa1  = 0x4033a6b9bd707687, /*  1.96512716674392571292e+01 */
f64sa2  = 0x4061350c526ae721, /*  1.37657754143519042600e+02 */
f64sa3  = 0x407b290dd58a1a71, /*  4.34565877475229228821e+02 */
f64sa4  = 0x40842b1921ec2868, /*  6.45387271733267880336e+02 */
f64sa5  = 0x407ad02157700314, /*  4.29008140027567833386e+02 */
f64sa6  = 0x405b28a3ee48ae2c, /*  1.08635005541779435134e+02 */
f64sa7  = 0x401a47ef8e484a93, /*  6.57024977031928170135e+00 */
f64sa8  = 0xbfaeeff2ee749a62, /* -6.04244152148580987438e-02 */
/*
 * Coefficients for approximation to  erfc in [1/.35,28]
 */
f64rb0  = 0xbf84341239e86f4a, /* -9.86494292470009928597e-03 */
f64rb1  = 0xbfe993ba70c285de, /* -7.99283237680523006574e-01 */
f64rb2  = 0xc031c209555f995a, /* -1.77579549177547519889e+01 */
f64rb3  = 0xc064145d43c5ed98, /* -1.60636384855821916062e+02 */
f64rb4  = 0xc083ec881375f228, /* -6.37566443368389627722e+02 */
f64rb5  = 0xc09004616a2e5992, /* -1.02509513161107724954e+03 */
f64rb6  = 0xc07e384e9bdc383f, /* -4.83519191608651397019e+02 */
f64sb1  = 0x403e568b261d5190, /*  3.03380607434824582924e+01 */
f64sb2  = 0x40745cae221b9f0a, /*  3.25792512996573918826e+02 */
f64sb3  = 0x409802eb189d5118, /*  1.53672958608443695994e+03 */
f64sb4  = 0x40a8ffb7688c246a, /*  3.19985821950859553908e+03 */
f64sb5  = 0x40a3f219cedf3be6, /*  2.55305040643316442583e+03 */
f64sb6  = 0x407da874e79fe763, /*  4.74528541206955367215e+02 */
f64sb7  = 0xc03670e242712d62; /* -2.24409524465858183362e+01 */

static sf_result32f sf_float32_erfc1(sf_float32 x, sf_fpu_state fpu) {
    sf_float32 s, P, Q, tmp1, tmp2;

    SF_F32_OP(fpu, s, sf_float32_sub, sf_float32_abs(x), 0x3f800000); /* s = fabsf(x) - 1; */
    { /* P = f32pa0+s*(f32pa1+s*(f32pa2+s*(f32pa3+s*(f32pa4+s*(f32pa5+s*f32pa6))))); */
        SF_F32_OP(fpu, P, sf_float32_mul, s, f32pa6);
        SF_F32_OP(fpu, P, sf_float32_add, f32pa5, P);
        SF_F32_OP(fpu, P, sf_float32_mul, s, P);
        SF_F32_OP(fpu, P, sf_float32_add, f32pa4, P);
        SF_F32_OP(fpu, P, sf_float32_mul, s, P);
        SF_F32_OP(fpu, P, sf_float32_add, f32pa3, P);
        SF_F32_OP(fpu, P, sf_float32_mul, s, P);
        SF_F32_OP(fpu, P, sf_float32_add, f32pa2, P);
        SF_F32_OP(fpu, P, sf_float32_mul, s, P);
        SF_F32_OP(fpu, P, sf_float32_add, f32pa1, P);
        SF_F32_OP(fpu, P, sf_float32_mul, s, P);
        SF_F32_OP(fpu, P, sf_float32_add, f32pa0, P);
    }
    { /* Q = 1+s*(f32qa1+s*(f32qa2+s*(f32qa3+s*(f32qa4+s*(f32qa5+s*f32qa6))))); */
        SF_F32_OP(fpu, Q, sf_float32_mul, s, f32qa6);
        SF_F32_OP(fpu, Q, sf_float32_add, f32qa5, Q);
        SF_F32_OP(fpu, Q, sf_float32_mul, s, Q);
        SF_F32_OP(fpu, Q, sf_float32_add, f32qa4, Q);
        SF_F32_OP(fpu, Q, sf_float32_mul, s, Q);
        SF_F32_OP(fpu, Q, sf_float32_add, f32qa3, Q);
        SF_F32_OP(fpu, Q, sf_float32_mul, s, Q);
        SF_F32_OP(fpu, Q, sf_float32_add, f32qa2, Q);
        SF_F32_OP(fpu, Q, sf_float32_mul, s, Q);
        SF_F32_OP(fpu, Q, sf_float32_add, f32qa1, Q);
        SF_F32_OP(fpu, Q, sf_float32_mul, s, Q);
        SF_F32_OP(fpu, Q, sf_float32_add, 0x3f800000, Q);
    }
    { /* return 1 - f32erx - P/Q; */
        SF_F32_OP(fpu, tmp1, sf_float32_sub, 0x3f800000, f32erx);
        SF_F32_OP(fpu, tmp2, sf_float32_div, P, Q);
        SF_F32_OP(fpu, tmp1, sf_float32_sub, tmp1, tmp2);
        return (sf_result32f) { tmp1, fpu };
    }
}

static sf_result32f sf_float32_erfc2(sf_float32 x, sf_fpu_state fpu) {
    sf_float32 s, R, S, z, tmp1, tmp2, tmp3;

    uint32_t ix = x & 0x7fffffff;

    if (ix < 0x3fa00000)  /* |x| < 1.25 */
        return sf_float32_erfc1(x, fpu);

    x = sf_float32_abs(x); /* x = fabsf(x); */
    { /* s = 1/(x*x); */
        SF_F32_OP(fpu, s, sf_float32_mul, x, x);
        SF_F32_OP(fpu, s, sf_float32_div, 0x3f800000, s);
    }
    if (ix < 0x4036db6d) {   /* |x| < 1/0.35 */
        { /* R = f32ra0+s*(f32ra1+s*(f32ra2+s*(f32ra3+s*(f32ra4+s*(f32ra5+s*(f32ra6+s*f32ra7)))))); */
            SF_F32_OP(fpu, R, sf_float32_mul, s, f32ra7);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra6, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra5, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra4, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra3, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra2, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra1, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32ra0, R);
        }
        { /* S = 1.0f+s*(f32sa1+s*(f32sa2+s*(f32sa3+s*(f32sa4+s*(f32sa5+s*(f32sa6+s*(f32sa7+s*f32sa8))))))); */
            SF_F32_OP(fpu, S, sf_float32_mul, s, f32sa8);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa7, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa6, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa5, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa4, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa3, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa2, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sa1, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, 0x3f800000, S);
        }
    } else {                 /* |x| >= 1/0.35 */
        { /* R = f32rb0+s*(f32rb1+s*(f32rb2+s*(f32rb3+s*(f32rb4+s*(f32rb5+s*f32rb6))))); */
            SF_F32_OP(fpu, R, sf_float32_mul, s, f32rb6);
            SF_F32_OP(fpu, R, sf_float32_add, f32rb5, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32rb4, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32rb3, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32rb2, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32rb1, R);
            SF_F32_OP(fpu, R, sf_float32_mul, s, R);
            SF_F32_OP(fpu, R, sf_float32_add, f32rb0, R);
        }
        { /* S = 1.0f+s*(f32sb1+s*(f32sb2+s*(f32sb3+s*(f32sb4+s*(f32sb5+s*(f32sb6+s*f32sb7)))))); */
            SF_F32_OP(fpu, S, sf_float32_mul, s, f32sb7);
            SF_F32_OP(fpu, S, sf_float32_add, f32sb6, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sb5, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sb4, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sb3, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sb2, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, f32sb1, S);
            SF_F32_OP(fpu, S, sf_float32_mul, s, S);
            SF_F32_OP(fpu, S, sf_float32_add, 0x3f800000, S);
        }
    }
    /* GET_FLOAT_WORD(ix, x); SET_FLOAT_WORD(z, ix&0xffffe000); */
    z = x & 0xffffe000;
    { /* return expf(-z*z - 0.5625f) * expf((z-x)*(z+x) + R/S)/x; */
        SF_F32_OP(fpu, tmp1, sf_float32_mul, sf_float32_neg(z), z);
        SF_F32_OP(fpu, tmp1, sf_float32_sub, tmp1, 0x3f100000);
        SF_F32_OP(fpu, tmp1, sf_float32_exp, tmp1);

        SF_F32_OP(fpu, tmp2, sf_float32_sub, z, x);
        SF_F32_OP(fpu, tmp3, sf_float32_add, z, x);
        SF_F32_OP(fpu, tmp2, sf_float32_mul, tmp2, tmp3);
        SF_F32_OP(fpu, tmp3, sf_float32_div, R, S);
        SF_F32_OP(fpu, tmp2, sf_float32_add, tmp2, tmp3);
        SF_F32_OP(fpu, tmp2, sf_float32_exp, tmp2);

        SF_F32_OP(fpu, tmp1, sf_float32_mul, tmp1, tmp2);
        SF_F32_OP(fpu, tmp1, sf_float32_div, tmp1, x);
        return (sf_result32f) { tmp1, fpu };
    }
}

sf_result32f sf_float32_erf(sf_float32 x, sf_fpu_state fpu) {
    sf_float32 r, s, z, y, tmp1, tmp2;
    uint32_t ix;
    int32_t sign, i32tmp;

    sign = (int32_t)(x >> 31);
    ix = x & 0x7fffffff;
    if (ix >= 0x7f800000) {
        /* erf(nan)=nan, erf(+-inf)=+-1 */
        { /* return 1-2*sign + 1/x; */
            i32tmp = 1 - 2 * sign;
            SF_F32_OP(fpu, tmp1, sf_int32_to_float32, i32tmp);
            SF_F32_OP(fpu, tmp2, sf_float32_div, 0x3f800000, x);
            return sf_float32_add(tmp1, tmp2, fpu);
        }
    }
    if (ix < 0x3f580000) {  /* |x| < 0.84375 */
        if (ix < 0x31800000) {  /* |x| < 2**-28 */
            /* avoid underflow */
            { /* return 0.125f*(8*x + f32efx8*x); */
                SF_F32_OP(fpu, tmp1, sf_float32_mul, 0x41000000, x);
                SF_F32_OP(fpu, tmp2, sf_float32_mul, f32efx8, x);
                SF_F32_OP(fpu, tmp1, sf_float32_add, tmp1, tmp2);
                return sf_float32_mul(0x3e000000, tmp1, fpu);
            }
        }
        SF_F32_OP(fpu, z, sf_float32_mul, x, x); /* z = x*x; */
        { /* r = f32pp0+z*(f32pp1+z*(f32pp2+z*(f32pp3+z*f32pp4))); */
            SF_F32_OP(fpu, r, sf_float32_mul, z, f32pp4);
            SF_F32_OP(fpu, r, sf_float32_add, f32pp3, r);
            SF_F32_OP(fpu, r, sf_float32_mul, z, r);
            SF_F32_OP(fpu, r, sf_float32_add, f32pp2, r);
            SF_F32_OP(fpu, r, sf_float32_mul, z, r);
            SF_F32_OP(fpu, r, sf_float32_add, f32pp1, r);
            SF_F32_OP(fpu, r, sf_float32_mul, z, r);
            SF_F32_OP(fpu, r, sf_float32_add, f32pp0, r);
        }
        { /* s = 1+z*(f32qq1+z*(f32qq2+z*(f32qq3+z*(f32qq4+z*f32qq5)))); */
            SF_F32_OP(fpu, s, sf_float32_mul, z, f32qq5);
            SF_F32_OP(fpu, s, sf_float32_add, f32qq4, s);
            SF_F32_OP(fpu, s, sf_float32_mul, z, s);
            SF_F32_OP(fpu, s, sf_float32_add, f32qq3, s);
            SF_F32_OP(fpu, s, sf_float32_mul, z, s);
            SF_F32_OP(fpu, s, sf_float32_add, f32qq2, s);
            SF_F32_OP(fpu, s, sf_float32_mul, z, s);
            SF_F32_OP(fpu, s, sf_float32_add, f32qq1, s);
            SF_F32_OP(fpu, s, sf_float32_mul, z, s);
            SF_F32_OP(fpu, s, sf_float32_add, 0x3f800000, s);
        }
        SF_F32_OP(fpu, y, sf_float32_div, r, s); /* y = r/s; */
        { /* return x + x*y; */
            SF_F32_OP(fpu, tmp1, sf_float32_mul, x, y);
            return sf_float32_add(x, tmp1, fpu);
        }
    }
    if (ix < 0x40c00000) { /* |x| < 6 */
        /* y = 1 - erfc2(ix,x); */
        SF_F32_OP(fpu, y, sf_float32_erfc2, x);
        SF_F32_OP(fpu, y, sf_float32_sub, 0x3f800000, y);
    } else {
        SF_F32_OP(fpu, y, sf_float32_sub, 0x3f800000, 0x3800000); /* y = 1 - 0x1p-120f; */
    }
    return (sf_result32f) { sign ? sf_float32_neg(y) : y, fpu };
}

static sf_result64f sf_float64_erfc1(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 s, P, Q, tmp1, tmp2;

    SF_F64_OP(fpu, s, sf_float64_sub, sf_float64_abs(x),
              0x3ff0000000000000); /* s = fabs(x) - 1; */

    { /* P = f64pa0+s*(f64pa1+s*(f64pa2+s*(f64pa3+s*(f64pa4+s*(f64pa5+s*f64pa6))))); */
        SF_F64_OP(fpu, P, sf_float64_mul, s, f64pa6);
        SF_F64_OP(fpu, P, sf_float64_add, f64pa5, P);
        SF_F64_OP(fpu, P, sf_float64_mul, s, P);
        SF_F64_OP(fpu, P, sf_float64_add, f64pa4, P);
        SF_F64_OP(fpu, P, sf_float64_mul, s, P);
        SF_F64_OP(fpu, P, sf_float64_add, f64pa3, P);
        SF_F64_OP(fpu, P, sf_float64_mul, s, P);
        SF_F64_OP(fpu, P, sf_float64_add, f64pa2, P);
        SF_F64_OP(fpu, P, sf_float64_mul, s, P);
        SF_F64_OP(fpu, P, sf_float64_add, f64pa1, P);
        SF_F64_OP(fpu, P, sf_float64_mul, s, P);
        SF_F64_OP(fpu, P, sf_float64_add, f64pa0, P);
    }

    { /* Q = 1+s*(f64qa1+s*(f64qa2+s*(f64qa3+s*(f64qa4+s*(f64qa5+s*f64qa6))))); */
        SF_F64_OP(fpu, Q, sf_float64_mul, s, f64qa6);
        SF_F64_OP(fpu, Q, sf_float64_add, f64qa5, Q);
        SF_F64_OP(fpu, Q, sf_float64_mul, s, Q);
        SF_F64_OP(fpu, Q, sf_float64_add, f64qa4, Q);
        SF_F64_OP(fpu, Q, sf_float64_mul, s, Q);
        SF_F64_OP(fpu, Q, sf_float64_add, f64qa3, Q);
        SF_F64_OP(fpu, Q, sf_float64_mul, s, Q);
        SF_F64_OP(fpu, Q, sf_float64_add, f64qa2, Q);
        SF_F64_OP(fpu, Q, sf_float64_mul, s, Q);
        SF_F64_OP(fpu, Q, sf_float64_add, f64qa1, Q);
        SF_F64_OP(fpu, Q, sf_float64_mul, s, Q);
        SF_F64_OP(fpu, Q, sf_float64_add, 0x3ff0000000000000, Q);
    }

    { /* return 1 - f64erx - P/Q; */
        SF_F64_OP(fpu, tmp1, sf_float64_sub, 0x3ff0000000000000, f64erx);
        SF_F64_OP(fpu, tmp2, sf_float64_div, P, Q);
        SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, tmp2);
        return (sf_result64f) { tmp1, fpu };
    }
}

static sf_result64f sf_float64_erfc2(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 s, z, R, S, tmp1, tmp2, tmp3;

    uint64_t ix = x & 0x7fffffff00000000;

    if (ix < 0x3ff4000000000000)  /* |x| < 1.25 */
        return sf_float64_erfc1(x, fpu);

    x = sf_float64_abs(x);
    { /* s = 1/(x*x); */
        SF_F64_OP(fpu, s, sf_float64_mul, x, x);
        SF_F64_OP(fpu, s, sf_float64_div, 0x3ff0000000000000, s);
    }
    if (ix < 0x4006db6d00000000) {  /* |x| < 1/.35 ~ 2.85714 */
        { /*R = f64ra0+s*(f64ra1+s*(f64ra2+s*(f64ra3+s*(f64ra4+s*(f64ra5+s*(f64ra6+s*f64ra7)))))); */
            SF_F64_OP(fpu, R, sf_float64_mul, s, f64ra7);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra6, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra5, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra4, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra3, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra2, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra1, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64ra0, R);
        }
        { /* S = 1.0+s*(f64sa1+s*(f64sa2+s*(f64sa3+s*(f64sa4+s*(f64sa5+s*(f64sa6+s*(f64sa7+s*f64sa8))))))); */
            SF_F64_OP(fpu, S, sf_float64_mul, s, f64sa8);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa7, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa6, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa5, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa4, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa3, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa2, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sa1, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, 0x3ff0000000000000, S);
        }
    } else {                /* |x| > 1/.35 */
        { /* R = f64rb0+s*(f64rb1+s*(f64rb2+s*(f64rb3+s*(f64rb4+s*(f64rb5+s*f64rb6))))); */
            SF_F64_OP(fpu, R, sf_float64_mul, s, f64rb6);
            SF_F64_OP(fpu, R, sf_float64_add, f64rb5, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64rb4, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64rb3, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64rb2, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64rb1, R);
            SF_F64_OP(fpu, R, sf_float64_mul, s, R);
            SF_F64_OP(fpu, R, sf_float64_add, f64rb0, R);
        }
        { /* S = 1.0+s*(f64sb1+s*(f64sb2+s*(f64sb3+s*(f64sb4+s*(f64sb5+s*(f64sb6+s*f64sb7)))))); */
            SF_F64_OP(fpu, S, sf_float64_mul, s, f64sb7);
            SF_F64_OP(fpu, S, sf_float64_add, f64sb6, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sb5, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sb4, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sb3, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sb2, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, f64sb1, S);
            SF_F64_OP(fpu, S, sf_float64_mul, s, S);
            SF_F64_OP(fpu, S, sf_float64_add, 0x3ff0000000000000, S);
        }
    }

    z = x & 0xffffffff00000000; /* SET_LOW_WORD(z,0); */
    { /* return exp(-z*z-0.5625)*exp((z-x)*(z+x)+R/S)/x; */
        SF_F64_OP(fpu, tmp1, sf_float64_mul, sf_float64_neg(z), z);
        SF_F64_OP(fpu, tmp1, sf_float64_sub, tmp1, 0x3fe2000000000000);
        SF_F64_OP(fpu, tmp1, sf_float64_exp, tmp1);

        SF_F64_OP(fpu, tmp2, sf_float64_sub, z, x);
        SF_F64_OP(fpu, tmp3, sf_float64_add, z, x);
        SF_F64_OP(fpu, tmp2, sf_float64_mul, tmp2, tmp3);
        SF_F64_OP(fpu, tmp3, sf_float64_div, R, S);
        SF_F64_OP(fpu, tmp2, sf_float64_add, tmp2, tmp3);
        SF_F64_OP(fpu, tmp2, sf_float64_exp, tmp2);

        SF_F64_OP(fpu, tmp1, sf_float64_mul, tmp1, tmp2);
        SF_F64_OP(fpu, tmp1, sf_float64_div, tmp1, x);
        return (sf_result64f) { tmp1, fpu };
    }
}

sf_result64f sf_float64_erf(sf_float64 x, sf_fpu_state fpu) {
    sf_float64 r, s, z, y, tmp1, tmp2;
    uint64_t ix;
    int64_t sign, i64tmp;

    sign = (int64_t)(x >> 63);
    ix = x & 0x7fffffff00000000;

    if (ix >= 0x7ff0000000000000) {
        /* erf(nan)=nan, erf(+-inf)=+-1 */
        { /* return 1-2*sign + 1/x; */
            i64tmp = 1 - 2 * sign;
            SF_F64_OP(fpu, tmp1, sf_int64_to_float64, i64tmp);
            SF_F64_OP(fpu, tmp2, sf_float64_div, 0x3ff0000000000000, x);
            return sf_float64_add(tmp1, tmp2, fpu);
        }
    }
    if (ix < 0x3feb000000000000) {  /* |x| < 0.84375 */
        if (ix < 0x3e30000000000000) {  /* |x| < 2**-28 */
            /* avoid underflow */
            { /* return 0.125*(8*x + f64efx8*x); */
                SF_F64_OP(fpu, tmp1, sf_float64_mul, 0x4020000000000000, x);
                SF_F64_OP(fpu, tmp2, sf_float64_mul, f64efx8, x);
                SF_F64_OP(fpu, tmp1, sf_float64_add, tmp1, tmp2);
                return sf_float64_mul(0x3fc0000000000000, tmp1, fpu);
            }
        }
        SF_F64_OP(fpu, z, sf_float64_mul, x, x); /* z = x*x; */
        { /* r = f64pp0+z*(f64pp1+z*(f64pp2+z*(f64pp3+z*f64pp4))); */
            SF_F64_OP(fpu, r, sf_float64_mul, z, f64pp4);
            SF_F64_OP(fpu, r, sf_float64_add, f64pp3, r);
            SF_F64_OP(fpu, r, sf_float64_mul, z, r);
            SF_F64_OP(fpu, r, sf_float64_add, f64pp2, r);
            SF_F64_OP(fpu, r, sf_float64_mul, z, r);
            SF_F64_OP(fpu, r, sf_float64_add, f64pp1, r);
            SF_F64_OP(fpu, r, sf_float64_mul, z, r);
            SF_F64_OP(fpu, r, sf_float64_add, f64pp0, r);
        }
        { /* s = 1.0+z*(f64qq1+z*(f64qq2+z*(f64qq3+z*(f64qq4+z*f64qq5)))); */
            SF_F64_OP(fpu, s, sf_float64_mul, z, f64qq5);
            SF_F64_OP(fpu, s, sf_float64_add, f64qq4, s);
            SF_F64_OP(fpu, s, sf_float64_mul, z, s);
            SF_F64_OP(fpu, s, sf_float64_add, f64qq3, s);
            SF_F64_OP(fpu, s, sf_float64_mul, z, s);
            SF_F64_OP(fpu, s, sf_float64_add, f64qq2, s);
            SF_F64_OP(fpu, s, sf_float64_mul, z, s);
            SF_F64_OP(fpu, s, sf_float64_add, f64qq1, s);
            SF_F64_OP(fpu, s, sf_float64_mul, z, s);
            SF_F64_OP(fpu, s, sf_float64_add, 0x3ff0000000000000, s);
        }
        SF_F64_OP(fpu, y, sf_float64_div, r, s); /* y = r/s; */
        { /* return x + x*y; */
            SF_F64_OP(fpu, tmp1, sf_float64_mul, x, y);
            return sf_float64_add(x, tmp1, fpu);
        }
    }

    if (ix < 0x4018000000000000) {  /* 0.84375 <= |x| < 6 */
        { /* y = 1 - erfc2(x, fpu); */
            SF_F64_OP(fpu, y, sf_float64_erfc2, x);
            SF_F64_OP(fpu, y, sf_float64_sub, 0x3ff0000000000000, y);
        }
    } else {
        SF_F64_OP(fpu, y, sf_float64_sub, 0x3ff0000000000000, 0x10000000000000); /* y = 1 - 0x1p-1022; */
    }

    return (sf_result64f) { sign ? sf_float64_neg(y) : y, fpu }; /* return sign ? -y : y; */
}
