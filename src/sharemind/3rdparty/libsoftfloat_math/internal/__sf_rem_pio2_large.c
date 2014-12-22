/* origin: musl-libc src/math/__rem_pio2_large.c */
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
/* origin: FreeBSD /usr/src/lib/msun/src/k_rem_pio2.c */
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
 * __rem_pio2_large(x,y,e0,nx,prec)
 * double x[],y[]; int e0,nx,prec;
 *
 * __rem_pio2_large return the last three digits of N with
 *              y = x - N*pi/2
 * so that |y| < pi/2.
 *
 * The method is to compute the integer (mod 8) and fraction parts of
 * (2/pi)*x without doing the full multiplication. In general we
 * skip the part of the product that are known to be a huge integer (
 * more accurately, = 0 mod 8 ). Thus the number of operations are
 * independent of the exponent of the input.
 *
 * (2/pi) is represented by an array of 24-bit integers in ipio2[].
 *
 * Input parameters:
 *      x[]     The input value (must be positive) is broken into nx
 *              pieces of 24-bit integers in double precision format.
 *              x[i] will be the i-th 24 bit of x. The scaled exponent
 *              of x[0] is given in input parameter e0 (i.e., x[0]*2^e0
 *              match x's up to 24 bits.
 *
 *              Example of breaking a double positive z into x[0]+x[1]+x[2]:
 *                      e0 = ilogb(z)-23
 *                      z  = scalbn(z,-e0)
 *              for i = 0,1,2
 *                      x[i] = floor(z)
 *                      z    = (z-x[i])*2**24
 *
 *
 *      y[]     ouput result in an array of double precision numbers.
 *              The dimension of y[] is:
 *                      24-bit  precision       1
 *                      53-bit  precision       2
 *                      64-bit  precision       2
 *                      113-bit precision       3
 *              The actual value is the sum of them. Thus for 113-bit
 *              precison, one may have to do something like:
 *
 *              long double t,w,r_head, r_tail;
 *              t = (long double)y[2] + (long double)y[1];
 *              w = (long double)y[0];
 *              r_head = t+w;
 *              r_tail = w - (r_head - t);
 *
 *      e0      The exponent of x[0]. Must be <= 16360 or you need to
 *              expand the ipio2 table.
 *
 *      nx      dimension of x[]
 *
 *      prec    an integer indicating the precision:
 *                      0       24  bits (single)
 *                      1       53  bits (double)
 *                      2       64  bits (extended)
 *                      3       113 bits (quad)
 *
 * External function:
 *      double scalbn(), floor();
 *
 *
 * Here is the description of some local variables:
 *
 *      jk      jk+1 is the initial number of terms of ipio2[] needed
 *              in the computation. The minimum and recommended value
 *              for jk is 3,4,4,6 for single, double, extended, and quad.
 *              jk+1 must be 2 larger than you might expect so that our
 *              recomputation test works. (Up to 24 bits in the integer
 *              part (the 24 bits of it that we compute) and 23 bits in
 *              the fraction part may be lost to cancelation before we
 *              recompute.)
 *
 *      jz      local integer variable indicating the number of
 *              terms of ipio2[] used.
 *
 *      jx      nx - 1
 *
 *      jv      index for pointing to the suitable ipio2[] for the
 *              computation. In general, we want
 *                      ( 2^e0*x[0] * ipio2[jv-1]*2^(-24jv) )/8
 *              is an integer. Thus
 *                      e0-3-24*jv >= 0 or (e0-3)/24 >= jv
 *              Hence jv = max(0,(e0-3)/24).
 *
 *      jp      jp+1 is the number of terms in PIo2[] needed, jp = jk.
 *
 *      q[]     double array with integral value, representing the
 *              24-bits chunk of the product of x and 2/pi.
 *
 *      q0      the corresponding exponent of q[0]. Note that the
 *              exponent for q[i] would be q0-24*i.
 *
 *      PIo2[]  double precision array, obtained by cutting pi/2
 *              into 24 bits chunks.
 *
 *      f[]     ipio2[] in floating point
 *
 *      iq[]    integer array by breaking up q[] in 24-bits chunk.
 *
 *      fq[]    final product of x*(2/pi) in fq[0],..,fq[jk]
 *
 *      ih      integer. If >0 it indicates q[] is >= 0.5, hence
 *              it also indicates the *sign* of the result.
 *
 */
/*
 * Constants:
 * The hexadecimal values are the intended ones for the following
 * constants. The decimal values may be used, provided that the
 * compiler will convert from decimal to binary accurately enough
 * to produce the hexadecimal values shown.
 */

#include "__sf_rem_pio2_large.h"
#include "assert.h"
#include "macros.h"
#include "../sf_floor.h"
#include "../sf_scalbn.h"


static const int64_t init_jk[] = { 3, 4, 4, 6 }; /* initial value for jk */

/*
 * Table of constants for 2/pi, 396 Hex digits (476 decimal) of 2/pi
 *
 *              integer array, contains the (24*i)-th to (24*i+23)-th
 *              bit of 2/pi after binary point. The corresponding
 *              floating value is
 *
 *                      ipio2[i] * 2^(-24(i+1)).
 *
 * NB: This table must have at least (e0-3)/24 + jk terms.
 *     For quad precision (e0 <= 16360, jk = 6), this is 686.
 */
static const int32_t ipio2[] = {
    0xa2f983, 0x6e4e44, 0x1529fc, 0x2757d1, 0xf534dd, 0xc0db62,
    0x95993c, 0x439041, 0xfe5163, 0xabdebb, 0xc561b7, 0x246e3a,
    0x424dd2, 0xe00649, 0x2eea09, 0xd1921c, 0xfe1deb, 0x1cb129,
    0xa73ee8, 0x8235f5, 0x2ebb44, 0x84e99c, 0x7026b4, 0x5f7e41,
    0x3991d6, 0x398353, 0x39f49c, 0x845f8b, 0xbdf928, 0x3b1ff8,
    0x97ffde, 0x05980f, 0xef2f11, 0x8b5a0a, 0x6d1f6d, 0x367ecf,
    0x27cb09, 0xb74f46, 0x3f669e, 0x5fea2d, 0x7527ba, 0xc7ebe5,
    0xf17b3d, 0x0739f7, 0x8a5292, 0xea6bfb, 0x5fb11f, 0x8d5d08,
    0x560330, 0x46fc7b, 0x6babf0, 0xcfbc20, 0x9af436, 0x1da9e3,
    0x91615e, 0xe61b08, 0x659985, 0x5f14a0, 0x68408d, 0xffd880,
    0x4d7327, 0x310606, 0x1556ca, 0x73a8c9, 0x60e27b, 0xc08c6b,
};

static const sf_float64 PIo2[] = {
    0x3ff921fb40000000, /* 1.57079625129699707031e+00 */
    0x3e74442d00000000, /* 7.54978941586159635335e-08 */
    0x3cf8469880000000, /* 5.39030252995776476554e-15 */
    0x3b78cc5160000000, /* 3.28200341580791294123e-22 */
    0x39f01b8380000000, /* 1.27065575308067607349e-29 */
    0x387a252040000000, /* 1.22933308981111328932e-36 */
    0x36e3822280000000, /* 2.73370053816464559624e-44 */
    0x3569f31d00000000  /* 2.16741683877804819444e-51 */
};

sf_result32i __sf_float64_rem_pio2_large(sf_float64 * x, sf_float64 * y,
        int64_t e0, int64_t nx, int64_t prec, sf_fpu_state fpu)
{
    int64_t jz, jx, jv, jp, jk, carry, iq[20], i, j, k, m, q0, ih;
    sf_int32 i32tmp;
    sf_int32 n;
    sf_float64 z, fw, f[20], fq[20], q[20], tmp;
    sf_flag ftmp;

    /* initialize jk*/
    jk = init_jk[prec];
    jp = jk;

    /* determine jx,jv,q0, note that 3>q0 */
    jx = nx-1;
    jv = (e0 - 3) / 24;
    if (jv < 0)
        jv = 0;
    q0 = e0 - 24 * (jv + 1);

    /* set up f[0] to f[jx+jk] where f[jx+jk] = ipio2[jv+jk] */
    j = jv - jx; m = jx + jk;
    for (i = 0; i <= m; i++, j++) {
        /* f[i] = j<0 ? 0.0 : (double)ipio2[j]; */
        f[i] = j < 0 ? 0x0 : sf_int32_to_float64(ipio2[j]);
    }

    /* compute q[0],q[1],...q[jk] */
    for (i = 0; i <= jk; i++) {
        for (j = 0, fw = 0x0; j <= jx; j++) {
            /* fw += x[j]*f[jx+i-j]; */
            SF_F64_OP(fpu, tmp, sf_float64_mul, x[j], f[jx+i-j]);
            SF_F64_OP(fpu, fw, sf_float64_add, fw, tmp);
        }
        q[i] = fw;
    }

    jz = jk;
recompute:
    /* distill q[] into iq[] reversingly */
    for (i = 0, j = jz, z = q[jz]; j > 0; i++, j--) {
        { /* fw = (double)(int32_t)(0x1p-24*z); */
            SF_F64_OP(fpu, fw, sf_float64_mul, 0x3e70000000000000, z);
            SF_I32_OP(fpu, i32tmp, sf_float64_to_int32_round_to_zero, fw);
            fw = sf_int32_to_float64(i32tmp);
        }
        { /* iq[i] = (int32_t)(z - 0x1p24*fw); */
            SF_F64_OP(fpu, tmp, sf_float64_mul, 0x4170000000000000, fw);
            SF_F64_OP(fpu, tmp, sf_float64_sub, z, tmp);
            SF_I64_OP(fpu, iq[i], sf_float64_to_int64_round_to_zero, tmp);
        }
        SF_F64_OP(fpu, z, sf_float64_add, q[j-1], fw); /* z = q[j-1]+fw; */
    }

    /* compute n */
    SF_F64_OP(fpu, z, sf_float64_scalbn, z, q0); /* z = scalbn(z,q0); */ /* actual value of z */
    { /* z -= 8.0*floor(z*0.125); */ /* trim off integer >= 8 */
        SF_F64_OP(fpu, tmp, sf_float64_mul, z, 0x3fc0000000000000);
        SF_F64_OP(fpu, tmp, sf_float64_floor, tmp);
        SF_F64_OP(fpu, tmp, sf_float64_mul, 0x4020000000000000, tmp);
        SF_F64_OP(fpu, z, sf_float64_sub, z, tmp);
    }
    SF_I32_OP(fpu, n, sf_float64_to_int32_round_to_zero, z); /* n  = (int32_t)z; */
    { /* z -= (double)n; */
        tmp = sf_int32_to_float64(n);
        SF_F64_OP(fpu, z, sf_float64_sub, z, tmp);
    }
    ih = 0;
    if (q0 > 0) {  /* need iq[jz-1] to determine n */
        i  = iq[jz-1] >> (24-q0); n += i;
        iq[jz-1] -= i << (24-q0);
        ih = iq[jz-1] >> (23-q0);
    } else if (q0 == 0) {
        ih = iq[jz-1] >> 23;
    } else {
        SF_FLAG_OP(fpu, ftmp, sf_float64_le, 0x3fe0000000000000, z);
        if (ftmp) { /* z >= 0.5 */
            ih = 2;
        }
    }

    if (ih > 0) {  /* q > 0.5 */
        n += 1; carry = 0;
        for (i = 0; i < jz; i++) {  /* compute 1-q */
            j = iq[i];
            if (carry == 0) {
                if (j != 0) {
                    carry = 1;
                    iq[i] = 0x1000000 - j;
                }
            } else {
                iq[i] = 0xffffff - j;
            }
        }
        if (q0 > 0) {  /* rare case: chance is 1 in 12 */
            switch (q0) {
            case 1:
                iq[jz-1] &= 0x7fffff;
                break;
            case 2:
                iq[jz-1] &= 0x3fffff;
                break;
            default:
                assert(0);
                break;
            }
        }
        if (ih == 2) {
            SF_F64_OP(fpu, z, sf_float64_sub, 0x3ff0000000000000, z); /* z = 1.0 - z; */
            if (carry != 0) {
                /* z -= scalbn(1.0,q0); */
                SF_F64_OP(fpu, tmp, sf_float64_scalbn, 0x3ff0000000000000, q0);
                SF_F64_OP(fpu, z, sf_float64_sub, z, tmp);
            }
        }
    }

    /* check if recomputation is needed */
    SF_FLAG_OP(fpu, ftmp, sf_float64_eq, z, 0x0);
    if (ftmp) { /* z == 0.0 */
        j = 0;
        for (i = jz-1; i >= jk; i--)
            j |= iq[i];
        if (j == 0) {  /* need recomputation */
            for (k = 1; iq[jk-k] == 0; k++);  /* k = no. of terms needed */

            for (i= jz+1; i <= jz+k; i++) {  /* add q[jz+1] to q[jz+k] */
                f[jx+i] = sf_int32_to_float64(ipio2[jv+i]); /* f[jx+i] = (double)ipio2[jv+i]; */
                for (j = 0, fw = 0x0; j <= jx; j++) {
                    /* fw += x[j]*f[jx+i-j]; */
                    SF_F64_OP(fpu, tmp, sf_float64_mul, x[j], f[jx+i-j]);
                    SF_F64_OP(fpu, fw, sf_float64_add, fw, tmp);
                }
                q[i] = fw;
            }
            jz += k;
            goto recompute;
        }
    }

    /* chop off zero terms */
    SF_FLAG_OP(fpu, ftmp, sf_float64_eq, z, 0x0);
    if (ftmp) { /* z == 0.0 */
        jz -= 1;
        q0 -= 24;
        while (iq[jz] == 0) {
            jz--;
            q0 -= 24;
        }
    } else { /* break z into 24-bit if necessary */
        SF_F64_OP(fpu, z, sf_float64_scalbn, z, -q0); /* z = scalbn(z,-q0); */
        SF_FLAG_OP(fpu, ftmp, sf_float64_le, 0x4170000000000000, z);
        if (ftmp) { /* z >= 0x1p24 */
            { /* fw = (double)(int32_t)(0x1p-24*z); */
                SF_F64_OP(fpu, tmp, sf_float64_mul, 0x3e70000000000000, z);
                SF_I32_OP(fpu, i32tmp, sf_float64_to_int32_round_to_zero, tmp);
                fw = sf_int32_to_float64(i32tmp);
            }
            { /* iq[jz] = (int32_t)(z - 0x1p24*fw); */
                SF_F64_OP(fpu, tmp, sf_float64_mul, 0x4170000000000000, fw);
                SF_F64_OP(fpu, tmp, sf_float64_sub, z, tmp);
                SF_I64_OP(fpu, iq[jz], sf_float64_to_int64_round_to_zero, tmp);
            }
            jz += 1;
            q0 += 24;
            SF_I64_OP(fpu, iq[jz], sf_float64_to_int64_round_to_zero, fw); /* iq[jz] = (int32_t)fw; */
        } else {
            SF_I64_OP(fpu, iq[jz], sf_float64_to_int64_round_to_zero, z); /* iq[jz] = (int32_t)z; */
        }
    }

    /* convert integer "bit" chunk to floating-point value */
    SF_F64_OP(fpu, fw, sf_float64_scalbn, 0x3ff0000000000000, q0); /* fw = scalbn(1.0,q0); */
    for (i = jz; i >= 0; i--) {
        { /* q[i] = fw*(double)iq[i]; */
            SF_F64_OP(fpu, tmp, sf_int64_to_float64, iq[i]);
            SF_F64_OP(fpu, q[i], sf_float64_mul, fw, tmp);
        }
        SF_F64_OP(fpu, fw, sf_float64_mul, fw, 0x3e70000000000000); /* fw *= 0x1p-24; */
    }

    /* compute PIo2[0,...,jp]*q[jz,...,0] */
    for (i = jz; i >= 0; i--) {
        for (fw = 0x0, k=0; k <= jp && k <= jz-i; k++) {
            /* fw += PIo2[k]*q[i+k]; */
            SF_F64_OP(fpu, tmp, sf_float64_mul, PIo2[k], q[i+k]);
            SF_F64_OP(fpu, fw, sf_float64_add, fw, tmp);
        }
        fq[jz-i] = fw;
    }

    /* compress fq[] into y[] */
    switch (prec) {
    case 0:
        fw = 0x0;
        for (i = jz; i >= 0; i--) {
            SF_F64_OP(fpu, fw, sf_float64_add, fw, fq[i]); /* fw += fq[i]; */
        }
        y[0] = ih == 0 ? fw : sf_float64_neg(fw);
        break;
    case 1:
    case 2:
        fw = 0x0;
        for (i = jz; i >= 0; i--) {
            SF_F64_OP(fpu, fw, sf_float64_add, fw, fq[i]); /* fw += fq[i]; */
        }
        y[0] = ih == 0 ? fw : sf_float64_neg(fw);
        SF_F64_OP(fpu, fw, sf_float64_sub, fq[0], fw); /* fw = fq[0]-fw; */
        for (i=1; i<=jz; i++) {
            SF_F64_OP(fpu, fw, sf_float64_add, fw, fq[i]); /* fw += fq[i]; */
        }
        y[1] = ih == 0 ? fw : sf_float64_neg(fw);
        break;
    case 3:  /* painful */
        for (i = jz; i > 0; i--) {
            SF_F64_OP(fpu, fw, sf_float64_add, fq[i-1], fq[i]); /* fw = fq[i-1]+fq[i]; */
            { /* fq[i] += fq[i-1]-fw; */
                SF_F64_OP(fpu, tmp, sf_float64_sub, fq[i-1], fw);
                SF_F64_OP(fpu, fq[i], sf_float64_add, fq[i], tmp);
            }
            fq[i-1] = fw;
        }
        for (i = jz; i > 1; i--) {
            SF_F64_OP(fpu, fw, sf_float64_add, fq[i-1], fq[i]); /* fw = fq[i-1]+fq[i]; */
            { /* fq[i] += fq[i-1]-fw; */
                SF_F64_OP(fpu, tmp, sf_float64_sub, fq[i-1], fw);
                SF_F64_OP(fpu, fq[i], sf_float64_add, fq[i], tmp);
            }
            fq[i-1] = fw;
        }
        for (fw = 0x0, i = jz; i >= 2; i--) {
            SF_F64_OP(fpu, fw, sf_float64_add, fw, fq[i]); /* fw += fq[i]; */
        }
        if (ih == 0) {
            y[0] = fq[0];
            y[1] = fq[1];
            y[2] = fw;
        } else {
            y[0] = sf_float64_neg(fq[0]); /* y[0] = -fq[0]; */
            y[1] = sf_float64_neg(fq[1]); /* y[1] = -fq[1]; */
            y[2] = sf_float64_neg(fw);    /* y[2] = -fw;    */
        }
    default:
        assert(0);
        break;
    }

    return (sf_result32i) { n & 0x7, fpu };
}
