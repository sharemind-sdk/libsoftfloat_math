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

/* origin: musl-libc src/math/ceil.c */
/* origin: musl-libc src/math/ceilf.c */
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

#include "sf_ceil.h"

#include "internal/macros.h"


static const sf_float64 toint = 0x4330000000000000; /* 1 / 2.22044604925031308085e-16 */

sf_result32f sf_float32_ceil(sf_float32 x, sf_fpu_state fpu) {
    int32_t e = (int32_t)(x >> 23 & 0xff) - 0x7f;
    uint32_t m;
    sf_float32 tmp;
    (void) tmp;

    if (e >= 23)
        return (sf_result32f) { x, fpu };
    if (e >= 0) {
        m = (uint32_t)0x007fffff >> e;
        if ((x & m) == 0)
            return (sf_result32f) { x, fpu };
        SF_F32_OP(fpu, tmp, sf_float32_add, x, 0x7b800000); /* FORCE_EVAL(x + 0x1p120f); */
        if (x >> 31 == 0)
            x += m;
        x &= ~m;
    } else {
        SF_F32_OP(fpu, tmp, sf_float32_add, x, 0x7b800000); /* FORCE_EVAL(x + 0x1p120f); */
        if (x >> 31) {
            x = 0x80000000; /* x = -0.0; */
        } else if (x << 1) {
            x = 0x3f800000; /* x = 1.0 */
        }
    }
    return (sf_result32f) { x, fpu };
}

sf_result64f sf_float64_ceil(sf_float64 x, sf_fpu_state fpu) {
    int64_t e = x >> 52 & 0x7ff;
    sf_float64 y, tmp;
    sf_flag ftmp;

    SF_FLAG_OP(fpu, ftmp, sf_float64_eq, x, 0x0);
    if (e >= 0x3ff + 52 || ftmp) /* e >= 0x3ff + 52 || x == 0 */
        return (sf_result64f) { x, fpu };
    /* y = int(x) - x, where int(x) is an integer neighbor of x */
    if (x >> 63) {
        /* y = x - toint + toint - x; */
        SF_F64_OP(fpu, y, sf_float64_sub, x, toint);
        SF_F64_OP(fpu, y, sf_float64_add, y, toint);
        SF_F64_OP(fpu, y, sf_float64_sub, y, x);
    } else {
        /* y = x + toint - toint - x; */
        SF_F64_OP(fpu, y, sf_float64_add, x, toint);
        SF_F64_OP(fpu, y, sf_float64_sub, y, toint);
        SF_F64_OP(fpu, y, sf_float64_sub, y, x);
    }
    /* special case because of non-nearest rounding modes */
    if (e <= 0x3ff - 1) {
        /* FORCE_EVAL(y); */
        return (sf_result64f) { x >> 63 ? 0x8000000000000000 : 0x3ff0000000000000, fpu }; /* return x >> 63 ? -0.0 : 1; */
    }
    SF_FLAG_OP(fpu, ftmp, sf_float64_lt, y, 0x0);
    if (ftmp) {/* y < 0 */
        /* return x + y + 1; */
        SF_F64_OP(fpu, tmp, sf_float64_add, x, y);
        return sf_float64_add(tmp, 0x3ff0000000000000, fpu);
    }

    return sf_float64_add(x, y, fpu); /* return x + y; */
}
