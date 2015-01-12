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

/* origin: musl-libc src/math/scalbn.c */
/* origin: musl-libc src/math/scalbnf.c */
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

#include <stdint.h>
#include "internal/macros.h"
#include "sf_scalbn.h"


sf_result32f sf_float32_scalbn(sf_float32 x, sf_int32 n, sf_fpu_state fpu) {
    sf_float32 u, y = x;

    if (n > 127) {
        SF_F32_OP(fpu, y, sf_float32_mul, y, 0x7f000000); /* y *= 0x1p127f; */
        n -= 127;
        if (n > 127) {
            SF_F32_OP(fpu, y, sf_float32_mul, y, 0x7f000000); /* y *= 0x1p127f; */
            n -= 127;
            if (n > 127)
                n = 127;
        }
    } else if (n < -126) {
        SF_F32_OP(fpu, y, sf_float32_mul, y, 0x800000); /* y *= 0x1p-126f; */
        n += 126;
        if (n < -126) {
            SF_F32_OP(fpu, y, sf_float32_mul, y, 0x800000); /* y *= 0x1p-126f; */
            n += 126;
            if (n < -126)
                n = -126;
        }
    }

    u = (sf_float32)(0x7f+n) << 23;
    return sf_float32_mul(y, u, fpu); /* x = y * u; */
}

sf_result64f sf_float64_scalbn(sf_float64 x, sf_int64 n, sf_fpu_state fpu) {
    sf_float64 u, y = x;

    if (n > 1023) {
        SF_F64_OP(fpu, y, sf_float64_mul, y, 0x7fe0000000000000); /* y *= 0x1p1023; */
        n -= 1023;
        if (n > 1023) {
            SF_F64_OP(fpu, y, sf_float64_mul, y, 0x7fe0000000000000); /* y *= 0x1p1023; */
            n -= 1023;
            if (n > 1023)
                n = 1023;
        }
    } else if (n < -1022) {
        SF_F64_OP(fpu, y, sf_float64_mul, y, 0x10000000000000); /* y *= 0x1p-1022; */
        n += 1022;
        if (n < -1022) {
            SF_F64_OP(fpu, y, sf_float64_mul, y, 0x10000000000000); /* y *= 0x1p-1022; */
            n += 1022;
            if (n < -1022)
                n = -1022;
        }
    }

    u = (sf_float64)(0x3ff+n) << 52;
    return sf_float64_mul(y, u, fpu); /* x = y * u; */
}
