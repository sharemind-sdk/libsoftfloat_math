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

#ifndef INTERNAL_MACROS_H
#define INTERNAL_MACROS_H

#define SF_OP(type,fpu,var,op,...) \
    do { \
        type __rv = op(__VA_ARGS__, (fpu)); \
        (var) = __rv.result; \
        (fpu) = __rv.fpu_state; \
    } while (0)

#define SF_FLAG_OP(fpu,var,op,...) \
    SF_OP(sf_resultFlag, (fpu), (var), (op), __VA_ARGS__);

#define SF_F32_OP(fpu,var,op,...) \
    SF_OP(sf_result32f, (fpu), (var), (op), __VA_ARGS__);

#define SF_I32_OP(fpu,var,op,...) \
    SF_OP(sf_result32i, (fpu), (var), (op), __VA_ARGS__);

#define SF_F64_OP(fpu,var,op,...) \
    SF_OP(sf_result64f, (fpu), (var), (op), __VA_ARGS__);

#define SF_I64_OP(fpu,var,op,...) \
    SF_OP(sf_result64i, (fpu), (var), (op), __VA_ARGS__);

#endif /* INTERNAL_MACROS_H */