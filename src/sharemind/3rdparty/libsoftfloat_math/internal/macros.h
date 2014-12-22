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
