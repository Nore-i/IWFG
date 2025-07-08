/* mex_stub.h – MINIMAL shim so code that #includes <mex.h> can
 * be compiled outside MATLAB.  Add more functions/macros here if
 * the rest of your source needs them.
 */
#ifndef MEX_STUB_H
#define MEX_STUB_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef ptrdiff_t  mwSize;
typedef ptrdiff_t  mwIndex;
typedef int32_t    int32_T;
typedef uint32_t   mxClassID;

/* Replace MATLAB memory helpers with libc versions */
#define mxMalloc      malloc
#define mxCalloc      calloc
#define mxRealloc     realloc
#define mxFree        free

/* Simple fatal-error stub */
static inline void mexErrMsgTxt(const char *msg)
{
    fprintf(stderr, "MEX-stub error: %s\n", msg);
    exit(EXIT_FAILURE);
}

#endif /* MEX_STUB_H */
