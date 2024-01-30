#include "tsPredSucc.cpp"
#include "tsBounds.cpp"
#include "tsInfo.cpp"
#include "tsUtil.cpp"
#include "tsSetOp.cpp"
#include "tsAri.cpp"
#include "tsStdFun.cpp"

size_t ciTestSetSize[2] = {8, 21};
unsigned int ciTestSet[][2][4] = {

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x000, 0x00000, 0x00000000 },
  { 0x0, 0x000, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x400, 0x8cccc, 0xcccccccd },
  { 0x1, 0x400, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x400, 0x26666, 0x66666666 },
  { 0x0, 0x000, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x401, 0x00000, 0x00000000 },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x000, 0x00000, 0x00000000 },
  { 0x0, 0x401, 0xacccc, 0xcccccccd }
},

{
  { 0x0, 0x400, 0x00000, 0x00000000 },
  { 0x0, 0x401, 0x46666, 0x66666666 }
},

#if defined(FILIB_EXTENDED)
{
  { 0x0, 0x7ff, 0x80000, 0x00000000 },
  { 0x0, 0x7ff, 0x80000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x000, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x000, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x1, 0x7fe, 0xfffff, 0xffffffff },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

#endif
};


size_t srTestSetSize[2] = {16, 31};
unsigned int srTestSet1[][2][4] = {

{
  { 0x1, 0x3ff, 0xb3333, 0x33333333 },
  { 0x0, 0x400, 0x33333, 0x33333333 }
},

{
  { 0x0, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x400, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x3ff, 0x19999, 0x9999999a },
  { 0x0, 0x401, 0x40000, 0x00000000 }
},

{
  { 0x1, 0x400, 0x00000, 0x00000000 },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x000, 0x00000, 0x00000000 },
  { 0x0, 0x400, 0x40000, 0x00000000 }
},

{
  { 0x1, 0x3ff, 0xe6666, 0x66666666 },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x400, 0xb3333, 0x33333333 },
  { 0x0, 0x400, 0xf3333, 0x33333333 }
},

{
  { 0x1, 0x3ff, 0x19999, 0x9999999a },
  { 0x0, 0x3ff, 0x33333, 0x33333333 }
},

{
  { 0x0, 0x3fd, 0xae147, 0xae147ae1 },
  { 0x0, 0x404, 0x50000, 0x00000000 }
},

{
  { 0x1, 0x402, 0x40000, 0x00000000 },
  { 0x1, 0x401, 0x40000, 0x00000000 }
},

{
  { 0x0, 0x400, 0x40000, 0x00000000 },
  { 0x0, 0x400, 0x40000, 0x00000000 }
},

{
  { 0x0, 0x000, 0x00000, 0x00000000 },
  { 0x0, 0x402, 0xa6666, 0x66666666 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x400, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x404, 0x50000, 0x00000000 },
  { 0x0, 0x404, 0x50000, 0x00000000 }
},

{
  { 0x0, 0x403, 0x18000, 0x00000000 },
  { 0x0, 0x403, 0x18000, 0x00000000 }
},

#if defined(FILIB_EXTENDED)
{
  { 0x0, 0x7ff, 0x80000, 0x00000000 },
  { 0x0, 0x7ff, 0x80000, 0x00000000 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x7ff, 0x80000, 0x00000000 },
  { 0x0, 0x7ff, 0x80000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x0, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x4a5, 0x11b0e, 0xc57e649a },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

#endif
};


unsigned int srTestSet2[][2][4] = {

{
  { 0x0, 0x442, 0x04356, 0x1a882930 },
  { 0x0, 0x445, 0xb1ae4, 0xd6e2ef50 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x000, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x401, 0x40000, 0x00000000 },
  { 0x0, 0x402, 0x40000, 0x00000000 }
},

{
  { 0x1, 0x404, 0x50000, 0x00000000 },
  { 0x1, 0x400, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x400, 0x00000, 0x00000000 },
  { 0x0, 0x401, 0x20000, 0x00000000 }
},

{
  { 0x1, 0x403, 0x40000, 0x00000000 },
  { 0x1, 0x3ff, 0xccccc, 0xcccccccd }
},

{
  { 0x0, 0x400, 0xb3333, 0x33333333 },
  { 0x0, 0x400, 0xc0000, 0x00000000 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x000, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x403, 0x18000, 0x00000000 },
  { 0x0, 0x404, 0x50000, 0x00000000 }
},

{
  { 0x1, 0x402, 0x40000, 0x00000000 },
  { 0x0, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x405, 0x90000, 0x00000000 },
  { 0x0, 0x405, 0x90000, 0x00000000 }
},

{
  { 0x1, 0x3fb, 0x99999, 0x9999999a },
  { 0x0, 0x402, 0xa6666, 0x66666666 }
},

{
  { 0x1, 0x3ff, 0x00000, 0x00000000 },
  { 0x0, 0x400, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x000, 0x00000, 0x00000000 },
  { 0x0, 0x000, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x3ff, 0x19999, 0x9999999a },
  { 0x1, 0x3ff, 0x19999, 0x9999999a }
},

{
  { 0x0, 0x403, 0x18000, 0x00000000 },
  { 0x0, 0x403, 0x18000, 0x00000000 }
},

#if defined(FILIB_EXTENDED)
{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x7ff, 0x80000, 0x00000000 },
  { 0x0, 0x7ff, 0x80000, 0x00000000 }
},

{
  { 0x0, 0x7ff, 0x80000, 0x00000000 },
  { 0x0, 0x7ff, 0x80000, 0x00000000 }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x1, 0x7fe, 0xfffff, 0xffffffff },
  { 0x1, 0x3ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x1, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x1, 0x7fe, 0xfffff, 0xffffffff }
},

{
  { 0x0, 0x7fe, 0xfffff, 0xffffffff },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

{
  { 0x1, 0x7ff, 0x00000, 0x00000000 },
  { 0x0, 0x7ff, 0x00000, 0x00000000 }
},

#endif
};


