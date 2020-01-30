/* mt19937_64_5.c (20200129) */

/*
   Copyright (C) 2020, Takuji Nishimura, All rights reserved.                 

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/* 
   This code implements the pseudorandom number generator proposed in the article [1].
   
   The length of the period of this generator is 2^19937-1.
   The number of nonzero coefficients in the characteristic polynomial is 5795.
   The recurrence relation of the generator is one that holds between five words of the sequence.
   
   [1] T. Nishimura, Tables of 64-bit Mersenne Twisters,
       ACM Transactions on Medeling and Computer Simulation, Vol. 10 (2000),
       Pages 348--357.
*/

#define NN 312
#define M0 63
#define M1 151
#define M2 224

#define MATRIX_A 0xB3815B624FC82E2FULL

#define UMASK 0xFFFFFFFF80000000ULL
#define LMASK 0x7FFFFFFFULL

#define MASK_B 0x599CFCBFCA660000ULL
#define MASK_C 0xFFFAAFFE00000000ULL
#define UU 26
#define SS 17
#define TT 33
#define LL 39

static unsigned long long mt[NN];
static int mti = NN+1;

void sgenrand(unsigned long long seed)
{
    unsigned long long ux, lx;
    for (mti=0; mti<NN; mti++) {
      seed = 9797719289936477ULL * seed + 1234567ULL;
      ux = seed & 0xFFFFFFFF00000000ULL;
      seed = 9797719289936477ULL * seed + 1234567ULL;
      lx = seed >> 32;
      mt[mti] = ux | lx;
      mt[mti] += (unsigned long long)(789 * (mti+1));
    }
    mti=0;
}

static void forward_state(void)
{
    int i;
    unsigned long long x;
    unsigned long long mag01[2] = {0ULL, MATRIX_A};

    if (mti==NN+1)
      sgenrand(987654321ULL);

    for (i=0; i<NN-M2; i++) {
      x = (mt[i]&UMASK) | (mt[i+1]&LMASK);
      mt[i] = (x >> 1) ^ mag01[(int)(x&1ULL)];
      mt[i] ^= mt[i+M0] ^ mt[i+M1] ^ mt[i+M2];
    }
    for (; i<NN-M1; i++) {
      x = (mt[i]&UMASK) | (mt[i+1]&LMASK);
      mt[i] = (x >> 1) ^ mag01[(int)(x&1ULL)];
      mt[i] ^= mt[i+M0] ^ mt[i+M1] ^ mt[i+M2-NN];
    }
    for (; i<NN-M0; i++) {
      x = (mt[i]&UMASK) | (mt[i+1]&LMASK);
      mt[i] = (x >> 1) ^ mag01[(int)(x&1ULL)];
      mt[i] ^= mt[i+M0] ^ mt[i+M1-NN] ^ mt[i+M2-NN];
    }
    for (; i<NN-1; i++) {
      x = (mt[i]&UMASK) | (mt[i+1]&LMASK);
      mt[i] = (x >> 1) ^ mag01[(int)(x&1ULL)];
      mt[i] ^= mt[i+M0-NN] ^ mt[i+M1-NN] ^ mt[i+M2-NN];
    }
    x = (mt[NN-1]&UMASK) | (mt[0]&LMASK);
    mt[NN-1] = (x >> 1) ^ mag01[(int)(x&1ULL)];
    mt[NN-1] ^= mt[M0-1] ^ mt[M1-1] ^ mt[M2-1];

    mti=0;
}

// generates a random number on [0,2^64-1]-interval
unsigned long long genrand_uint64(void)
{
  unsigned long long x;
  if (mti >= NN)
    forward_state();
  x = mt[mti++];
  x ^= (x >> UU);
  x ^= (x << SS) & MASK_B;
  x ^= (x << TT) & MASK_C;
  x ^= (x >> LL);
  return x;
}

// generates a random number on [0,1]-interval
double genrand_real1(void)
{
  unsigned long long x;
  if (mti >= NN)
    forward_state();
  x = mt[mti++];
  x ^= (x >> UU);
  x ^= (x << SS) & MASK_B;
  x ^= (x << TT) & MASK_C;
  x ^= (x >> LL);
  return (5.421010862427522170e-20) * (double)x;
}

// generates a random number on [0,1)-interval
double genrand_real2(void)
{
  unsigned long long x;
  if (mti >= NN)
    forward_state();
  x = mt[mti++];
  x ^= (x >> UU);
  x ^= (x << SS) & MASK_B;
  x ^= (x << TT) & MASK_C;
  x ^= (x >> LL);
  return (5.421010862427521568e-20) * (double)x;
}

// generates a random number on (0,1)-interval
// IEEE754
double genrand_real3(void)
{
  unsigned long long x;
  if (mti >= NN)
    forward_state();
  x = mt[mti++];
  x ^= (x >> UU);
  x ^= (x << SS) & MASK_B;
  x ^= (x << TT) & MASK_C;
  x ^= (x >> LL);
  return (2.220446049250313081e-16) * ((x>>12) + 0.5);
}
