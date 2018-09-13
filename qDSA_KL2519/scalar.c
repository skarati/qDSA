#include "scalar.h"

static const u32 m[31] = {0xD9, 0xFE, 0xC3, 0x8E, 0x83, 0x75, 0x4F, 0xCD,
		          0x90, 0x49, 0x33, 0x91, 0x39, 0xE4, 0xDD, 0xFD, 
			  0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 
			  0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};

static const u32 mu[32] = {0x27, 0x1,  0x3C, 0x71, 0x7C, 0x8A, 0xB0, 0x32, 
				     0x6F, 0xB6, 0xCC, 0x6E, 0xC6, 0x1B, 0x22, 0x2, 
				     0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0, 
				     0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x1};

static u32 lt(u32 a,u32 b) /* 16-bit inputs */
{
  unsigned int x = a;
  x -= (unsigned int) b; /* 0..65535: no; 4294901761..4294967295: yes */
  x >>= 31; /* 0: no; 1: yes */
  return x;
}

/* Reduce coefficients of r before calling reduce_add_sub */
static void reduce_add_sub(u16 r[31])
{
  u32 pb = 0;
  u32 b;
  u32 mask;
  int i;
  unsigned char t[31];

  for(i=0;i<31;i++) 
  {
    pb += m[i];
    b = lt(r[i],pb);
    t[i] = r[i]-pb+(b<<8);
    pb = b;
  }
  mask = b - 1;
  for(i=0;i<31;i++) 
    r[i] ^= mask & (r[i] ^ t[i]);
}

/* Reduce coefficients of x before calling barrett_reduce */
static void barrett_reduce(u16 r[31], const u32 x[62])
{
  /* See HAC, Alg. 14.42 */
  int i,j;
  u32 q2[64];
  u32 *q3 = q2 + 32;
  u32 r1[32];
  u32 r2[32];
  u32 carry;
  u32 pb = 0;
  u32 b;

  for (i = 0;i < 64;++i) q2[i] = 0;
  for (i = 0;i < 32;++i) r2[i] = 0;

  for(i=0;i<32;i++)
    for(j=0;j<32;j++)
      if(i+j >= 30) q2[i+j] += mu[i]*x[j+30];
  carry = q2[30] >> 8;
  q2[31] += carry;
  carry = q2[31] >> 8;
  q2[32] += carry;

  for(i=0;i<31;i++)r1[i] = x[i];
  for(i=0;i<31;i++)
    for(j=0;j<32;j++)
      if(i+j < 32) r2[i+j] += m[i]*q3[j];

  for(i=0;i<31;i++)
  {
    carry = r2[i] >> 8;
    r2[i+1] += carry;
    r2[i] &= 0xff;
  }

  for(i=0;i<31;i++) 
  {
    pb += r2[i];
    b = lt(r1[i],pb);
    r[i] = r1[i]-pb+(b<<8);
    pb = b;
  }

  /* XXX: Can it really happen that r<0?, See HAC, Alg 14.42, Step 3 
   * If so: Handle  it here!
   */

  reduce_add_sub(r);
  reduce_add_sub(r);
}

int  group_scalar_get31(u16 r[31], const unsigned char x[31])
{
  int i;
  for(i=0;i<31;i++)
    r[i] = x[i];
  return 0;
}

int  group_scalar_get62(u16 r[31], const unsigned char x[62])
{
  int i;
    u32 t[62];
  for(i=0;i<62;i++) t[i] = x[i];
    barrett_reduce(r, t);
  return 0;
}
/*
void group_scalar_pack(unsigned char r[GROUP_SCALAR_PACKEDBYTES], const u16 *x)
{
  int i;
  for(i=0;i<31;i++)
    r[i] = x->v[i];
}
*/

static void group_scalar_add(u16 r[31], const u16 x[31], const u16 y[31])
{
  int i, carry;
  for(i=0;i<31;i++) r[i] = x[i] + y[i];
  for(i=0;i<30;i++)
  {
    carry = r[i] >> 8;
    r[i+1] += carry;
    r[i] &= 0xff;
  }
  reduce_add_sub(r);
}

void group_scalar_sub(u16 r[31], const u16 x[31], const u16 y[31])
{
  u32 b = 0;
  u32 t;
  int i;
  u16 d[31];

  for(i=0;i<31;i++)
  {
    t = m[i] - y[i] - b;
    d[i] = t & 255;
    b = (t >> 8) & 1;
  }
  group_scalar_add(r,x,d);
}

void group_scalar_mul(u16 r[31], const u16 x[31], const u16 y[31])
{
  int i,j,carry;
  u32 t[62];
  for(i=0;i<62;i++)t[i] = 0;

  for(i=0;i<31;i++)
    for(j=0;j<31;j++)
      t[i+j] += x[i] * y[j];

  /* Reduce coefficients */
  for(i=0;i<61;i++)
  {
    carry = t[i] >> 8;
    t[i+1] += carry;
    t[i] &= 0xff;
  }

  barrett_reduce(r, t);
}

static void group_scalar_setzero(u16 r[31])
{
  int i;
  for(i=0;i<31;i++)
    r[i] = 0;
}

static void group_scalar_negate(u16 r[31], const u16 x[31])
{
  u16 t[31];
  group_scalar_setzero(t);
  group_scalar_sub(r,t,x);
}

void group_scalar_set_pos(u16 r[31])
{
    if ( r[0] & 1 ) { group_scalar_negate(r, r); }
}
