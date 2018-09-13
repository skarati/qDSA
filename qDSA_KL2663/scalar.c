#include "scalar.h"

static const u32 m[33] = {0x97, 0xB3, 0xEE, 0x11, 0x27, 0x30, 0xE6, 0x2A,
			  0xBA, 0x10, 0xA3, 0x8B, 0xD,  0xC7, 0xE2, 0x4A,
			  0x54, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 
			  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 
			  0x55};

static const u32 mu[34] = {0xB1, 0xAF, 0x9B, 0x5E, 0xA0, 0x4E, 0xE8, 0x7D,
			   0x74, 0x69, 0x44, 0x17, 0x86, 0x0,  0x7,  0x5E, 
			   0x9,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0, 
			   0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0, 
			   0x0, 0x3};

static u32 lt(u32 a,u32 b) /* 16-bit inputs */
{
  unsigned int x = a;
  x -= (unsigned int) b; /* 0..65535: no; 4294901761..4294967295: yes */
  x >>= 31; /* 0: no; 1: yes */
  return x;
}

/* Reduce coefficients of r before calling reduce_add_sub */
static void reduce_add_sub(u16 r[33])
{
  u32 pb = 0;
  u32 b;
  u32 mask;
  int i;
  unsigned char t[33];

  for(i=0;i<33;i++) 
  {
    pb += m[i];
    b = lt(r[i],pb);
    t[i] = r[i]-pb+(b<<8);
    pb = b;
  }
  mask = b - 1;
  for(i=0;i<33;i++) 
    r[i] ^= mask & (r[i] ^ t[i]);
}

/* Reduce coefficients of x before calling barrett_reduce */
static void barrett_reduce(u16 r[32], const u32 x[64])
{
  /* See HAC, Alg. 14.42 */
  int i,j;
  u32 q2[68];
  u32 *q3 = q2 + 34;
  u32 r1[34];
  u32 r2[34];
  u32 carry;
  u32 pb = 0;
  u32 b;

  for (i = 0;i < 68;++i) q2[i] = 0;
  for (i = 0;i < 34;++i) r2[i] = 0;

  for(i=0;i<34;i++)
    for(j=0;j<34;j++)
      if(i+j >= 32) q2[i+j] += mu[i]*x[j+32];
  carry = q2[32] >> 8;
  q2[33] += carry;
  carry = q2[33] >> 8;
  q2[34] += carry;

  for(i=0;i<34;i++)r1[i] = x[i];
  for(i=0;i<33;i++)
    for(j=0;j<34;j++)
      if(i+j < 34) r2[i+j] += m[i]*q3[j];

  for(i=0;i<33;i++)
  {
    carry = r2[i] >> 8;
    r2[i+1] += carry;
    r2[i] &= 0xff;
  }

  for(i=0;i<33;i++) 
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

int  group_scalar_get32(u16 r[33], const unsigned char x[33])
{
  int i;
  for(i=0;i<33;i++)
    r[i] = x[i];
  return 0;
}

int  group_scalar_get64(u16 r[33], const unsigned char x[66])
{
  int i;
    u32 t[66];
  for(i=0;i<66;i++) t[i] = x[i];
    barrett_reduce(r, t);
  return 0;
}
/*
void group_scalar_pack(unsigned char r[GROUP_SCALAR_PACKEDBYTES], const group_scalar *x)
{
  int i;
  for(i=0;i<32;i++)
    r[i] = x->v[i];
}
*/
static void group_scalar_add(u16 r[33], const u16 x[33], const u16 y[33])
{
  int i, carry;
  for(i=0;i<33;i++) r[i] = x[i] + y[i];
  for(i=0;i<32;i++)
  {
    carry = r[i] >> 8;
    r[i+1] += carry;
    r[i] &= 0xff;
  }
  reduce_add_sub(r);
}

void group_scalar_sub(u16 r[33], const u16 x[33], const u16 y[33])
{
  u32 b = 0;
  u32 t;
  int i;
  u16 d[33];

  for(i=0;i<33;i++)
  {
    t = m[i] - y[i] - b;
    d[i] = t & 255;
    b = (t >> 8) & 1;
  }
  group_scalar_add(r,x,d);
}

void group_scalar_mul(u16 r[33], const u16 x[33], const u16 y[33])
{
  int i,j,carry;
  u32 t[66];
  for(i=0;i<66;i++)t[i] = 0;

  for(i=0;i<33;i++)
    for(j=0;j<33;j++)
      t[i+j] += x[i] * y[j];

  /* Reduce coefficients */
  for(i=0;i<65;i++)
  {
    carry = t[i] >> 8;
    t[i+1] += carry;
    t[i] &= 0xff;
  }

  barrett_reduce(r, t);
}

static void group_scalar_setzero(u16 r[33])
{
  int i;
  for(i=0;i<33;i++)
    r[i] = 0;
}

static void group_scalar_negate(u16 r[33], const u16 x[33])
{
  u16 t[33];
  group_scalar_setzero(t);
  group_scalar_sub(r,t,x);
}

void group_scalar_set_pos(u16 r[33])
{
    if ( r[0] & 1 ) { group_scalar_negate(r, r); }
}
