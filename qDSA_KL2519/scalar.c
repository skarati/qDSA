#include "scalar.h"

static const u32 m[32] = {0xD9, 0xFE, 0xC3, 0x8E, 0x83, 0x75, 0x4F, 0xCD,
		          0x90, 0x49, 0x33, 0x91, 0x39, 0xE4, 0xDD, 0xFD, 
			  0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 
			  0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x00};


static const u64 m16[17] = {0xFED9, 0x8EC3, 0x7583, 0xCD4F, 0x4990, 0x9133, 0xE439, 0xFDDD, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFF,0x00};

static const u32 mu[32] = {0x27, 0x1,  0x3C, 0x71, 0x7C, 0x8A, 0xB0, 0x32, 
				     0x6F, 0xB6, 0xCC, 0x6E, 0xC6, 0x1B, 0x22, 0x2, 
				     0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0, 
				     0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x1};

static const u32 mu16[17] = {0x48C, 0x127, 0x713C, 0x8A7C, 0x32B0, 0xB66F, 0x6ECC, 0x1BC6, 0x222, 0x0, 0x0,0x0, 0x0, 0x0, 0x0, 0x0, 0x100};

static u64 lt(u64 a,u64 b) /* 16-bit inputs */
{
  u64 x = a;
  x -= (u64) b; /* 0..65535: no; 4294901761..4294967295: yes */
  x >>= 63; /* 0: no; 1: yes */
  return x;
}

/* Reduce coefficients of r before calling reduce_add_sub */
static void reduce_add_sub(u16 r[32])
{
  u32 pb = 0;
  u32 b;
  u32 mask;
  int i;
  unsigned char t[32];

  for(i=0;i<32;i++) 
  {
    pb += m[i];
    b = lt(r[i],pb);
    t[i] = r[i]-pb+(b<<8);
    pb = b;
  }
  mask = b - 1;
  for(i=0;i<32;i++) 
    r[i] ^= mask & (r[i] ^ t[i]);
}

static void reduce_add_sub16(u64 r[17])
{
  u64 pb = 0;
  u64 b;
  u64 mask;
  int i;
  u64 t[17];

  for(i=0;i<17;i++) 
  {
    pb += m[i];
    b = lt(r[i],pb);
    t[i] = r[i]-pb+(b<<16);
    pb = b;
  }
  mask = b - 1;
  for(i=0;i<17;i++) 
    r[i] ^= mask & (r[i] ^ t[i]);
}



/* Reduce coefficients of x before calling barrett_reduce */

static void barrett_reduce16(u16 r[32], const u32 x[64]){
	int i,j;

	u64 r16[17];
	u64 x16[32];
	
	for(i=0;i<32;i++){
		x16[i] = x[2*i] | ((u64)x[2*i+1] << 8);
	}

  	u64 q2[34];
  	u64 *q3 = q2 + 17;
  	u64 r1[17];
  	u64 r2[17];
  	u64 carry;
  	u64 pb = 0;
  	u64 b;

  	for (i = 0;i < 34;++i) 	q2[i] = 0;
  	for (i = 0;i < 17;++i) 	r2[i] = 0;

	for(i=0;i<17;i++)
    		for(j=0;j<17;j++)
      			if(i+j >= 15) q2[i+j] += mu16[i]*x16[j+15];

  	carry = q2[14] >> 16;
  	q2[15] += carry;
  	carry = q2[15] >> 16;
  	q2[16] += carry;
  	carry = q2[16] >> 16;
  	q2[17] += carry;

	for(i=0;i<17;i++)r1[i] = x16[i];
  	for(i=0;i<16;i++)
    		for(j=0;j<17;j++)
      			if(i+j < 17) r2[i+j] += m16[i]*q3[j];


  	for(i=0;i<16;i++){
    		carry = r2[i] >> 16;
    		r2[i+1] += carry;
    		r2[i] &= 0xffff;
  	}
	r2[i] &= 0xffff;

  	for(i=0;i<17;i++){
    		pb += r2[i];
    		b = lt(r1[i],pb);
    		r16[i] = r1[i]-pb+(b<<16);
    		pb = b;
 	}
	

	i=0;
	while(r16[16]>0){
  		reduce_add_sub16(r16);
		i++;
	}
	reduce_add_sub16(r16);
	
	for(i=0;i<16;i++){
		r[2*i] = r16[i] & 0xff;
		r[2*i+1] = r16[i] >> 8;
	}

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
    barrett_reduce16(r, t);
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

static void group_scalar_add(u16 r[32], const u16 x[32], const u16 y[32])
{
  int i, carry;
  for(i=0;i<32;i++) r[i] = x[i] + y[i];
  for(i=0;i<31;i++)
  {
    carry = r[i] >> 8;
    r[i+1] += carry;
    r[i] &= 0xff;
  }
  reduce_add_sub(r);
}

void group_scalar_sub(u16 r[32], const u16 x[32], const u16 y[32])
{
  u32 b = 0;
  u32 t=0;
  int i;
  u16 d[32];

  for(i=0;i<32;i++)
  {
    t += y[i];
    b = lt(x[i],t);
    r[i] = x[i]-t+(b<<8);
    t = b;
  }
 
}

void group_scalar_mul(u16 r[32], const u16 x[32], const u16 y[32])
{
  int i,j,carry;
  u32 t[64];
  for(i=0;i<64;i++)t[i] = 0;

  for(i=0;i<32;i++)
    for(j=0;j<32;j++)
      t[i+j] += x[i] * y[j];

  /* Reduce coefficients */
  for(i=0;i<63;i++)
  {
    carry = t[i] >> 8;
    t[i+1] += carry;
    t[i] &= 0xff;
  }

  barrett_reduce16(r, t);
}

static void group_scalar_setzero(u16 r[32])
{
  int i;
  for(i=0;i<32;i++)
    r[i] = 0;
}

static void group_scalar_negate(u16 r[32], const u16 x[32])
{
  u16 t[32];
  group_scalar_setzero(t);
  group_scalar_sub(r,t,x);
}

void group_scalar_set_pos(u16 r[32])
{
    if ( r[0] & 1 ) { group_scalar_negate(r, r); }
}
