#include "basics.h"
#include "kummer.h"
#include "scalar.c"
#include "hash.c"

const u8 alpha0=82;
const u8 alpha1=77;
const u8 delta=19;
gfe4x basenp;

void setup(gfe4x* basenp);
void keyGen(u8* d1, u8* d2, u8* pk);
void sign(u8* R, u8* s, u8* d1, u8* d2, u8* pk, u8 msg[32]);
int  verify(u8* R, u8* s, u8* pk, u8 msg[32]);
void add_gfe51(gfe51 *r, gfe51 *x, gfe51 *z);
void sub_gfe51(gfe51 *r, gfe51 *x, gfe51 *z);
void mul_gfe51_c(gfe51 *r, gfe51 *x, const u8 c);
void mul_gfe51_2(gfe51 *r, gfe51 *x);
void convert_KL2LE(gfe51 *x, gfe51 *z);
void convert_KL2LET2(gfe51 *x, gfe51 *z);
void clampedScalar8(u8 r[32]);
void clampedScalar16(u16 r[32]);

void setup(gfe4x* basenp){
	gfe work[4];

	convert_ctoi(&work[0],base);
	convert_ctoi(&work[1],base+32);
	gfe4_f_gfe_part1(basenp, work);
	gfe4x_hadamard(basenp, basenp);
	sq_gfe4(basenp, basenp);
	mulconst_gfe4(basenp, basenp, &BABA);
	gfe4x_hadamard(basenp, basenp);
	sq_gfe4(basenp, basenp);
	mulconst_gfe4(basenp, basenp, &ab11);
	gfe4_f_gfe_part2(basenp, work);

}

void clampedScalar8(u8 r[32]){
	int i,j;
	u64 r64[4];
	u64 r2[4], r3[4], carry;
	__int128_t r12[4];


	r[31] = r[31] & 0x13;
	r[31] = r[31] | 0x10;

	for(i=0;i<4;i++){
		r64[i] = r[i*8];
		for(j=1;j<8;j++){
			r64[i] = r64[i] | ((u64)r[i*8+j]<<(j*8));
		}
	}

	r2[0] = r64[0]<<2; 
	r2[1] = (r64[0]>>62) | (r64[1]<<2);
	r2[2] = (r64[1]>>62) | (r64[2]<<2);
	r2[3] = (r64[2]>>62) | (r64[3]<<2);

	r3[0] = r64[0]<<3;
	r3[1] = (r64[0]>>61) | (r64[1]<<3);
	r3[2] = (r64[1]>>61) | (r64[2]<<3);
	r3[3] = (r64[2]>>61) | (r64[3]<<3);
	
	r12[0] = (__int128_t)r2[0] + r3[0];
	r12[1] = (__int128_t)r2[1] + r3[1];
	r12[2] = (__int128_t)r2[2] + r3[2];
	r12[3] = (__int128_t)r2[3] + r3[3];
	
	carry = r12[0] >> 64; r12[1] += carry; r12[0] = r12[0] & 0xffffffffffffffff;
	carry = r12[1] >> 64; r12[2] += carry; r12[1] = r12[1] & 0xffffffffffffffff;
	carry = r12[2] >> 64; r12[3] += carry; r12[2] = r12[2] & 0xffffffffffffffff;
	
	for(i=0;i<4;i++){
		for(j=0;j<8;j++){
			r[8*i+j] = (u8)(r12[i] & 0xff);
			r12[i] = r12[i] >> 8;
		}
	}
	
}

void clampedScalar16(u16 r[32]){
	int i,j;
	u64 r64[4];
	u64 r2[4], r3[4], carry;
	__int128_t r12[4];


	r[31] = r[31] & 0x13;
	r[31] = r[31] | 0x10;

	for(i=0;i<4;i++){
		r64[i] = r[i*8];
		for(j=1;j<8;j++){
			r64[i] = r64[i] | ((u64)r[i*8+j]<<(j*8));
		}
	}

	r2[0] = r64[0]<<2;
	r2[1] = (r64[0]>>62) | (r64[1]<<2);
	r2[2] = (r64[1]>>62) | (r64[2]<<2);
	r2[3] = (r64[2]>>62) | (r64[3]<<2);

	r3[0] = r64[0]<<3;
	r3[1] = (r64[0]>>61) | (r64[1]<<3);
	r3[2] = (r64[1]>>61) | (r64[2]<<3);
	r3[3] = (r64[2]>>61) | (r64[3]<<3);
	
	r12[0] = (__int128_t)r2[0] + r3[0];
	r12[1] = (__int128_t)r2[1] + r3[1];
	r12[2] = (__int128_t)r2[2] + r3[2];
	r12[3] = (__int128_t)r2[3] + r3[3];
	
	carry = r12[0] >> 64; r12[1] += carry; r12[0] = r12[0] & 0xffffffffffffffff;
	carry = r12[1] >> 64; r12[2] += carry; r12[1] = r12[1] & 0xffffffffffffffff;
	carry = r12[2] >> 64; r12[3] += carry; r12[2] = r12[2] & 0xffffffffffffffff;
	
	for(i=0;i<4;i++){
		for(j=0;j<8;j++){
			r[8*i+j] = (u8)(r12[i] & 0xff);
			r12[i] = r12[i] >> 8;
		}
	}
	
}


void keyGen(u8* d1, u8* d2, u8* pk){
	int i;
	u8 d[32];
	u8 dd[64];
	gfe4x npi;
	gfe work[4];
	u8 temp[62];

	for(i=0;i<32;i++) d[i] = rand()%256;
	hash(dd, d, 32);
	memcpy(d1, dd, 32);
	memcpy(d2, dd+32, 32);
	clampedScalar8(d1);

	scalar_mult_fixed_base_compress_freeze(pk, basenp, d1);

	return;
}

void sign(u8* R, u8* s, u8* d1, u8* d2, u8* pk, u8 msg[32]){
	int i;
	u8 haship[64];
	u8 haship2[96];
	u8 hashop[64];
	u8 r8[32], h8[32];
	u16 r16[32], h16[32], s16[32];
	gfe4x npi;
	gfe work[4];
	u32 temp[64];
		
	memcpy(haship, d2, 32);
	memcpy(haship+32, msg, 32);
	hash(hashop, haship, 64);
	for(i=0;i<64;i++) temp[i] = hashop[i];

	barrett_reduce(r16, temp);
	clampedScalar16(r16);

	for(i=0;i<32;i++) r8[i] = (u8)r16[i];

	scalar_mult_fixed_base_compress_freeze(R, basenp, r8);

	memcpy(haship2, R, 32);
	memcpy(haship2+32, pk, 32);
	memcpy(haship2+64, msg, 32);
	hash(hashop, haship2, 96);	
	for(i=0;i<64;i++) temp[i] = hashop[i];
	barrett_reduce(h16, temp);
	clampedScalar16(h16);

	for(i=0;i<32;i++) s16[i] = d1[i];
	group_scalar_mul(s16, h16, s16);

    	group_scalar_sub(s16, r16, s16);
	for(i=0;i<32;i++) s[i] = (u8)s16[i];
	
	return;
}


int  verify(u8* R, u8* s, u8* pk, u8 msg[32]){
	int i;
	u8 haship2[96];
	u8 hashop[64];
	u8 h8[32];
	u16 h16[32];
	unsigned char base_rand[64];
	gfe51	xp51,zp51,xq51,zq51,xr51,zr51;
	gfe4x npi;
	gfe work[4];
	u32 temp[64];

	scalar_mult_fixed_base_decompress(&xp51,&zp51, basenp, s);

	memcpy(haship2, R, 32);
	memcpy(haship2+32, pk, 32);
	memcpy(haship2+64, msg, 32);
	hash(hashop, haship2, 96);	
	for(i=0;i<64;i++) temp[i] = hashop[i];
	barrett_reduce(h16, temp);
	for(i=0;i<32;i++) h8[i] = (u8)h16[i];
	clampedScalar8(h8);
	
	for(i=0;i<32;i++) base_rand[i]=pk[i];
	base_rand[32] = 1;
	for(i=33;i<64;i++) base_rand[i]=0;
	scalar_mult_var_base_decompress(&xq51, &zq51, base_rand, h8);

	convert_ctoi51(&xr51, R);
	zr51.v[0]=1;zr51.v[1]=0;zr51.v[2]=0;zr51.v[3]=0;zr51.v[4]=0;
	

	convert_KL2LE(&xp51, &zp51);
	convert_KL2LE(&xq51, &zq51);
	if((s[0]&1)==0)
		convert_KL2LET2(&xr51, &zr51);
	else convert_KL2LE(&xr51, &zr51);
		
	gfe51 bxx,bzz,bxz;
	gfe51 t1,t2,t3,t4,t5,t6,t7,t8;


	mul_gfe51(&t1, &xp51, &xq51);
	mul_gfe51(&t2, &zp51, &zq51);
	mul_gfe51(&t3, &t2, &mu51);
	mul_gfe51(&t4, &xp51, &zq51);
	mul_gfe51(&t5, &xq51, &zp51);

	sub_gfe51(&bxx,&t1,&t3);
	sq_gfe51(&bxx,&bxx);
	
	add_gfe51(&t7, &t4, &t5);
	add_gfe51(&bxz, &t1, &t3);
	mul_gfe51(&bxz, &bxz, &t7);
	mul_gfe51(&t7, &t4, &t5);
	mul_gfe51(&t7, &t7, &_2mu51_1);
	sub_gfe51(&bxz,&bxz,&t7);

	sub_gfe51(&bzz,&t4,&t5);
	sq_gfe51(&bzz,&bzz);
	


	sq_gfe51(&t1,&xr51);
	mul_gfe51(&t1, &t1, &bzz);

	mul_gfe51(&t2, &xr51, &zr51);
	mul_gfe51(&t2, &t2, &bxz);
	mul_gfe51_2(&t2, &t2);

	sq_gfe51(&t3,&zr51);
	mul_gfe51(&t3, &t3, &bxx);
	
	add_gfe51(&t1,&t1,&t3);
	sub_gfe51(&t1,&t1,&t2);

	makeUnique(&t1,&t1);

	
	t1.v[0]=t1.v[0]|t1.v[1];
	t1.v[0]=t1.v[0]|t1.v[2];
	t1.v[0]=t1.v[0]|t1.v[3];
	t1.v[0]=t1.v[0]|t1.v[4];

	if (t1.v[0]==0) return 0;
	else return 1;

}

void mul_gfe51_c(gfe51 *r, gfe51 *x, const u8 c){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=x->v[i]*c;

	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
	carry = (r->v[1])>>51; r->v[1]=r->v[1] & mask51; r->v[2]=r->v[2] + carry;
	carry = (r->v[2])>>51; r->v[2]=r->v[2] & mask51; r->v[3]=r->v[3] + carry;
	carry = (r->v[3])>>51; r->v[3]=r->v[3] & mask51; r->v[4]=r->v[4] + carry;
	carry = (r->v[4])>>51; r->v[4]=r->v[4] & mask51; 
	carry = carry*delta; r->v[0]=r->v[0] + carry;
	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
}

void mul_gfe51_2(gfe51 *r, gfe51 *x){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=x->v[i]<<1;
}

void add_gfe51(gfe51 *r, gfe51 *x, gfe51 *z){
	int i;
	for(i=0;i<5;i++) r->v[i]=x->v[i]+z->v[i];
}

void sub_gfe51(gfe51 *r, gfe51 *x, gfe51 *z){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=p51.v[i]+x->v[i]-z->v[i];

	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
	carry = (r->v[1])>>51; r->v[1]=r->v[1] & mask51; r->v[2]=r->v[2] + carry;
	carry = (r->v[2])>>51; r->v[2]=r->v[2] & mask51; r->v[3]=r->v[3] + carry;
	carry = (r->v[3])>>51; r->v[3]=r->v[3] & mask51; r->v[4]=r->v[4] + carry;
	carry = (r->v[4])>>51; r->v[4]=r->v[4] & mask51; 
	carry = carry*delta; r->v[0]=r->v[0] + carry;
	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
}


void convert_KL2LE(gfe51 *x, gfe51 *z){
	gfe51 t;
	
	mul_gfe51_c(x,x,alpha0);
	mul_gfe51_c(z,z,alpha1);
	sub_gfe51(z, x, z);
}


void convert_KL2LET2(gfe51 *x, gfe51 *z){
	gfe51 t;
	
	mul_gfe51_c(&t,x,alpha1);
	mul_gfe51_c(x,z,alpha0);
	sub_gfe51(z, x, &t);
}


