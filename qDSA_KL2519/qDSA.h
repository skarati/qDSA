//#include "shake256/keccak-tiny-unrolled.c"

#include "basics.h"
#include "kummer.h"
#include "scalar.c"
#include "hash.c"
//#include "measurement.h"

const u8 alpha0=81;
const u8 alpha1=20;
const u8 delta=9;
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
	d1[31] = d1[31]&0x7;
	d1[31] = d1[31]|0x4;
	d1[0] = d1[0] & 0xF8;

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
	barrett_reduce16(r16, temp);

	r16[31] = r16[31]&0x7;
	r16[31] = r16[31]|0x4;
	r16[0] = r16[0] & 0xF8;
	for(i=0;i<32;i++) r8[i] = (u8)r16[i];

	scalar_mult_fixed_base_compress_freeze(R, basenp, r8);
	
	memcpy(haship2, R, 32);
	memcpy(haship2+32, pk, 32);
	memcpy(haship2+64, msg, 32);
	hash(hashop, haship2, 96);	
	for(i=0;i<64;i++) temp[i] = hashop[i];

	barrett_reduce16(h16, temp);
	h16[31] = h16[31]&0x7;
	h16[31] = h16[31]|0x4;
	h16[0] = h16[0] & 0xF8;

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
	gfe51	xp51,zp51,xq51,zq51,xr51,zr51,xr51_2,zr51_2;
	gfe4x npi;
	gfe work[4];
	u32 temp[64];
	u8 v1=0,v2=0;

	scalar_mult_fixed_base_decompress(&xp51,&zp51, basenp, s);

	memcpy(haship2, R, 32);
	memcpy(haship2+32, pk, 32);
	memcpy(haship2+64, msg, 32);
	hash(hashop, haship2, 96);	
	for(i=0;i<64;i++) temp[i] = hashop[i];

	barrett_reduce16(h16, temp);

	h16[31] = h16[31]&0x7;
	h16[31] = h16[31]|0x4;
	h16[0] = h16[0] & 0xF8;
	for(i=0;i<32;i++) h8[i] = (u8)h16[i];
		
	
	for(i=0;i<32;i++) base_rand[i]=pk[i];
	base_rand[32] = 1;
	for(i=33;i<64;i++) base_rand[i]=0;
	scalar_mult_var_base_decompress(&xq51, &zq51, base_rand, h8);

	convert_ctoi51(&xr51,R);
	zr51.v[0]=1;zr51.v[1]=0;zr51.v[2]=0;zr51.v[3]=0;zr51.v[4]=0;

	convert_KL2LE(&xp51, &zp51);
	convert_KL2LE(&xq51, &zq51);
	if((s[0]&1)==0)
		convert_KL2LET2(&xr51, &zr51);
	else
		convert_KL2LE(&xr51, &zr51);
	
	
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

	u64 carry;
	carry = (t1.v[0])>>51; t1.v[0]=t1.v[0] & mask51; t1.v[1]=t1.v[1] + carry;
	carry = (t1.v[1])>>51; t1.v[1]=t1.v[1] & mask51; t1.v[2]=t1.v[2] + carry;
	carry = (t1.v[2])>>51; t1.v[2]=t1.v[2] & mask51; t1.v[3]=t1.v[3] + carry;
	carry = (t1.v[3])>>51; t1.v[3]=t1.v[3] & mask51; t1.v[4]=t1.v[4] + carry;
	carry = (t1.v[4])>>47; t1.v[4]=t1.v[4] & mask47; 
	carry = carry*delta; t1.v[0]=t1.v[0] + carry;
	carry = (t1.v[0])>>51; t1.v[0]=t1.v[0] & mask51; t1.v[1]=t1.v[1] + carry;

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
	carry = (r->v[4])>>47; r->v[4]=r->v[4] & mask47; 
	carry = carry*delta; r->v[0]=r->v[0] + carry;
	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
}

void mul_gfe51_2(gfe51 *r, gfe51 *x){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=x->v[i]<<1;

	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
	carry = (r->v[1])>>51; r->v[1]=r->v[1] & mask51; r->v[2]=r->v[2] + carry;
	carry = (r->v[2])>>51; r->v[2]=r->v[2] & mask51; r->v[3]=r->v[3] + carry;
	carry = (r->v[3])>>51; r->v[3]=r->v[3] & mask51; r->v[4]=r->v[4] + carry;
	carry = (r->v[4])>>47; r->v[4]=r->v[4] & mask47; 
	carry = carry*delta; r->v[0]=r->v[0] + carry;
	carry = (r->v[0])>>51; r->v[0]=r->v[0] & mask51; r->v[1]=r->v[1] + carry;
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
	carry = (r->v[4])>>47; r->v[4]=r->v[4] & mask47; 
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


