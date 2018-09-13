//#include "shake256/keccak-tiny-unrolled.c"

#include "basics.h"
#include "kummer.h"
#include "scalar.c"
#include "hash.c"

const u16 alpha0=260;
const u16 alpha1=139;
const u16 delta=3;
gfe4x basenp;

void setup(gfe4x* basenp);
void keyGen(u8* d1, u8* d2, u8* pk);
void sign(u8* R, u8* s, u8* d1, u8* d2, u8* pk, u8 msg[32]);
int  verify(u8* R, u8* s, u8* pk, u8 msg[32]);
void add_gfe54(gfe54 *r, gfe54 *x, gfe54 *z);
void sub_gfe54(gfe54 *r, gfe54 *x, gfe54 *z);
void mul_gfe54_c(gfe54 *r, gfe54 *x, const u16 c);
void mul_gfe54_2(gfe54 *r, gfe54 *x);
void convert_KL2LE(gfe54 *x, gfe54 *z);

void setup(gfe4x* basenp){
	gfe work[4];

	convert_ctoi(&work[0],base);
	convert_ctoi(&work[1],base+34);
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
	u8 d[33];
	u8 dd[64];
	gfe4x npi;
	gfe work[4];
	u8 temp[66];
	

	for(i=0;i<32;i++) d[i] = rand()%256;
	hash(dd, d, 32);
	memcpy(d1, dd, 33);
	memcpy(d2, dd+33, 31);
	dd[31] = rand()%256;
	d1[31] = d1[31] & 85;
	d1[30] = d1[30] & 63;

	scalar_mult_fixed_base_compress_freeze(pk, basenp, d1);

	return;
}

void sign(u8* R, u8* s, u8* d1, u8* d2, u8* pk, u8 msg[32]){
	int i;
	u8 haship[64];
	u8 haship2[100];
	u8 hashop[64];
	u8 r8[33], h8[33];
	u16 r16[33], h16[33], s16[33];
	gfe4x npi;
	gfe work[4];
	u32 temp[66];
		
	memcpy(haship, d2, 32);
	memcpy(haship+32, msg, 32);
	hash(hashop, haship, 64);
	for(i=0;i<64;i++) temp[i] = hashop[i];
	temp[64] = 0; temp[65] = 0;
	barrett_reduce(r16, temp);
	for(i=0;i<33;i++) r8[i] = (u8)r16[i];
	scalar_mult_fixed_base_compress_freeze(R, basenp, r8);


	memcpy(haship2, R, 34);
	memcpy(haship2+34, pk, 34);
	memcpy(haship2+68, msg, 32);
	hash(hashop, haship2, 100);	
	for(i=0;i<64;i++) temp[i] = hashop[i];
	temp[64] = 0; temp[65] = 0;
	barrett_reduce(h16, temp);
	group_scalar_set_pos(h16);
	
	for(i=0;i<33;i++) s16[i] = d1[i];
	group_scalar_mul(s16, h16, s16);
    	group_scalar_sub(s16, r16, s16);
	for(i=0;i<33;i++) s[i] = (u8)s16[i];
	

	return;
}


int  verify(u8* R, u8* s, u8* pk, u8 msg[32]){
	int i;
	u8 haship2[100];
	u8 hashop[64];
	u8 h8[33];
	u16 h16[33];
	unsigned char base_rand[68];
	gfe54	xp54,zp54,xq54,zq54,xr54,zr54;
	gfe4x npi;
	gfe work[4];
	u32 temp[66];

	scalar_mult_fixed_base_decompress(&xp54,&zp54, basenp, s);

	memcpy(haship2, R, 34);
	memcpy(haship2+34, pk, 34);
	memcpy(haship2+68, msg, 32);
	hash(hashop, haship2, 100);	
	for(i=0;i<64;i++) temp[i] = hashop[i];
	temp[64] = 0; temp[65] = 0;
	barrett_reduce(h16, temp);
	group_scalar_set_pos(h16);
	for(i=0;i<33;i++) h8[i] = (u8)h16[i];
	
	for(i=0;i<34;i++) base_rand[i]=pk[i];
	base_rand[34] = 1;
	for(i=35;i<68;i++) base_rand[i]=0;
	scalar_mult_var_base_decompress(&xq54, &zq54, base_rand, h8);

	convert_ctoi54(&xr54, R);
	zr54.v[0]=1;zr54.v[1]=0;zr54.v[2]=0;zr54.v[3]=0;zr54.v[4]=0;

	convert_KL2LE(&xp54, &zp54);
	convert_KL2LE(&xq54, &zq54);
	convert_KL2LE(&xr54, &zr54);
		
	gfe54 bxx,bzz,bxz;
	gfe54 t1,t2,t3,t4,t5,t6,t7,t8;


	mul_gfe54(&t1, &xp54, &xq54);
	mul_gfe54(&t2, &zp54, &zq54);
	mul_gfe54(&t3, &t2, &mu54);
	mul_gfe54(&t4, &xp54, &zq54);
	mul_gfe54(&t5, &xq54, &zp54);

	sub_gfe54(&bxx,&t1,&t3);
	sq_gfe54(&bxx,&bxx);
	
	add_gfe54(&t7, &t4, &t5);
	add_gfe54(&bxz, &t1, &t3);
	mul_gfe54(&bxz, &bxz, &t7);
	mul_gfe54(&t7, &t4, &t5);
	mul_gfe54(&t7, &t7, &_2mu54_1);
	sub_gfe54(&bxz,&bxz,&t7);

	sub_gfe54(&bzz,&t4,&t5);
	sq_gfe54(&bzz,&bzz);
	
	sq_gfe54(&t1,&xr54);
	mul_gfe54(&t1, &t1, &bzz);

	mul_gfe54(&t2, &xr54, &zr54);
	mul_gfe54(&t2, &t2, &bxz);
	mul_gfe54_2(&t2, &t2);

	sq_gfe54(&t3,&zr54);
	mul_gfe54(&t3, &t3, &bxx);
	
	add_gfe54(&t1,&t1,&t3);
	sub_gfe54(&t1,&t1,&t2);

	makeUnique(&t1,&t1);

	t1.v[0]=t1.v[0]|t1.v[1];
	t1.v[0]=t1.v[0]|t1.v[2];
	t1.v[0]=t1.v[0]|t1.v[3];
	t1.v[0]=t1.v[0]|t1.v[4];

	if (t1.v[0]==0) return 0;
	else return 1;

}

void mul_gfe54_c(gfe54 *r, gfe54 *x, const u16 c){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=x->v[i]*c;

	carry = (r->v[0])>>54; r->v[0]=r->v[0] & mask54; r->v[1]=r->v[1] + carry;
	carry = (r->v[1])>>54; r->v[1]=r->v[1] & mask54; r->v[2]=r->v[2] + carry;
	carry = (r->v[2])>>54; r->v[2]=r->v[2] & mask54; r->v[3]=r->v[3] + carry;
	carry = (r->v[3])>>54; r->v[3]=r->v[3] & mask54; r->v[4]=r->v[4] + carry;
	carry = (r->v[4])>>50; r->v[4]=r->v[4] & mask50; 
	carry = carry*delta; r->v[0]=r->v[0] + carry;
	carry = (r->v[0])>>54; r->v[0]=r->v[0] & mask54; r->v[1]=r->v[1] + carry;
}

void mul_gfe54_2(gfe54 *r, gfe54 *x){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=x->v[i]<<1;
}

void add_gfe54(gfe54 *r, gfe54 *x, gfe54 *z){
	int i;
	for(i=0;i<5;i++) r->v[i]=x->v[i]+z->v[i];
}

void sub_gfe54(gfe54 *r, gfe54 *x, gfe54 *z){
	int i;
	u64 carry;

	for(i=0;i<5;i++) r->v[i]=p54.v[i]+x->v[i]-z->v[i];

	carry = (r->v[0])>>54; r->v[0]=r->v[0] & mask54; r->v[1]=r->v[1] + carry;
	carry = (r->v[1])>>54; r->v[1]=r->v[1] & mask54; r->v[2]=r->v[2] + carry;
	carry = (r->v[2])>>54; r->v[2]=r->v[2] & mask54; r->v[3]=r->v[3] + carry;
	carry = (r->v[3])>>54; r->v[3]=r->v[3] & mask54; r->v[4]=r->v[4] + carry;
	carry = (r->v[4])>>50; r->v[4]=r->v[4] & mask50; 
	carry = carry*delta; r->v[0]=r->v[0] + carry;
	carry = (r->v[0])>>54; r->v[0]=r->v[0] & mask54; r->v[1]=r->v[1] + carry;
}


void convert_KL2LE(gfe54 *x, gfe54 *z){
	gfe54 t;
	
	mul_gfe54_c(&t,x,alpha1);
	mul_gfe54_c(x,z,alpha0);
	sub_gfe54(z, x, &t);
}


