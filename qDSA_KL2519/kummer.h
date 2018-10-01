#ifndef KUMMER_H_
#define KUMMER_H_
#include "vecmult.h"
#include "mult.h"
#include "gfe51.h"


#define sq_gfe51(x, y) nsq_gfe51(x, y, 1)
inline void invert(gfe51 *op, gfe51 *m);

inline void gfe4x_permute(gfe4x *op,gfe4x *ip, int hl);
inline void gfe4x_hadamardUnreduced(gfe4x *op, gfe4x *ip);
inline void gfe4x_hadamard(gfe4x *op, gfe4x *ip);
inline void mul_gfe4(gfe4x *r64, gfe4x *m, gfe4x *n);
inline void mul_gfe4_expand(gfe4x *r64, gfe4x *m, gfe4x *n);
inline void sq_gfe4(gfe4x *r64, gfe4x *m);
inline void mulconst_gfe4(gfe4x *r64, gfe4x *a, const vec *b);
inline void mulconst_gfe4Unreduced(gfe4x *r64, gfe4x *a, const vec *b);
inline void scalar_mult_fixed_base_compress_freeze(unsigned char op[32], gfe4x base, unsigned char n[32]);
inline void scalar_mult_fixed_base_decompress(gfe51 *x51, gfe51 *z51, gfe4x base, unsigned char n[32]);
inline void scalar_mult_var_base_compress_freeze(unsigned char op[32], unsigned char base_rand[64], unsigned char n[32]);
inline void scalar_mult_var_base_decompress(gfe51 *x51, gfe51 *z51, unsigned char base_rand[64], unsigned char n[32]);


inline void scalar_mult_fixed_base_compress_freeze(unsigned char op[32], gfe4x base, unsigned char n[32]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe re[4],x,z,temp,xinvz;
	gfe51	x51,z51,t51;
	
	np = base;
	gfe4_t_gfe(&np, re);
	bit = 0;
	j = 3;
	i=31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1) {i--; j=7;}
	}
	//k = 1;
	
  	for(;i>=0;i--){
		for(;j>=0;j--){
			bit = (n[i]>>j) & 1;
			gfe4x_hadamard(&np, &np);
			gfe4x_permute(&npt,&np,bit);
			mul_gfe4_expand(&np, &np, &npt);
			mulconst_gfe4Unreduced(&np, &np, &BABA);
			gfe4x_hadamard(&np, &np);
			sq_gfe4(&np, &np);
			mulconst_gfe4Unreduced(&np, &np, &abxz[bit]);
		}

		j=7;
	}
	
	VECREDUCEPARTB2519((&np)->v[0],(&np)->v[1],(&np)->v[2],(&np)->v[3],(&np)->v[4],(&np)->v[5],(&np)->v[6],(&np)->v[7],(&np)->v[8]);
	gfe4_t_gfe(&np, re);
	//x = re[0];
	//z = re[1];
	//pack51(&x,&x51);
	//pack51(&z,&z51);
	pack51(&re[0],&x51);
	pack51(&re[1],&z51);	
	invert(&t51,&z51);
	mul_gfe51(&t51,&x51,&t51);
	REDUCEPARTB2519_51(t51.v[0],t51.v[1],t51.v[2],t51.v[3],t51.v[4]);
	makeUnique(&t51,&t51);
	convert_i51toc(&t51,op);
	
}

inline void scalar_mult_fixed_base_decompress(gfe51 *x51, gfe51 *z51, gfe4x base, unsigned char n[32]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe re[4],x,z,temp,xinvz;
	gfe51 t51;
	
	np = base;
	gfe4_t_gfe(&np, re);
	bit = 0;
	j = 3;
	i=31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1) {i--; j=7;}
	}
	//k = 1;
	
  	for(;i>=0;i--){
		for(;j>=0;j--){
			bit = (n[i]>>j) & 1;
			gfe4x_hadamard(&np, &np);
			gfe4x_permute(&npt,&np,bit);
			mul_gfe4_expand(&np, &np, &npt);
			mulconst_gfe4Unreduced(&np, &np, &BABA);
			gfe4x_hadamard(&np, &np);
			sq_gfe4(&np, &np);
			mulconst_gfe4Unreduced(&np, &np, &abxz[bit]);
		}

		j=7;
	}
	
	VECREDUCEPARTB2519((&np)->v[0],(&np)->v[1],(&np)->v[2],(&np)->v[3],(&np)->v[4],(&np)->v[5],(&np)->v[6],(&np)->v[7],(&np)->v[8]);
	gfe4_t_gfe(&np, re);
	//x = re[0];
	//z = re[1];
	//pack51(&x,x51);
	//pack51(&z,z51);
	pack51(&re[0],x51);
	pack51(&re[1],z51);	

	return;	
}

inline void scalar_mult_var_base_compress_freeze(unsigned char op[32], unsigned char base_rand[64], unsigned char n[32]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe work[4],re[4],x,z,temp,xinvz;
	gfe4x temp2,pabxz[2],pzxba;
	gfe51	x51,z51,t51;
		
	convert_ctoi(&work[0],base_rand);
	convert_ctoi(&work[1],base_rand+32);
	set_base_point(pabxz,work);
	gfe4_f_gfe_part1(&np, work);
	gfe4x_hadamard(&np, &np);
	sq_gfe4(&np, &np);
	mulconst_gfe4(&np, &np, &BABA);
	gfe4x_hadamard(&np, &np);
	sq_gfe4(&np, &np);
	mulconst_gfe4(&np, &np, &ab11);
	gfe4_f_gfe_part2(&np, work);

	bit = 0;
	j = 3;
	i=31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1) {i--; j=7;}
	}
	//k = 1;
		
  	for(;i>=0;i--){
    		for(;j>=0;j--){
			bit = (n[i]>>j) & 1;
			gfe4x_hadamardUnreduced(&np, &np);
			gfe4x_permute(&npt,&np,bit);
			mul_gfe4_expand(&np, &np, &npt);
			mulconst_gfe4Unreduced(&np, &np, &BABA);
			gfe4x_hadamard(&np, &np);
			sq_gfe4(&np, &np);
			mul_gfe4(&np, &np, &pabxz[bit]);
		}
		j=7;
	}
	gfe4_t_gfe(&np, re);
	//x = re[0];
	//z = re[1];
	//pack51(&x,&x51);
	//pack51(&z,&z51);
	pack51(&re[0],&x51);
	pack51(&re[1],&z51);	
	invert(&t51,&z51);
	mul_gfe51(&t51,&x51,&t51);
	REDUCEPARTB2519_51(t51.v[0],t51.v[1],t51.v[2],t51.v[3],t51.v[4]);
	makeUnique(&t51,&t51);
	convert_i51toc(&t51,op);

	return;
}


inline void scalar_mult_var_base_decompress(gfe51 *x51, gfe51 *z51, unsigned char base_rand[64], unsigned char n[32]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe work[4],re[4],x,z,temp,xinvz;
	gfe4x temp2,pabxz[2],pzxba;
	gfe51	t51;
		
	convert_ctoi(&work[0],base_rand);
	convert_ctoi(&work[1],base_rand+32);
	set_base_point(pabxz,work);
	gfe4_f_gfe_part1(&np, work);
	gfe4x_hadamard(&np, &np);
	sq_gfe4(&np, &np);
	mulconst_gfe4(&np, &np, &BABA);
	gfe4x_hadamard(&np, &np);
	sq_gfe4(&np, &np);
	mulconst_gfe4(&np, &np, &ab11);
	gfe4_f_gfe_part2(&np, work);

	bit = 0;
	j = 3;
	i=31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1) {i--; j=7;}
	}
	//k = 1;
		
  	for(;i>=0;i--){
    		for(;j>=0;j--){
			bit = (n[i]>>j) & 1;
			gfe4x_hadamardUnreduced(&np, &np);
			gfe4x_permute(&npt,&np,bit);
			mul_gfe4_expand(&np, &np, &npt);
			mulconst_gfe4Unreduced(&np, &np, &BABA);
			gfe4x_hadamard(&np, &np);
			sq_gfe4(&np, &np);
			mul_gfe4(&np, &np, &pabxz[bit]);
		}
		j=7;
	}
	gfe4_t_gfe(&np, re);
	//x = re[0];
	//z = re[1];
	//pack51(&x,x51);
	//pack51(&z,z51);
	pack51(&re[0],x51);
	pack51(&re[1],z51);	
	
	return;
}


inline void gfe4x_permute(gfe4x *op,gfe4x *ip, int hl){
	int i,a,b;
	int hl0,hl1,hl2,hl3;
	hl0 = 4*hl; hl1 = hl0+1; hl2 = hl0+2; hl3 = hl0+3;
	vec perm_mask = _mm256_set_epi32(hl3,hl2,hl1,hl0,hl3,hl2,hl1,hl0);
	op->v[0] = _mm256_permutevar8x32_epi32(ip->v[0],perm_mask);
	op->v[1] = _mm256_permutevar8x32_epi32(ip->v[1],perm_mask);
	op->v[2] = _mm256_permutevar8x32_epi32(ip->v[2],perm_mask);
	op->v[3] = _mm256_permutevar8x32_epi32(ip->v[3],perm_mask);
	op->v[4] = _mm256_permutevar8x32_epi32(ip->v[4],perm_mask);
	op->v[5] = _mm256_permutevar8x32_epi32(ip->v[5],perm_mask);
	op->v[6] = _mm256_permutevar8x32_epi32(ip->v[6],perm_mask);
	op->v[7] = _mm256_permutevar8x32_epi32(ip->v[7],perm_mask);
	op->v[8] = _mm256_permutevar8x32_epi32(ip->v[8],perm_mask);
}
/*
inline void mul_gfe(gfe *r64, gfe *m, gfe *n){

	u64 t[17];

	MULT9(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8],
	n->v[0], n->v[1], n->v[2], n->v[3], n->v[4], n->v[5], n->v[6], n->v[7], n->v[8]);

	REDUCE2519(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8],t);

}


inline void sq_gfe(gfe *r64, gfe *m){

	u64 t[17];

	SQ9(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8]);

	REDUCE2519(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8],t);
}
*/

inline void mul_gfe4(gfe4x *r64, gfe4x *m, gfe4x *n) {

	vec t[17]; 

	VECMULT9(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8],
	n->v[0], n->v[1], n->v[2], n->v[3], n->v[4], n->v[5], n->v[6], n->v[7], n->v[8]);

	VECREDUCE2519(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8], t);
}

inline void mul_gfe4_expand(gfe4x *r64, gfe4x *m, gfe4x *n) {

	vec t[17]; 

	VECMULT9(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8],
	n->v[0], n->v[1], n->v[2], n->v[3], n->v[4], n->v[5], n->v[6], n->v[7], n->v[8]);

	VECREDUCE2519expand(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8], t);
}


inline void sq_gfe4(gfe4x *r64, gfe4x *m){

	vec t[17];

	VECSQ9(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8]);

	VECREDUCE2519(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8], t);
}


const vec hadamardoffset[9] =   {{0, 0x9FFFFFA61,0, 0x9FFFFFA61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x9FFFFFF61,0, 0x9FFFFFF61},
					{0, 0x4FFFFFF61,0, 0x4FFFFFF61}};
inline void gfe4x_hadamard(gfe4x *op, gfe4x *ip){
  	int i;
	vec t[9];
	vec temp1,temp2,temp3,temp;

	temp = ip->v[0];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[0] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[0]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[1];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[1] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[1]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[2];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[2] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[2]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[3];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[3] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[3]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[4];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[4] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[4]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[5];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[5] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[5]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[6];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[6] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[6]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[7];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[7] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[7]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[8];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[8] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[8]),_mm256_xor_si256(temp3,plusminusplusminus));

	VECREDUCEPARTB2519(op->v[0],op->v[1],op->v[2],op->v[3],op->v[4],op->v[5],op->v[6],op->v[7],op->v[8]);
	
}



const vec hadamardoffsetUnreduced[9] =   {{0, 0x1FFFFFEF,0, 0x1FFFFFEF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0x1FFFFFFF,0, 0x1FFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF}};

inline void gfe4x_hadamardUnreduced(gfe4x *op, gfe4x *ip){
  	int i;
	vec t[9];
	vec temp1,temp2,temp3,temp;

	temp = ip->v[0];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[0] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[0]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[1];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[1] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[1]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[2];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[2] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[2]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[3];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[3] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[3]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[4];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[4] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[4]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[5];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[5] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[5]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[6];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[6] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[6]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[7];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[7] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[7]),_mm256_xor_si256(temp3,plusminusplusminus));

	temp = ip->v[8];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[8] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[8]),_mm256_xor_si256(temp3,plusminusplusminus));
}



inline void mulconst_gfe4Unreduced(gfe4x *r64, gfe4x *a, const vec *b){
	
	r64->v[0] = _mm256_mul_epi32(a->v[0],*b);
	r64->v[1] = _mm256_mul_epi32(a->v[1],*b);
	r64->v[2] = _mm256_mul_epi32(a->v[2],*b);
	r64->v[3] = _mm256_mul_epi32(a->v[3],*b);
	r64->v[4] = _mm256_mul_epi32(a->v[4],*b);
	r64->v[5] = _mm256_mul_epi32(a->v[5],*b);
	r64->v[6] = _mm256_mul_epi32(a->v[6],*b);
	r64->v[7] = _mm256_mul_epi32(a->v[7],*b);
	r64->v[8] = _mm256_mul_epi32(a->v[8],*b);
}

inline void mulconst_gfe4(gfe4x *r64, gfe4x *a, const vec *b){
	
	r64->v[0] = _mm256_mul_epi32(a->v[0],*b);
	r64->v[1] = _mm256_mul_epi32(a->v[1],*b);
	r64->v[2] = _mm256_mul_epi32(a->v[2],*b);
	r64->v[3] = _mm256_mul_epi32(a->v[3],*b);
	r64->v[4] = _mm256_mul_epi32(a->v[4],*b);
	r64->v[5] = _mm256_mul_epi32(a->v[5],*b);
	r64->v[6] = _mm256_mul_epi32(a->v[6],*b);
	r64->v[7] = _mm256_mul_epi32(a->v[7],*b);
	r64->v[8] = _mm256_mul_epi32(a->v[8],*b);

	VECREDUCEPARTB2519(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8]);
}

/*
inline unsigned int gt0(gfe *m){
	u32 r=0;
	int i;
	for(i=0;i<9;i++){
		r = r | m->v[i];
	}
	
	return r;
}

inline unsigned int ugtv(gfe *m, gfe *n){
	u32 r=0;

	int i;
	for(i=8;i>=0;i--){
		if (m->v[i] > n->v[i]){
			r = 1;
			break;
		}
		else if (m->v[i] < n->v[i])			
			break;
		else continue;
	}
	
	return r;
}

const unsigned int l28=0x10000000;
inline void div2(gfe *m){
	u32 temp;

	temp = m->v[8] & 1; m->v[8] = m->v[8] >> 1; m->v[7] = m->v[7] + temp*l28;
	temp = m->v[7] & 1; m->v[7] = m->v[7] >> 1; m->v[6] = m->v[6] + temp*l28;
	temp = m->v[6] & 1; m->v[6] = m->v[6] >> 1; m->v[5] = m->v[5] + temp*l28;
	temp = m->v[5] & 1; m->v[5] = m->v[5] >> 1; m->v[4] = m->v[4] + temp*l28;
	temp = m->v[4] & 1; m->v[4] = m->v[4] >> 1; m->v[3] = m->v[3] + temp*l28;
	temp = m->v[3] & 1; m->v[3] = m->v[3] >> 1; m->v[2] = m->v[2] + temp*l28;
	temp = m->v[2] & 1; m->v[2] = m->v[2] >> 1; m->v[1] = m->v[1] + temp*l28;
	temp = m->v[1] & 1; m->v[1] = m->v[1] >> 1; m->v[0] = m->v[0] + temp*l28;
	temp = m->v[0] & 1; m->v[0] = m->v[0] >> 1; 
	
}

static const unsigned int b28=0x10000000;

inline void mul2(gfe *m){
	u32 temp;

	m->v[0] = m->v[0] << 1; temp = m->v[0] >> 28; m->v[0] = m->v[0] & 0xfffffff;
	m->v[1] = m->v[1] << 1;	m->v[1] = m->v[1] + temp; temp = m->v[1] >> 28; m->v[1] = m->v[1] & 0xfffffff;
	m->v[2] = m->v[2] << 1; m->v[2] = m->v[2] + temp; temp = m->v[2] >> 28; m->v[2] = m->v[2] & 0xfffffff;
	m->v[3] = m->v[3] << 1; m->v[3] = m->v[3] + temp; temp = m->v[3] >> 28; m->v[3] = m->v[3] & 0xfffffff;
	m->v[4] = m->v[4] << 1; m->v[4] = m->v[4] + temp; temp = m->v[4] >> 28; m->v[4] = m->v[4] & 0xfffffff;
	m->v[5] = m->v[5] << 1; m->v[5] = m->v[5] + temp; temp = m->v[5] >> 28; m->v[5] = m->v[5] & 0xfffffff;
	m->v[6] = m->v[6] << 1; m->v[6] = m->v[6] + temp; temp = m->v[6] >> 28; m->v[6] = m->v[6] & 0xfffffff;
	m->v[7] = m->v[7] << 1; m->v[7] = m->v[7] + temp; temp = m->v[7] >> 28; m->v[7] = m->v[7] & 0xfffffff;
	m->v[8] = m->v[8] << 1; m->v[8] = m->v[8] + temp;

}


inline void addp(gfe *m, gfe *p){
	u32 temp;

	m->v[8] = m->v[8] + p->v[8];
	m->v[7] = m->v[7] + p->v[7];
	m->v[6] = m->v[6] + p->v[6];
	m->v[5] = m->v[5] + p->v[5];
	m->v[4] = m->v[4] + p->v[4];
	m->v[3] = m->v[3] + p->v[3];
	m->v[2] = m->v[2] + p->v[2];
	m->v[1] = m->v[1] + p->v[1];
	m->v[0] = m->v[0] + p->v[0];

}

inline void subuv(gfe *op, gfe *m, gfe *n){
	u32 b=0;
	u32 t;
	
	if(m->v[0] >= n->v[0])
		op->v[0] = m->v[0] - n->v[0];
	else{
		op->v[0] = (m->v[0] + 0x10000000) - n->v[0];
		b = 1;
	}
	t = n->v[1]+b;
	if(m->v[1] >= t){
		op->v[1] = m->v[1] - t;
		b = 0;
	}else{
		op->v[1] = (m->v[1] + 0x10000000) - t;
		b = 1;
	}

	t = n->v[2]+b;
	if(m->v[2] >= t){
		op->v[2] = m->v[2] - t;
		b = 0;
	}else{
		op->v[2] = (m->v[2] + 0x10000000) - t;
		b = 1;
	}
	t = n->v[3]+b;
	if(m->v[3] >= t){
		op->v[3] = m->v[3] - t;
		b = 0;
	}else{
		op->v[3] = (m->v[3] + 0x10000000) - t;
		b = 1;
	}
	t = n->v[4]+b;
	if(m->v[4] >= t){
		op->v[4] = m->v[4] - t;
		b = 0;
	}else{
		op->v[4] = (m->v[4] + 0x10000000) - t;
		b = 1;
	}
	t = n->v[5]+b;
	if(m->v[5] >= t){
		op->v[5] = m->v[5] - t;
		b = 0;
	}else{
		op->v[5] = (m->v[5] + 0x10000000) - t;
		b = 1;
	}
	t = n->v[6]+b;
	if(m->v[6] >= t){
		op->v[6] = m->v[6] - t;
		b = 0;
	}else{
		op->v[6] = (m->v[6] + 0x10000000) - t;
		b = 1;
	}
	t = n->v[7]+b;
	if(m->v[7] >= t){
		op->v[7] = m->v[7] - t;
		b = 0;
	}else{
		op->v[7] = (m->v[7] + 0x10000000) - t;
		b = 1;
	}
	op->v[8] = m->v[8] - n->v[8] - b;

}

inline void adduv(gfe *op, gfe *m, gfe *n){
	u32 temp;

	op->v[8] = m->v[8] + n->v[8];
	op->v[7] = m->v[7] + n->v[7];
	op->v[6] = m->v[6] + n->v[6];
	op->v[5] = m->v[5] + n->v[5];
	op->v[4] = m->v[4] + n->v[4];
	op->v[3] = m->v[3] + n->v[3];
	op->v[2] = m->v[2] + n->v[2];
	op->v[1] = m->v[1] + n->v[1];
	op->v[0] = m->v[0] + n->v[0];

	temp = op->v[0] >> 28; op->v[1] = op->v[1] + temp;     op->v[0] = op->v[0] & 0xfffffff;
	temp = op->v[1] >> 28; op->v[2] = op->v[2] + temp;     op->v[1] = op->v[1] & 0xfffffff;
	temp = op->v[2] >> 28; op->v[3] = op->v[3] + temp;     op->v[2] = op->v[2] & 0xfffffff;
	temp = op->v[3] >> 28; op->v[4] = op->v[4] + temp;     op->v[3] = op->v[3] & 0xfffffff;
	temp = op->v[4] >> 28; op->v[5] = op->v[5] + temp;     op->v[4] = op->v[4] & 0xfffffff;
	temp = op->v[5] >> 28; op->v[6] = op->v[6] + temp;     op->v[5] = op->v[5] & 0xfffffff;
	temp = op->v[6] >> 28; op->v[7] = op->v[7] + temp;     op->v[6] = op->v[6] & 0xfffffff;
	temp = op->v[7] >> 28; op->v[8] = op->v[8] + temp;     op->v[7] = op->v[7] & 0xfffffff;

}


inline void invert(gfe *op, gfe *m){
	gfe r = {0,0,0,0,0,0,0,0,0};
	gfe s = {1,0,0,0,0,0,0,0,0};
	int i,k = 0;
	gfe a = p;
	gfe b = *m;
	
	while (gt0(&b)){
		if ((a.v[0]&1) == 0){
			div2(&a);	
			mul2(&s);
		}else if((b.v[0]&1) == 0){
			div2(&b);
			mul2(&r);
		}else if(ugtv(&a, &b)){
			subuv(&a,&a,&b);
			div2(&a);
			adduv(&r,&r,&s);
			mul2(&s);
		}else {
			subuv(&b,&b,&a);
			div2(&b);
			adduv(&s,&r,&s);
			mul2(&r);
		}k = k+1;
	}

	if (ugtv(&r,&p))
		subuv(&r,&r,&p);


	for (i=1;i<=(k-251);i++){	
		if ((r.v[0]&1) == 0){
			div2(&r);
		}else{
			addp(&r,&p);		
			div2(&r);
		}
	}
	subuv(&r,&p,&r);

	mul_gfe(op, &r, &adjust_invert);
}
*/
inline void invert(gfe51 *op, gfe51 *m){
	int i;
	gfe51 t;
	gfe51 x2, x3, x4, x5, x8, x9, x11, x_5_0, x_10_0, x_15_0, x_30_0;
	gfe51 x_31_1, x_32_2, x_62_1, x_63_2, x_124_1, x_125_2, x_248_2;
	
	/* 2  */		sq_gfe51(&x2,m);
				sq_gfe51(&x4,&x2);
	/* 5 */			mul_gfe51(&x5, &x4, m);
	/* 8  */		sq_gfe51(&t,&x4);
	/* 9  */		mul_gfe51(&x9, &t, m);
	/* 11 */		mul_gfe51(&x11, &x9, &x2);
	/* 22 */		sq_gfe51(&t,&x11);
	/* 2^5-1 */		mul_gfe51(&x_5_0, &x9, &t);

	/*2^6 - 2*/		sq_gfe51(&x_10_0,&x_5_0);
	/*2^10 - 2^5*/		nsq_gfe51(&x_10_0,&x_10_0,4); //for(i=0;i<4;i++) sq_gfe(&x_10_0,&x_10_0);
	/*2^10 - 1 */		mul_gfe51(&x_10_0, &x_10_0, &x_5_0);

	/*2^11 - 2*/		sq_gfe51(&x_15_0,&x_10_0);
	/*2^15 - 2^5*/		nsq_gfe51(&x_15_0,&x_15_0,4); //for(i=0;i<4;i++) sq_gfe(&x_15_0,&x_15_0);
	/*2^15 - 1 */		mul_gfe51(&x_15_0, &x_15_0, &x_5_0);


	/*2^16 - 2*/		sq_gfe51(&x_30_0,&x_15_0);
	/*2^30 - 2^15*/		nsq_gfe51(&x_30_0,&x_30_0,14);//for(i=0;i<14;i++) sq_gfe(&x_30_0,&x_30_0);
	/*2^30 - 1 */		mul_gfe51(&x_30_0, &x_15_0, &x_30_0);

	/*2^31 - 2*/		sq_gfe51(&x_31_1,&x_30_0);

	/*2^32 - 2^2*/		sq_gfe51(&x_32_2,&x_31_1);
	/*2^33 - 2^3*/		sq_gfe51(&t,&x_32_2);
	/*2^62 - 2^32*/		nsq_gfe51(&t,&t,29);//for(i=0;i<29;i++) sq_gfe(&t,&t);
	/*2^62 - 2^2 */		mul_gfe51(&t, &t, &x_32_2);
	/*2^62 - 2 */		mul_gfe51(&x_62_1, &t, &x2);

	/*2^63 - 2^2*/		sq_gfe51(&x_63_2,&x_62_1);
	/*2^64 - 2^3*/		sq_gfe51(&t,&x_63_2);
	/*2^124 - 2^63*/	nsq_gfe51(&t,&t,60);//for(i=0;i<60;i++) sq_gfe(&t,&t);
	/*2^124 - 2^2 */	mul_gfe51(&t, &t, &x_63_2);
	/*2^124 - 2 */		mul_gfe51(&x_124_1, &t, &x2);
	
	/*2^125 - 2^2*/		sq_gfe51(&x_125_2,&x_124_1);
	/*2^126 - 2^3*/		sq_gfe51(&t,&x_125_2);
	/*2^248 - 2^125*/	nsq_gfe51(&t,&t,122);//for(i=0;i<122;i++) sq_gfe(&t,&t);
	/*2^248 - 2^2 */	mul_gfe51(&t, &t, &x_125_2);
	/*2^248 - 2 */		mul_gfe51(&t, &t, &x2);
	
	/*2^251-2^4 */		nsq_gfe51(&t,&t,3); //for(i=0;i<3;i++) sq_gfe(&t,&t);;
	
	/*2^251 - 11*/		mul_gfe51(op, &t, &x5);

}



#endif
