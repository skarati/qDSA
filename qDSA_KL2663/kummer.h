#ifndef KUMMER_H_
#define KUMMER_H_
#include "vecmult.h"
#include "mult.h"
#include "gfe54.h"

#define sq_gfe54(x, y) nsq_gfe54(x, y, 1)

inline void invert(gfe54 *op, gfe54 *m);
inline void gfe4x_permute(gfe4x *op,gfe4x *ip, int hl);
inline void gfe4x_hadamard(gfe4x *op, gfe4x *ip);
inline void gfe4x_hadamardUnreduced(gfe4x *op, gfe4x *ip);
inline void mul_gfe4(gfe4x *r64, gfe4x *m, gfe4x *n);
inline void mul_gfe4_expand(gfe4x *r64, gfe4x *m, gfe4x *n);
inline void sq_gfe4(gfe4x *r64, gfe4x *m);
inline void mulconst_gfe4(gfe4x *r64, gfe4x *a, const vec *b);
inline void mulconst_gfe4Unreduced(gfe4x *r64, gfe4x *a, const vec *b);
inline void scalar_mult_fixed_base_compress_freeze(unsigned char op[34], gfe4x base, unsigned char n[33]);
inline void scalar_mult_fixed_base_decompress(gfe54 *x54, gfe54 *z54, gfe4x base, unsigned char n[33]);
inline void scalar_mult_var_base_compress_freeze(unsigned char op[34], unsigned char base_rand[68], unsigned char n[33]);
inline void scalar_mult_var_base_decompress(gfe54 *x54, gfe54 *z54, unsigned char base_rand[68], unsigned char n[33]);

inline void scalar_mult_fixed_base_compress_freeze(unsigned char op[34], gfe4x base, unsigned char n[33]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe re[4],x,z,xinvz;
	gfe54 x54,z54,temp;
	vec mulconst;
	vec db;

	np = base;
	gfe4_t_gfe(&np, re);
	bit = 0;
	i =32;
	j=7;
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
			mul_gfe4(&np, &np, &npt);
			mulconst_gfe4Unreduced(&np, &np, &BABA);
			gfe4x_hadamard(&np, &np);
			sq_gfe4(&np, &np);
			db = _mm256_set_epi32(-bit,-bit,-bit,-bit,-bit,-bit,-bit,-bit);
			mulconst = _mm256_blendv_epi8(abxz[0], abxz[1],db);
			mulconst_gfe4Unreduced(&np, &np, &mulconst);
		}

		j=7;
	}
	VECREDUCEPARTB2663((&np)->v[0],(&np)->v[1],(&np)->v[2],(&np)->v[3],(&np)->v[4],(&np)->v[5],(&np)->v[6],(&np)->v[7],(&np)->v[8],(&np)->v[9]);
	gfe4_t_gfe(&np, re);
	pack54(&re[0],&x54);
	pack54(&re[1],&z54);
	invert(&temp,&z54);
	mul_gfe54(&temp,&x54,&temp);
	REDUCEPARTB2663(temp.v[0],temp.v[1],temp.v[2],temp.v[3],temp.v[4]);
	makeUnique(&temp,&temp);
	convert_i54toc(&temp,op);

	return ;

}

inline void scalar_mult_fixed_base_decompress(gfe54 *x54, gfe54 *z54, gfe4x base, unsigned char n[33]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe re[4],x,z,xinvz;
	gfe54 temp;
	vec mulconst;
	vec db;
	
	np = base;
	gfe4_t_gfe(&np, re);

	bit = 0;
	i =32;
	j=7;
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
			mul_gfe4(&np, &np, &npt);
			mulconst_gfe4Unreduced(&np, &np, &BABA);
			gfe4x_hadamard(&np, &np);
			sq_gfe4(&np, &np);
			db = _mm256_set_epi32(-bit,-bit,-bit,-bit,-bit,-bit,-bit,-bit);
			mulconst = _mm256_blendv_epi8(abxz[0], abxz[1],db);
			mulconst_gfe4Unreduced(&np, &np, &mulconst);
		}

		j=7;
	}
	VECREDUCEPARTB2663((&np)->v[0],(&np)->v[1],(&np)->v[2],(&np)->v[3],(&np)->v[4],(&np)->v[5],(&np)->v[6],(&np)->v[7],(&np)->v[8],(&np)->v[9]);
	gfe4_t_gfe(&np, re);
	pack54(&re[0],x54);
	pack54(&re[1],z54);
	
	return ;

}



inline void scalar_mult_var_base_compress_freeze(unsigned char op[34], unsigned char base_rand[68], unsigned char n[33]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe work[4],re[4],x,z,xinvz;
	gfe4x temp2,pabxz[2],pzxba;
	gfe54 x54,z54,temp;
	gfe4x mulconst;
	vec db;
		
	convert_ctoi(&work[0],base_rand);
	convert_ctoi(&work[1],base_rand+34);
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
	i =32;
	j=7;
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
			db = _mm256_set_epi32(-bit,-bit,-bit,-bit,-bit,-bit,-bit,-bit);
			for(k=0;k<10;k++) 
				mulconst.v[k] = _mm256_blendv_epi8(pabxz[0].v[k],pabxz[1].v[k],db);
			mul_gfe4(&np, &np, &mulconst);
		}
		j=7;
	}
	VECREDUCEPARTB2663((&np)->v[0],(&np)->v[1],(&np)->v[2],(&np)->v[3],(&np)->v[4],(&np)->v[5],(&np)->v[6],(&np)->v[7],(&np)->v[8],(&np)->v[9]);
	gfe4_t_gfe(&np, re);
	pack54(&re[0],&x54);
	pack54(&re[1],&z54);
	invert(&temp,&z54);
	mul_gfe54(&temp,&x54,&temp);
	REDUCEPARTB2663(temp.v[0],temp.v[1],temp.v[2],temp.v[3],temp.v[4]);
	makeUnique(&temp,&temp);
	convert_i54toc(&temp,op);

	return ;
}

inline void scalar_mult_var_base_decompress(gfe54 *x54, gfe54 *z54, unsigned char base_rand[68], unsigned char n[33]){
	int bit, i, j, k;
	gfe4x np,npt;
	gfe work[4],re[4],x,z,xinvz;
	gfe4x temp2,pabxz[2],pzxba;
	gfe54 temp;
	gfe4x mulconst;
	vec db;
		
	convert_ctoi(&work[0],base_rand);
	convert_ctoi(&work[1],base_rand+34);
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
	i =32;
	j=7;
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
			db = _mm256_set_epi32(-bit,-bit,-bit,-bit,-bit,-bit,-bit,-bit);
			for(k=0;k<10;k++) 
				mulconst.v[k] = _mm256_blendv_epi8(pabxz[0].v[k],pabxz[1].v[k],db);
			mul_gfe4(&np, &np, &mulconst);
		}
		j=7;
	}
	VECREDUCEPARTB2663((&np)->v[0],(&np)->v[1],(&np)->v[2],(&np)->v[3],(&np)->v[4],(&np)->v[5],(&np)->v[6],(&np)->v[7],(&np)->v[8],(&np)->v[9]);
	gfe4_t_gfe(&np, re);
	pack54(&re[0],x54);
	pack54(&re[1],z54);

	return ;
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
	op->v[9] = _mm256_permutevar8x32_epi32(ip->v[9],perm_mask);
}

inline void mul_gfe4(gfe4x *r64, gfe4x *m, gfe4x *n) {

	vec t[19]; 

	VECMULT10(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8], m->v[9],
	n->v[0], n->v[1], n->v[2], n->v[3], n->v[4], n->v[5], n->v[6], n->v[7], n->v[8], n->v[9]);

	VECREDUCE2663(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8],r64->v[9], t);
}

inline void mul_gfe4_expand(gfe4x *r64, gfe4x *m, gfe4x *n) {

	vec t[19]; 

	VECMULT10(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8], m->v[9],
	n->v[0], n->v[1], n->v[2], n->v[3], n->v[4], n->v[5], n->v[6], n->v[7], n->v[8], n->v[9]);

	VECREDUCE2663expand(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8],r64->v[9], t);
}

inline void sq_gfe4(gfe4x *r64, gfe4x *m){

	vec t[19];

	VECSQ10(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],
	m->v[0], m->v[1], m->v[2], m->v[3], m->v[4], m->v[5], m->v[6], m->v[7], m->v[8], m->v[9]);

	VECREDUCE2663(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8],r64->v[9], t);
}


const vec hadamardoffset[10] =   {{0, 0x7FFFFFD001,0, 0x7FFFFFD001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFFF001,0, 0x7FFFFFF001},
					{0, 0x7FFFFF001,0, 0x7FFFFF001}};

inline void gfe4x_hadamard(gfe4x *op, gfe4x *ip){
  	int i;
	vec t[10];
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

	temp = ip->v[9];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[9] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffset[9]),_mm256_xor_si256(temp3,plusminusplusminus));

	VECREDUCEPARTB2663(op->v[0],op->v[1],op->v[2],op->v[3],op->v[4],op->v[5],op->v[6],op->v[7],op->v[8],op->v[9]);
	
}

const vec hadamardoffsetUnreduced[10] =   {{0, 0xFFFFFFB,0, 0xFFFFFFB},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFFF,0, 0xFFFFFFF},
					{0, 0xFFFFFF,0, 0xFFFFFF}};

inline void gfe4x_hadamardUnreduced(gfe4x *op, gfe4x *ip){
  	int i;
	vec t[10];
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
	
	temp = ip->v[9];
	temp2 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(2,2,0,0));
	temp3 = _mm256_permute4x64_epi64 (temp,_MM_SHUFFLE(3,3,1,1));
	op->v[9] = _mm256_add_epi64(_mm256_add_epi64(temp2,hadamardoffsetUnreduced[9]),_mm256_xor_si256(temp3,plusminusplusminus));
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
	r64->v[9] = _mm256_mul_epi32(a->v[9],*b);
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
	r64->v[9] = _mm256_mul_epi32(a->v[9],*b);

	VECREDUCEPARTB2663(r64->v[0],r64->v[1],r64->v[2],r64->v[3],r64->v[4],r64->v[5],r64->v[6],r64->v[7],r64->v[8],r64->v[9]);
}

inline void invert(gfe54 *op, gfe54 *m) {
	int i;
	gfe54 x2;
	gfe54 x_2_1,x_3_2,x_4_1,x_7_3,x_8_1,x_8_4,x_16_1,x_32_1,x_64_1,x_128_1,x_256_1,x_264_1,x_266_5;
	
	/* 2  */		sq_gfe54(&x2,m);
	/*2^2-1*/		mul_gfe54(&x_2_1, &x2, m);
	/*2^3-2*/		sq_gfe54(&x_3_2,&x_2_1);
	/*2^4-2^2*/		sq_gfe54(&x_4_1,&x_3_2);
	/*2^4-1*/		mul_gfe54(&x_4_1, &x_4_1, &x_2_1);
								
	/*2^5-2*/		sq_gfe54(&x_7_3,&x_4_1);
	/*2^7-2^3*/		nsq_gfe54(&x_7_3,&x_7_3,2);
	/*2^8-2^4*/		sq_gfe54(&x_8_1,&x_7_3);
	/*2^8-1*/		mul_gfe54(&x_8_1, &x_8_1, &x_4_1);
								
								

	/*2^9-2*/		sq_gfe54(&x_16_1,&x_8_1);
	/*2^16-2^8*/		nsq_gfe54(&x_16_1,&x_16_1,7);
	/*2^16-1*/		mul_gfe54(&x_16_1, &x_8_1, &x_16_1);

	/*2^17-2*/		sq_gfe54(&x_32_1,&x_16_1);
	/*2^32-2^16*/		nsq_gfe54(&x_32_1,&x_32_1,15);
	/*2^32-1*/		mul_gfe54(&x_32_1, &x_32_1, &x_16_1);

	/*2^33-2*/		sq_gfe54(&x_64_1,&x_32_1);
	/*2^64-2^32*/		nsq_gfe54(&x_64_1,&x_64_1,31);
	/*2^64-1*/		mul_gfe54(&x_64_1, &x_32_1, &x_64_1);

	/*2^65-2*/		sq_gfe54(&x_128_1,&x_64_1);
	/*2^128-2^64*/		nsq_gfe54(&x_128_1,&x_128_1,63);
	/*2^128-1*/		mul_gfe54(&x_128_1, &x_128_1, &x_64_1);

	/*2^129-2*/		sq_gfe54(&x_256_1,&x_128_1);
	/*2^256-2^128*/		nsq_gfe54(&x_256_1,&x_256_1,127);
	/*2^256-1*/		mul_gfe54(&x_256_1, &x_128_1, &x_256_1);


	/*2^257-2*/		sq_gfe54(&x_264_1,&x_256_1);
	/*2^263-2^7*/		nsq_gfe54(&x_264_1,&x_264_1,6);

	/*2^263-2^3*/		mul_gfe54(&x_266_5, &x_264_1, &x_7_3);
	/*2^263-2*/		mul_gfe54(&x_266_5, &x_266_5, &x_3_2);
	/*2^263-1*/		mul_gfe54(&x_266_5, &x_266_5, m);
	/*2^266-2^3*/		nsq_gfe54(&x_266_5,&x_266_5,3);

	/*2^266-5*/		mul_gfe54(op, &x_266_5, &x_2_1);

}


#endif
