#ifndef GROUP_SCALAR_H
#define GROUP_SCALAR_H

#define GROUP_SCALAR_PACKEDBYTES 32

//typedef struct 
//{
//  uint16_t v[32]; 
//}
//group_scalar;

int group_scalar_get32(u16 r[32], const unsigned char x[32]);
int group_scalar_get64(u16 r[32], const unsigned char x[64]);
//void group_scalar_pack(unsigned char s[GROUP_SCALAR_PACKEDBYTES], const group_scalar *r);

void group_scalar_sub(u16 r[32], const u16 x[32], const u16 y[32]);
void group_scalar_mul(u16 r[32], const u16 x[32], const u16 y[32]);
void group_scalar_set_pos(u16 r[32]);

#endif
