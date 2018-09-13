#ifndef GROUP_SCALAR_H
#define GROUP_SCALAR_H

#define GROUP_SCALAR_PACKEDBYTES 31
typedef uint16_t u16;
//typedef struct 
//{
//  uint16_t v[32]; 
//}
//group_scalar;

int group_scalar_get31(u16 r[31], const unsigned char x[31]);
int group_scalar_get62(u16 r[31], const unsigned char x[62]);
//void group_scalar_pack(unsigned char s[GROUP_SCALAR_PACKEDBYTES], const u16 r[31]);

void group_scalar_sub(u16 r[31], const u16 x[31], const u16 y[31]);
void group_scalar_mul(u16 r[31], const u16 x[31], const u16 y[31]);
void group_scalar_set_pos(u16 r[31]);

#endif
