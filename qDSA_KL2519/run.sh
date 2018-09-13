clear
gcc -m64 -mavx2 -O3 -fomit-frame-pointer test.c consts.S mul_gfe51.S nsq_gfe51.S
./a.out
rm a.out
