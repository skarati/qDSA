basics.h: contains the basic macros.
mult.h: contains the macros for multiplication, squaring and reduction.
vecmult.h: contains the macros for SIMD multiplication, squaring and reduction.
kummer.h: contains the functions for the kummer ladder.
measurement.h: contains the macros for measuring the code.
main.c: the C file used to measure the performance.

The .S files are used for 5-limb multiplication and squaring and are modified versions of the corresponding Sandy2x code available from
https://bench.cr.yp.to/supercop/supercop-20160910.tar.xz.

check_curve.m: magma code for determining the various parameters of the Kummer line.
check_curve.txt: output of check_curve.m

Compile command:
gcc -m64 -mavx2 -O3 -fomit-frame-pointer main.c mul_gfe54.S nsq_gfe54.S consts.S
