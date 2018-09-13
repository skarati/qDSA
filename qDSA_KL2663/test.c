//#include "shake256/keccak-tiny-unrolled.c"

#include "qDSA.h"
#include "measurement.h"
int main(){
	int i;
	u8 d1[33],d2[32];
	unsigned char pk[34];
	u8 msg[32];
	u8 R[34];
	u8 s[33];

	setup(&basenp);
	
	//MEASURE({
		keyGen(d1, d2, pk);
	//});
	//printf("Total CPU cycles for KeyGen Min: %.2f.\n", RDTSC_clk_min);
	//printf("Total CPU cycles for KeyGen Median: %.2f.\n", RDTSC_clk_median);
	//printf("Total CPU cycles for KeyGen Max: %.2f.\n", RDTSC_clk_max);

	for(i=0;i<32;i++) msg[i] = rand()%256;
	//MEASURE({
		sign(R, s, d1, d2, pk, msg);
	//});
	//printf("Total CPU cycles for Sign Min: %.2f.\n", RDTSC_clk_min);
	//printf("Total CPU cycles for Sign Median: %.2f.\n", RDTSC_clk_median);
	//printf("Total CPU cycles for Sign Max: %.2f.\n", RDTSC_clk_max);


	//MEASURE({
		printf("\nverify=%d\n",verify(R, s, pk, msg));
	//});
	//printf("Total CPU cycles for Verify Min: %.2f.\n", RDTSC_clk_min);
	//printf("Total CPU cycles for Verify Median: %.2f.\n", RDTSC_clk_median);
	//printf("Total CPU cycles for Verify Max: %.2f.\n", RDTSC_clk_max);


	return 0;
}

