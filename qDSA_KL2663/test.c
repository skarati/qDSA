#include "qDSA.h"
#include "measurement.h"
#include "time.h"

int main(){
	int i, v=0;
	u8 d1[33];
        u8 d2[32];
	unsigned char pk[34];
	u8 msg[32];
	u8 R[34];
	u8 s[33];

	setup(&basenp);
	
	srand(time(NULL));
	for(int j=0;j<100;j++){
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
		//verify(R, s, pk, msg);
	//});
		v += verify(R, s, pk, msg);
		printf("\nverify = %d\n",v);
		//printf("\nverify = %d\n",verify(R, s, pk, msg));
	//printf("Total CPU cycles for Verify Min: %.2f.\n", RDTSC_clk_min);
	//printf("Total CPU cycles for Verify Median: %.2f.\n", RDTSC_clk_median);
	//printf("Total CPU cycles for Verify Max: %.2f.\n", RDTSC_clk_max);
	}
	return 0;
}

