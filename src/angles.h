#ifndef ANGLES_H
#define ANGLES_H

#include <math.h>
#include "types.h"

#define GOOD_ENOUGH_PI 3.14159
#define ARCMIN_CONV_RATE 0.000291
#define RAD2DEG_CONV_RATE 57.296
#define RAD2QUARTDEG_CONV_RATE 229.184

void preprocess_cart(galactic_coordinate coords_in[],cart_coordinate coords_out[], int n_coords)
{
	int i;
	for (i = 0;i < n_coords; i++) {
		float phi 	= coords_in[i].right_ascension*ARCMIN_CONV_RATE;
		float theta 	= GOOD_ENOUGH_PI/2 - coords_in[i].declination*ARCMIN_CONV_RATE;
		
		float sin_theta = sinf(theta);

		coords_out[i].x = sin_theta*cosf(phi);
		coords_out[i].y = sin_theta*sinf(phi);
		coords_out[i].z = cosf(theta);
	}
}

void preprocess_cart_SoA(galactic_coordinate coords_in[],cart_coordinate_SoA coords_out, int n_coords)
{
        int i;
        for (i = 0;i < n_coords; i++) {
                float phi       = coords_in[i].right_ascension*ARCMIN_CONV_RATE;
                float theta     = GOOD_ENOUGH_PI/2 - coords_in[i].declination*ARCMIN_CONV_RATE;

                float sin_theta = sinf(theta);
                coords_out.x[i] = sin_theta*cosf(phi);
                coords_out.y[i] = sin_theta*sinf(phi);
                coords_out.z[i] = cosf(theta);
        }
}


void preprocess_cart_SoA32(galactic_coordinate coords_in[],cart_coordinate_SoA32 coords_out[], int n_coords)
{
        int i;
        for (i = 0;i < n_coords; i++) {
                float phi       = coords_in[i].right_ascension*ARCMIN_CONV_RATE;
                float theta     = GOOD_ENOUGH_PI/2 - coords_in[i].declination*ARCMIN_CONV_RATE;

                float sin_theta = sinf(theta);
		int warp = i/32;
		int pos = i%32;
                coords_out[warp].x[pos] = sin_theta*cosf(phi);
                coords_out[warp].y[pos] = sin_theta*sinf(phi);
                coords_out[warp].z[pos] = cosf(theta);
        }
}



void check_sum(BIN_TYPE *bins, uint64_t expect){
	uint64_t sum = 0;
	for (int i = 0;i < 360;i++){
		sum+=bins[i];
	}
	if (sum == expect){
		printf("Sum test:   SUCCESS\n");
	}else{
		printf("Sum test failed, sum was: %lld expected %lld\nOff by: %lld\n",sum,expect,expect-sum);
	}
}
void check_values(BIN_TYPE *bins){
	uint64_t expected[5] = {
		397712,
		1185699,
		1948514,
		2675650,
		3386036
	};
	int err = 0;
	for (int i = 0;i < 5 ;i++){
		if (bins[i] != expected[i]){
			printf("bin %d: %lld, expected %lld. off by: %lld\n",i,bins[i],expected[i],expected[i]-bins[i]);
			err = 1;
		}
	}
	if (err == 0){
		printf("Value test: SUCCESS\n");
	}
}

void calc_kolmogorov(double *kolm_smir, BIN_TYPE *dr, BIN_TYPE *dd, BIN_TYPE *rr){
	for (int i = 0; i < 360; i++){
		kolm_smir[i] = ((double)dd[i] - 2*(double)dr[i] + (double)rr[i])/(double)rr[i];
		if (kolm_smir[i] > 10){
			printf("Bin %d suspiciously high: (%lld - 2*(%lld) + %lld)/%lld = %lf\n",
				i,
				dd[i],dr[i],rr[i],rr[i],kolm_smir[i]);
		}
	}
}

void check_kolmogorov(double *kolm_smir){
	for (int i = 0; i < 360; i++){
		if (kolm_smir[i] > .6 || kolm_smir[i] < -.6){
			printf("Bin %d not randomly distributed: omega_i = %lf\n",i,kolm_smir[i]);
		}
	}
}

#endif
