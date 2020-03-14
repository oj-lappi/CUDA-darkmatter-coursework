#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>

typedef struct galactic_coordinate {
	float right_ascension;
	float declination;
} galactic_coordinate;

typedef struct cart_coordinate {
	float x;
	float y;
	float z;
} cart_coordinate;

typedef struct cart_coordinate_SoA {
        float *x;
        float *y;
        float *z;
} cart_coordinate_SoA;

typedef struct cart_coordinate_SoA32 {
        float x[32];
        float y[32];
        float z[32];
} cart_coordinate_SoA32;

typedef unsigned long long BIN_TYPE;
#endif
