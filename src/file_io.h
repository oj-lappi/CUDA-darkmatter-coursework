#include "types.h"
#include "angles.h"
#include <stdlib.h>
#include <stdio.h>
#define EMPTY_FILE_ERROR 1
#define LINE_FORMAT_ERROR  2

int buf_size = 30;

int read_dataset_FILE(galactic_coordinate dataset[], FILE *fp, int lines)
{
	galactic_coordinate c;
	rewind(fp);
	//char *buf = (char *)calloc(buf_size, sizeof(char));
	char buf[buf_size];
	if(fgets(buf, buf_size, fp) == NULL)
                return EMPTY_FILE_ERROR;

	int n = 0;
	int i;
	for (i = 0; i < lines;){
		n = sscanf(buf,"%f %f",&c.right_ascension,&c.declination);
		if(n < 2){
			if( i != 0)
				return LINE_FORMAT_ERROR;
		}else if (n == 2){
			dataset[i] = c;
			i++;
		}else{
			return n;
		}
		if (fgets(buf, buf_size, fp) == NULL)
			break;
	}
	return 0;

}

int read_dataset(galactic_coordinate dataset[], char *filename, int lines)
{
	FILE *fp = fopen(filename, "r");
        if (fp == NULL) {
                perror("Error opening file");
                return 1;
        }
	int r = read_dataset_FILE(dataset,fp,lines);
	if (fclose(fp) != 0) {
                perror("Error closing file");
                return 1;
        }

        return r;
}

int write_kolmogorov(char *filename, double kolm_smir[], int n_bins)
{
	printf("Opening file %s for writing\n",filename);
	FILE *fp = fopen(filename,"w");
	if (fp == NULL){
		//ERROR
		printf("error opening\n");
		return 1;
	}
	printf("writing...\n");
	fputs("Bin\tCount\n", fp);
	for (int i = 0;i < n_bins;i++){
		fprintf(fp, "%-3d\t%'lf\n", i, kolm_smir[i]);
	}
	if(fclose(fp) != 0){
		//ERROR
		printf("error closing\n");
		return 1;
	}
	return 0;

}


int write_bins(char *filename, BIN_TYPE bins[], int n_bins)
{
	printf("Opening file %s for writing\n",filename);
	FILE *fp = fopen(filename,"w");
	if (fp == NULL){
		//ERROR
		printf("error opening\n");
		return 1;
	}
	printf("writing...\n");
	fputs("Bin\tCount\n", fp);
	for (int i = 0;i < n_bins;i++){
		fprintf(fp, "%-3d\t%'lld\n", i, bins[i]);
	}
	if(fclose(fp) != 0){
		//ERROR
		printf("error closing\n");
		return 1;
	}
	return 0;

}

void read_all_coords(char **argv, cart_coordinate *coords, uint32_t set_dim ){
	galactic_coordinate *raw_coords = (galactic_coordinate *)calloc(set_dim*2,sizeof(galactic_coordinate));
	read_dataset(raw_coords, argv[1], set_dim);
	read_dataset(raw_coords+set_dim, argv[2], set_dim);
	preprocess_cart(raw_coords, coords, set_dim*2);
}

void read_all_coords_SoA(char **argv, cart_coordinate_SoA coords, uint32_t set_dim ){
	galactic_coordinate *raw_coords = (galactic_coordinate *)calloc(set_dim*2,sizeof(galactic_coordinate));
	read_dataset(raw_coords, argv[1], set_dim);
	read_dataset(raw_coords+set_dim, argv[2], set_dim);
	preprocess_cart_SoA(raw_coords, coords, set_dim*2);
}

void read_all_coords_SoA32(char **argv, cart_coordinate_SoA32 *coords, uint32_t set_dim ){
	galactic_coordinate *raw_coords = (galactic_coordinate *)calloc(set_dim*2,sizeof(galactic_coordinate));
	read_dataset(raw_coords, argv[1], set_dim);
	read_dataset(raw_coords+set_dim, argv[2], set_dim);
	preprocess_cart_SoA32(raw_coords, coords, set_dim*2);
}
