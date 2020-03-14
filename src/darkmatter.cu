/*
A program to calculate a statistic for galaxy distributions
*/
#include <stdint.h>
#include <stdio.h>
#include "file_io.h"
#include "types.h"
#include "angles.h"
#include "cuda_macros.h"
#include <cuda_profiler_api.h>
#include <nvToolsExt.h>

#define SET_DIM  100000
#define SAME_DIM 4999950000

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

#define BIN_DIM 360

#define ANGLE_PER_THREAD 5

__global__
void calc_DR_bins(
	cart_coordinate *D,
	cart_coordinate *R,
	BIN_TYPE *bins
)
{
	uint64_t t_id_in_block = threadIdx.y*blockDim.x+threadIdx.x;
	uint64_t x_angle_start = blockDim.x*blockIdx.x*ANGLE_PER_THREAD;
	uint64_t y_angle_start = blockDim.y*blockIdx.y*ANGLE_PER_THREAD;
	
	__shared__ uint32_t block_bins[BIN_DIM];
	__shared__ cart_coordinate shD[ANGLE_PER_THREAD*32];
	__shared__ cart_coordinate shR[ANGLE_PER_THREAD*32];
	if (threadIdx.y == 31){
		for (int i = threadIdx.x;i < 32*ANGLE_PER_THREAD && x_angle_start +i < 100000;i+=32){
			shD[i] = D[x_angle_start+i];
		}
	}else if (threadIdx.y == 30){
		for (int i = threadIdx.x;i < 32*ANGLE_PER_THREAD && y_angle_start +i < 100000;i+=32){
			shR[i] = R[y_angle_start+i];
		}
	}else if (t_id_in_block < BIN_DIM){
        	block_bins[t_id_in_block] = 0;
	}
	__syncthreads();
	for(int i = threadIdx.x;i <  32*ANGLE_PER_THREAD && x_angle_start +i < 100000;i+=32){
		for(int j = threadIdx.y;j < 32*ANGLE_PER_THREAD && y_angle_start + j < 100000;j+=32){
			uint16_t bin_id = (uint16_t)(acosf(
				shD[i].x*shR[j].x+
				shD[i].y*shR[j].y+
				shD[i].z*shR[j].z
			)*RAD2QUARTDEG_CONV_RATE);	
			bin_id = MIN(bin_id,359);
			atomicAdd_block(&block_bins[bin_id],1);
		}
	}

	__syncthreads();

	if (t_id_in_block < BIN_DIM){
		atomicAdd(&bins[t_id_in_block],block_bins[t_id_in_block]);
        }
}

int main(int argc, char *argv[]) {
	if (argc < 3){printf("Too few arguments, need 2\n");return 1;}
	cudaFree(0);
	nvtxRangePush(__FUNCTION__);
	nvtxMark("Allocate");
	cart_coordinate *real_coords;
	BIN_TYPE *h_dr_bins;
	BIN_TYPE *h_dd_bins;
	BIN_TYPE *h_rr_bins;
	
	real_coords = (cart_coordinate*)calloc(100000*2,sizeof(cart_coordinate));
	h_dr_bins = (BIN_TYPE *)calloc(BIN_DIM,sizeof(BIN_TYPE));
	h_rr_bins = (BIN_TYPE *)calloc(BIN_DIM,sizeof(BIN_TYPE));
	h_dd_bins = (BIN_TYPE *)calloc(BIN_DIM,sizeof(BIN_TYPE));

	cart_coordinate *d_real_coords;
	BIN_TYPE *d_dr_bins;
	BIN_TYPE *d_dd_bins;
	BIN_TYPE *d_rr_bins;
	nvtxRangePop();

	cudaMalloc(&d_real_coords,SET_DIM*2*sizeof(cart_coordinate));
	cudaMalloc(&d_dr_bins,BIN_DIM*sizeof(BIN_TYPE));
	cudaMalloc(&d_dd_bins,BIN_DIM*sizeof(BIN_TYPE));
	cudaMalloc(&d_rr_bins,BIN_DIM*sizeof(BIN_TYPE));
	
	cart_coordinate *d_fake_coords = d_real_coords+SET_DIM;
	unsigned int grid_len =((3125)+ANGLE_PER_THREAD-1)/ANGLE_PER_THREAD;
	dim3 grid_dims(grid_len,grid_len);
	dim3 block_dims(32,32);

	cudaStream_t stream[2];
	for (int i = 0; i < 2; ++i)
		cudaStreamCreate(&stream[i]);
		
	//Read one file, start DD
	nvtxRangePush(__FUNCTION__);
	nvtxMark("Read real coords");

	galactic_coordinate *raw_coords = (galactic_coordinate *)calloc(SET_DIM,sizeof(galactic_coordinate));
	read_dataset(raw_coords,argv[1],SET_DIM);
	preprocess_cart(raw_coords,real_coords,SET_DIM);
	nvtxRangePop();

	cudaMemcpyAsync(d_real_coords,real_coords,SET_DIM*sizeof(cart_coordinate),cudaMemcpyHostToDevice,stream[0]);
	gpuErrchk(cudaPeekAtLastError());
 	calc_DR_bins<<<grid_dims, block_dims,0,stream[0]>>>(d_real_coords,d_real_coords,d_dd_bins);
	gpuErrchk(cudaPeekAtLastError());
	
	//While the first kernel runs, read second file, start RR
	nvtxRangePush(__FUNCTION__);
	nvtxMark("Read random coords");

	read_dataset(raw_coords,argv[2],SET_DIM);
	preprocess_cart(raw_coords,real_coords+SET_DIM,SET_DIM);

	cudaMemcpyAsync(d_fake_coords,real_coords+SET_DIM,SET_DIM*sizeof(cart_coordinate),cudaMemcpyHostToDevice,stream[1]);
	gpuErrchk(cudaPeekAtLastError());
	
	nvtxRangePop();
	
	cudaStreamSynchronize(stream[1]);
 	calc_DR_bins<<<grid_dims, block_dims,0,stream[1]>>>(d_fake_coords,d_fake_coords,d_rr_bins);
	gpuErrchk(cudaPeekAtLastError());
	
	//Finally, run DR
 	calc_DR_bins<<<grid_dims, block_dims,0,stream[0]>>>(d_real_coords,d_fake_coords,d_dr_bins);
	gpuErrchk(cudaPeekAtLastError());

	cudaMemcpyAsync(h_dd_bins,d_dd_bins,BIN_DIM*sizeof(BIN_TYPE),cudaMemcpyDeviceToHost);
	cudaMemcpyAsync(h_rr_bins,d_rr_bins,BIN_DIM*sizeof(BIN_TYPE),cudaMemcpyDeviceToHost);
	cudaMemcpyAsync(h_dr_bins,d_dr_bins,BIN_DIM*sizeof(BIN_TYPE),cudaMemcpyDeviceToHost);
	
        for (int i = 0; i < 2; ++i)
                cudaStreamDestroy(stream[i]);

	double kolm_smir[BIN_DIM];
	nvtxRangePush(__FUNCTION__);
	nvtxMark("Calculate statistic");
	calc_kolmogorov(kolm_smir, h_dr_bins, h_dd_bins, h_rr_bins);
	nvtxRangePop();
	nvtxRangePush(__FUNCTION__);
	nvtxMark("Write statistic");
	write_kolmogorov("out/distribution.tsv",kolm_smir,BIN_DIM);
	nvtxRangePop();

	check_sum(h_dr_bins,10000000000);
	check_sum(h_rr_bins,10000000000);
	check_sum(h_dd_bins,10000000000);
	check_values(h_dr_bins);
	check_kolmogorov(kolm_smir);

	cudaFree(d_dr_bins);
	cudaFree(d_dd_bins);
	cudaFree(d_rr_bins);
	cudaFree(d_real_coords);
}

