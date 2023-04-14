#ifndef __CONVOLUTION__
#define __CONVOLUTION__ 1

#if defined(__cplusplus)
extern "C"{
#endif   //(__cplusplus)

#include "DFT.h"

typedef struct OpaqueConvolutionStruct Convolution;

typedef void (*conv_output_callback_t)(void* SELF, dft_sample_t* output, int num_samples);

Convolution* conv_new_uniform_partitioning (dft_sample_t* ir, int ir_len, int latency /*power of 2*/, conv_output_callback_t output_callback, void* callback_self);

Convolution* conv_new_double_fdl           (dft_sample_t* ir, int ir_len, int latency /*power of 2*/, conv_output_callback_t output_callback, void* callback_self);

Convolution* conv_new_multiple_partitioning(dft_sample_t* ir, int ir_len, int num_segments, int block_sizes[], int blocks_per_segment[], conv_output_callback_t output_callback, void* callback_self);

Convolution* conv_new_optimum_partitioning (dft_sample_t* ir, int ir_len, int latency /*power of 2*/, conv_output_callback_t output_callback, void* callback_self);

void         conv_process           (Convolution* self, dft_sample_t* input, int num_samples);
int          conv_get_latency       (Convolution* self);
int          conv_get_num_segments  (Convolution* self);
int          conv_get_ir_length     (Convolution* self);
 /*user makes sure theses arrays are of length conv_get_num_segments()*/
void         conv_get_block_sizes   (Convolution* self, int* block_sizes, int* blocks_per_segment);
void         conv_print_partitioning(Convolution* self);

Convolution* conv_destroy           (Convolution* self);

#if defined(__cplusplus)
}
#endif   //(__cplusplus)

#endif   // __CONVOLUTION__
