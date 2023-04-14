# Partitioned-Convolution
ANSI C Implementation of Guillermo García -  Optimal Filter Partition for Efficient Convolution with Short Input/Output Delay. 

This library performs fast convolution between an infinite-length (audio) sequence and a a finite-length impulse response, including very long ones up to many seconds in duration. It has arbitrarily low latency. It is suitable for state of the art real time convolution-reverb and echo-cancellation. It supports uniform partitioning, two-stage (double FDL) partitioning, and general non-uniform partitioning. It implements the Viterbi algorighm descrbed in the cited paper for finding the optimal partitioning. Additionally it does not rely on FFTW, but instead includes a variety of FFT implementations, including the one described in Sorensen et al; Real Valued Fast Fourier Transform Algorithms (1987), which is by far the fastest on any system I have tested. This library is suitable for performing convolution on Teensy4 microcontrollers (more on that soon).

## Getting Started Code Snipit
```c
#include "Convolution.h"

#define LATENCY 256                  //audio_samples, must be power of 2
#define INPUT_BUFFER_SIZE BLOCK_SIZE //could be any length
#define IR_LENGTH  (44100 * 3)       //any length

void output_callback(void* SELF, dft_sample_t* output, int num_samples);

/*--------------------------------------------------------------------*/
int main(void)
{
  dft_sample_t* ir_buffer[IR_LENGTH];
  dft_sample_t* in_buffer[INPUT_BUFFER_SIZE]

  /* you have to fill ir_buffer with the impulse response audio samples here */

  /* create a convolver with uniform partitioning */
  //Convolution* conv = conv_new_uniform_partitioning(ir_buffer, IR_LENGTH, LATENCY, output_callback, NULL);
  
  /* create a 2-stage convolver, the size of the second stage will be chosen for you so as to be optimum */
  //Convolution* conv = conv_new_double_fdl(ir_buffer, IR_LENGTH, LATENCY, output_callback, NULL);
  
  /* create a convolver with the optimum partitioning, which may be uniform or non-uniform and could have several stages */
  Convolution* conv = conv_new_optimum_partitioning(ir_buffer, IR_LENGTH, LATENCY, output_callback, NULL);

  /* create a partitoning that you specify */
  /*
  const int num_segments = 3;
  int block_sizes[num_segments] = {256, 2048, 16384};
  int blocks_per_segment[num_segments] = {8, 7, 7};
  Convolution* conv = conv_new_multiple_partitioning(ir_buffer, IR_LENGTH, num_segments, block_sizes, blocks_per_segment, output_callback, NULL);
  */
  
  
  for(;;)
    {
      /* you have to fill in_buffer with samples here */
      conv_process(conv, in_buffer, INPUT_BUFFER_SIZE);
    }
  
  conv = conv_destroy(conv);
}


/*--------------------------------------------------------------------*/
void output_callback(void* SELF, dft_sample_t* output, int num_samples)
{
  //*output contains num_samples samples of the the convolved audio
}
```

## Constructors, Destructor
##### Overview

##### New Uniform Partitioning
```c
Convolution*      conv_new_uniform_partitioning  (dft_sample_t*          ir, 
                                                  int                    ir_len, 
                                                  int                    latency, 
                                                  conv_output_callback_t output_callback, 
                                                  void*                  callback_self);
```
This function creates a new convolution object with uniform partitioning.

Args:
* ir: a buffer containing the impulse response. It may be of any length. dft_sample_t is float. Only Mono IRs are supported.
* ir_len: the number of samples in the impulse response.
* latency: the desired latecncy of the algorithm. This will also be the block size used in the convolution. Decreasing the latency increases the amount of time it takes the algorithm to process a given amount of input. This value must be a power of 2.
* output_callback: a function that you supply that will be called each time the algorithm has finished processing a bit of the input data.
* callback_self: a pointer to anything, or NULL, that will be passed back as the first argument to output_callback.

Returns: A fully initalized Convolution object, or NULL on failure.

##### New Double Frequency Delay Line
```c
Convolution*      conv_new_double_fdl            (dft_sample_t*          ir, 
                                                  int                    ir_len, 
                                                  int                    latency, 
                                                  conv_output_callback_t output_callback, 
                                                  void*                  callback_self);


```
This function creates a new convolution object with non-uniform partitioning. It will contain a head segment whose block size is equal to latency, and a tail segment whose block size is chosen to be optimal. The number of blocks in each stage is chosen to be optimal.

Args:
* ir: a buffer containing the impulse response. It may be of any length. dft_sample_t is float. Only Mono IRs are supported.
* ir_len: the number of samples in the impulse response.
* latency: the desired latecncy of the algorithm. This will also be the block size used in the head segment of the convolution. Decreasing the latency increases the amount of time it takes the algorithm to process a given amount of input. This value must be a power of 2.
* output_callback: a function that you supply that will be called each time the algorithm has finished processing a bit of the input data.
* callback_self: a pointer to anything, or NULL, that will be passed back as the first argument to output_callback.

Returns: A fully initalized Convolution object, or NULL on failure.

##### New Multiple Partitioning with Non-Uniform Block Sizes
```c
Convolution*      conv_new_multiple_partitioning (dft_sample_t*          ir, 
                                                  int                    ir_len, 
                                                  int                    num_segments, 
                                                  int                    block_sizes[], 
                                                  int                    blocks_per_segment[], 
                                                  conv_output_callback_t output_callback, 
                                                  void*                  callback_self);


```
This function creates a new convolution object with non-uniform partitioning. It will contain a head segment whose block size is equal to latency, and a tail segment whose block size is chosen to be optimal. The number of blocks in each stage is chosen to be optimal.

Args:
* ir: a buffer containing the impulse response. It may be of any length. dft_sample_t is float. Only Mono IRs are supported.
* ir_len: the number of samples in the impulse response.
* num_segments: the number of segments in the partitioning
* block_sizes: an array of length num_segments, specifying the block sizes for each segment. The block sizes must all be powers of 2 and must increase in size
* blocks_per_segment: an array of length num_segments specifying the number of blocks in each segment. Note that a block_size of length N must not occur eariler than N samples into the impulse response. e.g. One block of size 256 followed by a block of size 512 is not valid. 2 blocks of 256 followed by one of 512 is.
* output_callback: a function that you supply that will be called each time the algorithm has finished processing a bit of the input data.
* callback_self: a pointer to anything, or NULL, that will be passed back as the first argument to output_callback.

Returns: A fully initalized Convolution object, or NULL on failure.
