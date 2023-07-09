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
These functions create a Convolution object with the specified impulse response, partitioning, and latency. These values are static and cannot be changed after the object has been created. After creation, the object can then be used to process input audio data.

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
* latency: the desired latecncy of the algorithm, in number of samples. This will also be the block size used in the convolution. Decreasing the latency increases the amount of time it takes the algorithm to process a given amount of input. This value must be a power of 2.
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
* latency: the desired latecncy of the algorithm, in number of samples. This will also be the block size used in the head segment of the convolution. Decreasing the latency increases the amount of time it takes the algorithm to process a given amount of input. This value must be a power of 2.
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
This function creates a new convolution object with non-uniform partitioning. You specify how many segments, how many blocks per segment, and the block sizes. Not all partitionings are valid, and it is unlikely that you would choose a good one just by guessing, so in general you wouldn't use this unless you know what you are doing. Usually you would use conv_new_optimum_partitioning() instead.

Args:
* ir: a buffer containing the impulse response. It may be of any length. dft_sample_t is float. Only Mono IRs are supported.
* ir_len: the number of samples in the impulse response.
* num_segments: the number of segments in the partitioning
* block_sizes: an array of length num_segments, specifying the block sizes for each segment. The block sizes must all be powers of 2 and must increase in size. The first block size specifies the latency of the algorithm in samples.
* blocks_per_segment: an array of length num_segments specifying the number of blocks in each segment. Note that a block_size of length N must not occur eariler than N samples into the impulse response. e.g. One block of size 256 followed by a block of size 512 is not valid. 2 blocks of 256 followed by one of 512 is.
* output_callback: a function that you supply that will be called each time the algorithm has finished processing a bit of the input data.
* callback_self: a pointer to anything, or NULL, that will be passed back as the first argument to output_callback.

Returns: A fully initalized Convolution object, or NULL on failure.

##### New Optimum Partitoning
```c
Convolution*      conv_new_optimum_partitioning  (dft_sample_t*          ir, 
                                                  int                    ir_len, 
                                                  int                    latency, 
                                                  conv_output_callback_t output_callback, 
                                                  void*                  callback_self);


```
This implements the Viterbi algorithm in the op cit García paper to find the optimum partitioning. It will then return a Convolution object with that partitoning. For low latencies and long IR, the convolution obtained by this method will be much faster than double FDL or uniform block size. Be aware, however, that the Viterbi algorithm that finds the optimum partitioning has high space and time complexity. It will use 2N^2 bytes of space and has a time complexity of O(N^2) where N is ceil(ir_len / latency); N should not be greater than 2^16. In some cases it might make more sense to pre-compute vaules and recall them from a list. In the García paper, the author reports data on partitionings where N=131072 (IR length 131072, Latency 1), which even with some smart use of space requires some 68 GB of memory to compute according to their algorithm, so I'm not sure how it was computed in the year 2002. 

Args:
* ir: a buffer containing the impulse response. It may be of any length. dft_sample_t is float. Only Mono IRs are supported.
* ir_len: the number of samples in the impulse response.
* latency: the desired latecncy of the algorithm, in number of samples. This will also be the block size used in the head segment of the convolution. Decreasing the latency increases the amount of time it takes the algorithm to process a given amount of input. This value must be a power of 2.
* output_callback: a function that you supply that will be called each time the algorithm has finished processing a bit of the input data.
* callback_self: a pointer to anything, or NULL, that will be passed back as the first argument to output_callback.

Returns: A fully initalized Convolution object, or NULL on failure.


##### Convolution Destory
```c
Convolution*                                      conv_destroy          (Convolution* self);

```
Destroy the specified convolution object and free all of the memory associated with it

Args:
* self: The Convolution object to be destroyed

Returns: NULL. Call like conv = conv_destroy(conv); to ensure that your object isn't pointing at freed memory.

## Processing Audio
##### Overview
Use the previously-created Convolution object to process input data and perform the desired convolution.

##### Process
```c
void              conv_process                    (Convolution*  self, 
                                                   dft_sample_t* input, 
                                                   int           num_samples);
                                                   
```
Process audio data, perform the desired convolution.

Args:
* self: A fully initalized Convolution object created with conv_new_optimum_partitioning or similar.
* input: a float buffer containing the audio data that should be convolved with the previously specified impulse response.
* num_samples: the number of samples in input. The input buffer can be of any size. It usually makes sense for the input to be the length of the latency specified when creating the object, but it dosen't need to be. It can be longer or shorter. Consequently one call to this function could result in the completion of zero, one or many buffers of convolved output data. Therefore the convolution is not performed 'in-place', and this function does not return any output samples. Instead, the callback specified when creating the Convolution object will be called whenever output data is ready, which could be zero, one or several times when this function is called.
* callback_self: a pointer to anything, or NULL, that will be passed back as the first argument to output_callback.

Returns: void

## Getting Information
##### Overview
Get information about a Convolution object. This is especially useful if the object was created with an optimum partitioning such that the exact details are not known to the caller.

##### Get Latency
```c
int               conv_get_latency                (Convolution* self);
                                                   
```
Get the latency, in samples, of the algorithm

Args:
* self: A fully initalized Convolution object created with conv_new_optimum_partitioning or similar.

Returns: the latency of the specified object

##### Get IR Length
```c
int               ir_length           (Convolution* self);
                                                   
```
Get the length in samples of the impulse response. Note that it might be longer than you initially thought, because the end of it might have been zero-padded.

Args:
* self: A fully initalized Convolution object created with conv_new_optimum_partitioning or similar.

Returns: the number of samples in the (possibly) zero-padded impulse response

##### Get Number of Segments
```c
int               conv_get_num_segments           (Convolution* self);
                                                   
```
Get the number of segments in the partitioning

Args:
* self: A fully initalized Convolution object created with conv_new_optimum_partitioning or similar.

Returns: the number of segments of the specified object, e.g. 1 for uniform partitioning, 2 for a double_FDL, and potentially more for an optimum non-uniform partitioning.


##### Get Block Sizes
```c
void              conv_get_block_sizes           (Convolution* self, 
                                                  int*         block_sizes, 
                                                  int*         blocks_per_segment);
                                                   
```
Get the power-of-two block sizes, and number of blocks per segment, for all of the segments in the partitioning.

Args:
* self: A fully initalized Convolution object created with conv_new_optimum_partitioning or similar.
* block_sizes: pointer to an empty array of length conv_get_num_segments(self), where the requested block sizes will be written
* blocks_per_segment: pointer to an empty array of length conv_get_num_segments(self), where the requested number of blocks per segment will be written

Returns: void

##### Print Partitioning
```c
void              conv_print_partitioning        (Convolution* self);
                                                   
```
Convenience function for debugging that prints information about the partitioning to stdout. This function may change or be removed without warning and it isn't recommended to be used in production code.

Args:
* self: A fully initalized Convolution object created with conv_new_optimum_partitioning or similar.

Returns: void

