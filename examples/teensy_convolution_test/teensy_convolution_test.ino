//gcc *.c ../common/*.c ../../lib/*.c -O2

#include "Convolution.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define BLOCK_SIZE 256
#define INPUT_SIZE BLOCK_SIZE
#define IR_LENGTH  (22050)
#define OUTPUT_LENGTH (44100)

unsigned output_length_so_far = 0;

void output_callback(void* SELF, dft_sample_t* output, int num_samples);

void fill_buffer(dft_sample_t *r1, int buff_size)
{
  int i;
  
  for(i=0; i<buff_size; i++)
    r1[i] = 2 * (random() / (double)RAND_MAX) - 1;
}

void setup(void)
{
  Serial.begin(1000000);

  dft_sample_t* ir_buffer = (dft_sample_t*)calloc(IR_LENGTH, sizeof(*ir_buffer));
  if(ir_buffer == NULL)
    {Serial.printf("unable to create ir buffer\r\n"); return;}
    
  dft_sample_t* in_buffer = (dft_sample_t*)calloc(INPUT_SIZE, sizeof(*in_buffer));
  if(ir_buffer == NULL)
    {Serial.printf("unable to create in buffer\r\n"); return;}

  fill_buffer(ir_buffer, IR_LENGTH);

  Convolution* conv = conv_new_uniform_partitioning(ir_buffer, IR_LENGTH, BLOCK_SIZE, output_callback, NULL);
  //Convolution* conv = conv_new_double_fdl(ir_buffer, IR_LENGTH, BLOCK_SIZE, output_callback, NULL);
  //Convolution* conv = conv_new_optimum_partitioning(ir_buffer, IR_LENGTH, BLOCK_SIZE, output_callback, NULL);

  /*
  const int num_segments = 3;
  int block_sizes[num_segments] = {256, 2048, 16384};
  int blocks_per_segment[num_segments] = {8, 7, 7};
  Convolution* conv = conv_new_multiple_partitioning(ir_buffer, IR_LENGTH, num_segments, block_sizes, blocks_per_segment, output_callback, NULL);
  */

  free(ir_buffer);
  
  uint32_t start, end;
  start = micros();
  
  while(output_length_so_far < OUTPUT_LENGTH)
    {
      fill_buffer (in_buffer, INPUT_SIZE);
      conv_process(conv, in_buffer, INPUT_SIZE);
    }
  
  end = micros();
  
  float secs = output_length_so_far / 44100.0;
  float ms   = (end-start)/1000.0;
  Serial.printf("Processed %.1f seconds of audio in %.1f milliseconds (%.1f x realtime)\r\n", secs, ms, secs*1000.0 / ms);
  

  free(in_buffer);
  conv = conv_destroy(conv);
}

void loop(void)
{
  
}
void output_callback(void* SELF, dft_sample_t* output, int num_samples)
{
  output_length_so_far += num_samples;
}
