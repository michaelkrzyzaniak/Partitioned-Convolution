//gcc *.c ../common/*.c ../../lib/*.c -O2

#include "../../lib/Convolution.h"
#include "../common/MKAiff.h"
#include "../common/Timestamp.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//2048
#define BLOCK_SIZE 256
#define INPUT_SIZE BLOCK_SIZE

void output_callback(void* SELF, dft_sample_t* output, int num_samples);
void process_entire_input(Convolution* conv, MKAiff* input_aiff, MKAiff* output_aiff, dft_sample_t* in_buffer, int out_duration);

/*--------------------------------------------------------------------*/
int main(void)
{
  MKAiff* ir_aiff = aiffWithContentsOfFile("../audio/ir_mono_44_1k.wav");
  MKAiff* in_aiff = aiffWithContentsOfFile("../audio/input.wav");
  
  Convolution* conv;
  
  if(ir_aiff == NULL)
    {fprintf(stderr, "unable to open ir aiff\r\n"); exit(-1);}
  if(in_aiff == NULL)
    {fprintf(stderr, "unable to open input aiff\r\n"); exit(-1);}
  
  aiffChangeGain(ir_aiff, 0.25);
  aiffChangeGain(in_aiff, 0.25);
  
  int ir_rate = aiffSampleRate(ir_aiff);
  int in_rate = aiffSampleRate(in_aiff);
  int ir_channels = aiffNumChannels(ir_aiff);
  int in_channels = aiffNumChannels(in_aiff);
  int ir_duration = aiffDurationInFrames(ir_aiff); //131072 value used in Garcia paper
  int in_duration = aiffDurationInFrames(in_aiff);
  int out_duration = ir_duration + in_duration + BLOCK_SIZE;
  
  if(ir_channels > 1)
    {
      fprintf(stderr, "Only mono IRs are currently supported. Mixing %i channels down to 1", ir_channels);
      aiffMakeMono(ir_aiff);
      ir_channels = 1;
    }

  if(in_channels > 1)
    {
      fprintf(stderr, "Only mono input is currently supported. Mixing %i channels down to 1", in_channels);
      aiffMakeMono(in_aiff);
      in_channels = 1;
    }

  if(ir_rate != in_rate)
    {
      fprintf(stderr, "Sample rate mismatch. Input: %i, IR: %i. Resampling IR to %i\r\n", in_rate, ir_rate, in_rate);
      aiffResample(ir_aiff, in_rate, aiffInterpLinear);
      ir_rate = in_rate;
    }
  
  MKAiff* output_aiff;// = aiffWithDurationInFrames(1, in_rate, 24, out_duration + in_rate);
  //if(output_aiff == NULL)
    //{fprintf(stderr, "unable to create output aiff\r\n"); exit(-1);}

  dft_sample_t* ir_buffer = calloc(ir_duration, sizeof(*ir_buffer));
  if(ir_buffer == NULL)
    {fprintf(stderr, "unable to create ir buffer\r\n"); exit(-1);}
    
  dft_sample_t* in_buffer = calloc(INPUT_SIZE, sizeof(*ir_buffer));
  if(ir_buffer == NULL)
    {fprintf(stderr, "unable to create in buffer\r\n"); exit(-1);}

  aiffRewindPlayheadToBeginning(ir_aiff);
  aiffReadFloatingPointSamplesAtPlayhead (ir_aiff, ir_buffer, ir_duration, aiffYes);
  ir_aiff = aiffDestroy(ir_aiff);
  
  fprintf(stderr, "Single FDL:\r\n");
  output_aiff = aiffWithDurationInFrames(1, in_rate, 24, out_duration + in_rate);
  conv = conv_new_uniform_partitioning(ir_buffer, ir_duration, BLOCK_SIZE, output_callback, output_aiff);
  process_entire_input(conv, in_aiff, output_aiff, in_buffer, out_duration);
  aiffSaveWaveWithFilename(output_aiff, "single_fdl.wav");
  output_aiff = aiffDestroy(output_aiff);
  conv = conv_destroy(conv);
  
  fprintf(stderr, "Optimum Double FDL:\r\n");
  output_aiff = aiffWithDurationInFrames(1, in_rate, 24, out_duration + in_rate);
  conv = conv_new_double_fdl(ir_buffer, ir_duration, BLOCK_SIZE, output_callback, output_aiff);
  process_entire_input(conv, in_aiff, output_aiff, in_buffer, out_duration);
  aiffSaveWaveWithFilename(output_aiff, "double_fdl.wav");
  output_aiff = aiffDestroy(output_aiff);
  conv = conv_destroy(conv);
  
  fprintf(stderr, "Optimum Multiple FDL:\r\n");
  output_aiff = aiffWithDurationInFrames(1, in_rate, 24, out_duration + in_rate);
  conv = conv_new_optimum_partitioning(ir_buffer, ir_duration, BLOCK_SIZE, output_callback, output_aiff);
  process_entire_input(conv, in_aiff, output_aiff, in_buffer, out_duration);
  aiffSaveWaveWithFilename(output_aiff, "multiple_fdl.wav");
  output_aiff = aiffDestroy(output_aiff);
  conv = conv_destroy(conv);
  
  /*
  const int num_segments = 3;
  int block_sizes[num_segments] = {256, 2048, 16384};
  int blocks_per_segment[num_segments] = {8, 7, 7};
  Convolution* conv = conv_new_multiple_partitioning(ir_buffer, ir_duration, num_segments, block_sizes, blocks_per_segment, output_callback, output_aiff);
  */
  
  
  free(ir_buffer);
  in_aiff = aiffDestroy(in_aiff);
  output_aiff = aiffDestroy(output_aiff);
  free(in_buffer);
}

/*--------------------------------------------------------------------*/
void process_entire_input(Convolution* conv, MKAiff* in_aiff, MKAiff* output_aiff, dft_sample_t* in_buffer, int out_duration)
{
  aiffRewindPlayheadToBeginning(in_aiff);
   
  timestamp_microsecs_t start, end;
  start = timestamp_get_current_time();
  
  while(aiffDurationInFrames(output_aiff) < out_duration)
    {
      memset(in_buffer, 0, INPUT_SIZE * sizeof(dft_sample_t));
      aiffReadFloatingPointSamplesAtPlayhead (in_aiff, in_buffer, INPUT_SIZE, aiffYes);
      conv_process     (conv, in_buffer, INPUT_SIZE);
    }
  
  end = timestamp_get_current_time();
  
  conv_print_partitioning(conv);
  float secs = aiffDurationInSeconds(output_aiff);
  float ms   = (end-start)/1000.0;
  fprintf(stderr, "Processed %.1f seconds of audio in %.1f milliseconds (%.1f x realtime)\r\n\n", secs, ms, secs*1000.0 / ms);
}

/*--------------------------------------------------------------------*/
void output_callback(void* SELF, dft_sample_t* output, int num_samples)
{
  MKAiff* output_aiff = (MKAiff*)SELF;
  aiffAppendFloatingPointSamples(output_aiff, output, num_samples, aiffFloatSampleType);
}
