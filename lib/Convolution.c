#include "Convolution.h"

#include <string.h> //memcpy
#include <stdlib.h> //calloc
#include <math.h>   //ceil

#include <stdio.h> //TESTING ONLY


//beware double evaluation
#define CONV_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define CONV_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define CONV_K 1.5

/*--------------------------------------------------------------------*/
typedef struct OpaqueUniformPartitioningStruct
{
  int filter_length;
  int block_length;
  int num_blocks;
  int fft_n;
  
  dft_sample_t* input_buffer;
  dft_sample_t* output_buffer;
  dft_sample_t** filter_blocks;
  dft_sample_t** delay_line;
  
  int delay_line_index;
  int num_samples_received;
  
  int segment_offset; //for mixing the output with other segments

  void* parent;
}Convolution_Uniform_Partitioning;

/*--------------------------------------------------------------------*/
struct OpaqueConvolutionStruct
{
  Convolution_Uniform_Partitioning** filter_segments;
  int num_filter_segments;
  conv_output_callback_t output_ready_callback;
  void* output_ready_callback_self;
  dft_sample_t* output_buffer;
  int output_buffer_length;
  int output_ring_index;
  float cost;
};

/*--------------------------------------------------------------------*/
Convolution_Uniform_Partitioning* conv_new_uniform_partitioning_private(dft_sample_t* ir, int ir_len, int block_size /*power of 2*/, int segment_offset,  void* parent);
Convolution_Uniform_Partitioning* conv_destroy_uniform_paritioning(Convolution_Uniform_Partitioning* self);
void conv_process_uniform_partitioning(Convolution_Uniform_Partitioning* self, dft_sample_t* input, int num_samples);
void conv_output_ready(Convolution* SELF, dft_sample_t* output, int num_samples, int segment_offset);
//void conv_calculate_total_filter_length(Convolution* self);

/*--------------------------------------------------------------------*/
Convolution_Uniform_Partitioning* conv_new_uniform_partitioning_private(dft_sample_t* ir, int ir_len, int block_size /*power of 2*/, int segment_offset,  void* parent)
{
  Convolution_Uniform_Partitioning* self = (Convolution_Uniform_Partitioning*) calloc(1, sizeof(*self));
  if(self != NULL)
    {
      int i;
      
      /* scan the tail of the ir for silence here? */
      /* normalize the ir here? */
      
      self->parent = parent;
      self->segment_offset = segment_offset;
      
      self->block_length = dft_smallest_power_of_2_at_least_as_great_as(block_size);
      self->num_blocks = ceil(ir_len / block_size);
      if(ir_len % block_size)
        ++self->num_blocks;
      self->filter_length = self->num_blocks * self->block_length;
      self->fft_n = 2 * self->block_length;
    
      self->input_buffer  = (dft_sample_t*)calloc(self->block_length, sizeof(dft_sample_t));
      self->output_buffer = (dft_sample_t*)calloc(self->fft_n, sizeof(dft_sample_t));
      if((self->input_buffer == NULL) || (self->output_buffer == NULL))
        return conv_destroy_uniform_paritioning(self);
    
      self->filter_blocks = (dft_sample_t**)calloc(self->num_blocks, sizeof(dft_sample_t*));
      self->delay_line    = (dft_sample_t**)calloc(self->num_blocks, sizeof(dft_sample_t*));
      if((self->filter_blocks == NULL) || (self->delay_line == NULL))
          return conv_destroy_uniform_paritioning(self);

      /*
        break ir into blocks of length self->block_length.
        0-pad to size self->fft_n (2*self->block_length)
        then take the DFT of each
      */
      for(i=0; i<self->num_blocks; i++)
        {
          self->filter_blocks[i] = (dft_sample_t*)calloc(self->fft_n, sizeof(dft_sample_t));
          self->delay_line[i]    = (dft_sample_t*)calloc(self->fft_n, sizeof(dft_sample_t));
          if((self->filter_blocks[i] == NULL) || (self->delay_line[i] == NULL))
            return conv_destroy_uniform_paritioning(self);
          
          int num_samples = CONV_MIN(ir_len, self->block_length);
          memcpy(self->filter_blocks[i], ir, num_samples*sizeof(dft_sample_t));
          
          rdft_real_forward_dft(self->filter_blocks[i], self->fft_n);
          //rdft_bit_reverse_indices(self->filter_blocks[i], self->fft_n);
          
          ir_len -= num_samples;
          ir += num_samples;
        }
    }
  
  return self;
}

/*--------------------------------------------------------------------*/
Convolution_Uniform_Partitioning* conv_destroy_uniform_paritioning(Convolution_Uniform_Partitioning* self)
{
  if(self != NULL)
    {
      int i;
      if(self->input_buffer != NULL)
        free(self->input_buffer);
      if(self->output_buffer != NULL)
        free(self->output_buffer);

      if(self->filter_blocks != NULL)
        {
          for(i=0; i<self->num_blocks; i++)
            if(self->filter_blocks[i] != NULL)
              free(self->filter_blocks[i]);
     
          free(self->filter_blocks);
        }
     
      if(self->delay_line != NULL)
        {
          for(i=0; i<self->num_blocks; i++)
            if(self->delay_line[i] != NULL)
              free(self->delay_line[i]);
     
          free(self->delay_line);
        }
     
      free(self);
    }
  return (Convolution_Uniform_Partitioning*) NULL;
}

/*--------------------------------------------------------------------*/
void conv_process_uniform_partitioning(Convolution_Uniform_Partitioning* self, dft_sample_t* input, int num_samples)
{
  while(num_samples > 0)
    {
      int samples_needed = self->block_length - self->num_samples_received;
      int samples_used = CONV_MIN(samples_needed, num_samples);
      
      memcpy(self->input_buffer+self->num_samples_received, input, samples_used*sizeof(dft_sample_t));
      self->num_samples_received += samples_used;
    
      if(self->num_samples_received == self->block_length)
        {
          int i, block;
          
          dft_sample_t* real = self->delay_line[self->delay_line_index];
          memcpy(real+self->block_length, self->input_buffer, self->block_length*sizeof(dft_sample_t));
        
          rdft_real_forward_dft(real, self->fft_n);
          memset(self->output_buffer, 0, self->fft_n * sizeof(dft_sample_t));
        
          for(block=0; block<self->num_blocks; block++)
            {
              dft_sample_t* a = self->delay_line[(self->delay_line_index + block) % self->num_blocks];
              dft_sample_t* b = self->filter_blocks[block];
              dft_sample_t* c = self->output_buffer;
              int N = self->fft_n;
              
              c[0] += a[0] * b[0];
              c[self->block_length] += a[self->block_length] * b[self->block_length];
              
              for(i=1; i<self->block_length; i++)
                {
                  c[i]   += a[i]*b[i]   - a[N-i]*b[N-i];
                  c[N-i] += a[i]*b[N-i] + b[i]*a[N-i];
                }
            }
        
          rdft_real_inverse_dft(self->output_buffer, self->fft_n);
          //self->output_ready_callback(self->output_ready_callback_self, self->output_buffer+self->block_length, self->block_length);

          conv_output_ready(self->parent, self->output_buffer+self->block_length, self->block_length, self->segment_offset);
        
          self->num_samples_received = 0;
          //++self->delay_line_index;
          self->delay_line_index += self->num_blocks - 1;
          self->delay_line_index %= self->num_blocks;
          
          //store first half of of current window in sencond half of next window for overlap-save
          memcpy(self->delay_line[self->delay_line_index], self->input_buffer, self->block_length*sizeof(dft_sample_t));
        }
        
      num_samples -= samples_used;
      input += samples_used;
    }
}

/*--------------------------------------------------------------------*/
void conv_output_ready(Convolution* SELF, dft_sample_t* output, int num_samples, int segment_offset)
{
  Convolution* self = (Convolution*) SELF;
  
  if(self->num_filter_segments == 1)
    {
      self->output_ready_callback(self->output_ready_callback_self, output, num_samples);
      return;
    }
    
  else
    {
      int i;
      int start = self->output_ring_index;
      if(segment_offset > 0)
        start += segment_offset - num_samples;
      
      
      for(i=0; i<num_samples; i++)
        self->output_buffer[(start + i) % self->output_buffer_length] += output[i];
      
      if(segment_offset == 0)
        {
          self->output_ready_callback(self->output_ready_callback_self, self->output_buffer + self->output_ring_index, num_samples);
          memset(self->output_buffer + self->output_ring_index, 0, num_samples * sizeof(dft_sample_t));
          self->output_ring_index += num_samples;
          self->output_ring_index %= self->output_buffer_length;
        }
    }
}

/*--------------------------------------------------------------------*/
Convolution* conv_new_private(int num_segments, conv_output_callback_t output_callback, void* callback_self)
{
  Convolution* self = (Convolution*) calloc(1, sizeof(*self));
  if(self != NULL)
    {
      self->output_ready_callback = output_callback;
      self->output_ready_callback_self = callback_self;
      self->num_filter_segments = num_segments;
      self->filter_segments = (Convolution_Uniform_Partitioning**)calloc(self->num_filter_segments, sizeof(*self->filter_segments));
      if(self->filter_segments == NULL)
        return conv_destroy(self);
    }
  return self;
}

/*--------------------------------------------------------------------*/
Convolution* conv_new_uniform_partitioning(dft_sample_t* ir, int ir_len, int latency /*power of 2*/, conv_output_callback_t output_callback, void* callback_self)
{
  Convolution* self = conv_new_private(1, output_callback, callback_self);
  if(self != NULL)
    {
      self->filter_segments[0] = conv_new_uniform_partitioning_private(ir, ir_len, latency, 0, self);
      if(self->filter_segments[0] == NULL)
        return conv_destroy(self);
      self->cost = 4 * CONV_K * log2(2*latency) + 4*self->filter_segments[0]->num_blocks;
    }
  return self;
}

/*--------------------------------------------------------------------*/
float double_fdl_cost(float k, float N, float B, float T)
{
  return 4*k*log2(2*N) + 4*B/(double)N + 4*k*log2(2*B) + 4*((T/(float)B)-1);
}

/*--------------------------------------------------------------------*/
Convolution* conv_new_double_fdl(dft_sample_t* ir, int ir_len, int latency /*power of 2*/, conv_output_callback_t output_callback, void* callback_self)
{
  int T = latency * ceil(ir_len / (float)latency);
  float k = CONV_K;
  float n = latency;
  float kn_over_2_ln_2 = k*n/1.3863;
  float B = -kn_over_2_ln_2 + sqrt((kn_over_2_ln_2*kn_over_2_ln_2) + (T * latency));
  int b_upper = dft_smallest_power_of_2_at_least_as_great_as(B);
  int b_lower = b_upper >> 1;
  int b_upper_cost = double_fdl_cost(k, n, b_upper, T);
  int b_lower_cost = double_fdl_cost(k, n, b_lower, T);
  int b_final = (b_upper_cost < b_lower_cost) ? b_upper : b_lower;
  //if(b_final > 2048) b_final = 2048;
  
  int head_length = b_final;
  int tail_length = T - head_length;
  
  int block_sizes[2]        = {latency, b_final};
  int blocks_per_segment[2] = {head_length/latency, tail_length/b_final};
  
  Convolution* self = conv_new_multiple_partitioning(ir, ir_len, 2, block_sizes, blocks_per_segment, output_callback, callback_self);
  if(self != NULL)
    {
      self->cost = (b_upper_cost < b_lower_cost) ? b_upper_cost : b_lower_cost;
    }
  return self;
}

/*--------------------------------------------------------------------*/
Convolution* conv_new_multiple_partitioning(dft_sample_t* ir, int ir_len, int num_segments, int block_sizes[], int blocks_per_segment[], conv_output_callback_t output_callback, void* callback_self)
{
  Convolution* self = conv_new_private(num_segments, output_callback, callback_self);
  if(self != NULL)
    {
      int offset = 0;
      int filter_len = 0;
      int ir_segment_len;
      int i;
      
      for(i=0; i<num_segments; i++)
        {
          filter_len = blocks_per_segment[i] * block_sizes[i];
          ir_segment_len = CONV_MIN(filter_len, ir_len);
          if(ir_segment_len < 0) ir_segment_len = 0;
          self->filter_segments[i] = conv_new_uniform_partitioning_private(ir+offset, ir_segment_len, block_sizes[i], offset, self);
          offset += filter_len;
          ir_len -= filter_len;
          if(self->filter_segments[i] == NULL)
            return conv_destroy(self);
        }
        
      self->output_buffer_length = offset - filter_len;
      self->output_buffer = calloc(self->output_buffer_length, sizeof(dft_sample_t));
      if(self->output_buffer == NULL)
        return conv_destroy(self);
    }
  return self;
}

/*--------------------------------------------------------------------*/
int conv_get_p(int i)
{
  int p = 0;
  i += 1;
  
  while(i > 0)
    {
      i >>= 1;
      ++p;
    }
  
  return 1 << (p-1);
}

/*--------------------------------------------------------------------*/
#define IX(vals, cols, row, col) vals[(row) * (cols) + (col)]
Convolution* conv_new_optimum_partitioning(dft_sample_t* ir, int ir_len, int latency, conv_output_callback_t output_callback, void* callback_self)
{
  if (ir_len <= 0) return NULL;
  
  float k = CONV_K;
  int maxnum_segments = 32; //usually 4 is more than enough
  int t, T = ir_len / latency;
  if(ir_len % latency) ++T;
  int T_minus_1 = T - 1;
  int num_segments;
  int block_sizes[maxnum_segments];
  int blocks_per_segment[maxnum_segments];
  
  //Matrix*  cost_matrix = matrix_new(T, T);
  unsigned short* cost_matrix = (unsigned short*)calloc(T*T, sizeof(*cost_matrix));
  
  if(cost_matrix == NULL)
    return NULL;
    //{perror("unable to init matrix\r\n"); return;}
  
  int i, j, p1, q1, p2, q2;
  
  unsigned short min_cost, min_index = 0;
  
  int accumulated_cost = min_cost = 4*k*log2(2*latency) + 4;
  IX(cost_matrix, T, 0, 0) = accumulated_cost;

  accumulated_cost += 4;
  IX(cost_matrix, T, 0, 1) = accumulated_cost;
  IX(cost_matrix, T, T_minus_1-1, T_minus_1-1) = 0;

  //calculate transition costs and indices
  //put the costs in the upper traingle at index y, x,
  //put the indices in the lower triangle at index T-1-y, T-1-x
  for(t=2; t<T; t++)
    {
      p2 = q2 = 1;
      
      for(i=0; i<t; i++)
        {
          min_cost = 0xFFFF;
  
          p1 = q1 = 1;
          
          for(j=0; j<t-1; j++)
            {
              accumulated_cost = IX(cost_matrix, T, j, t-1);
              
              //subsequent fractions of the same block
              if((p1 == p2) && ((q1+1) == q2))
                accumulated_cost += 0;
                
              //append another block of the same size
              else if((p1 == p2) && (q2==1))
                accumulated_cost += 4;
                
              //start a new larger block
              else if((p1 < p2) && (q2==1))
                accumulated_cost += 4*k*log2(2*latency*p2) + 4;
              
              //illegal transition
              else
                accumulated_cost = 0xFFFF;
              
              if(accumulated_cost < min_cost)
                {
                  min_cost = accumulated_cost;
                  min_index = j;
                }
                
              ++q1; if(q1 > p1) {p1*=2; q1=1;}
            }
          
          IX(cost_matrix, T, i, t) = min_cost;
          IX(cost_matrix, T, T_minus_1-i, T_minus_1-t) = min_index;
    
          ++q2; if(q2 > p2) {p2*=2; q2=1;}
        }
    }
  
  //find min cost in the last column
  if(T > 1)
    min_cost = 0xFFFF;

  for(i=0; i<T-1; i++)
    {
      accumulated_cost = IX(cost_matrix, T, i, T_minus_1);
      if(accumulated_cost < min_cost)
        {
          min_cost = accumulated_cost;
          min_index = i;
        }
    }
  
  //trace backwards to get the min indices
  //store the winners in the first row of matrix
  //since we are done with it
  for(t=T_minus_1; t>0; t--)
    {
      IX(cost_matrix, T, 0, T_minus_1-t) = min_index;
      min_index = IX(cost_matrix, T, T_minus_1-min_index, T_minus_1-t);
    }
    

  //forward pass to get final block sizes and blocks per segment
  num_segments = 1;
  block_sizes[0] = latency;
  blocks_per_segment[0] = 1;
  
  if(T > 1)
    ++blocks_per_segment[0];
  
  p1 = q1 = 1;
  
  for(t=2; t<T; t++)
    {
      min_index = IX(cost_matrix, T, 0, T_minus_1-t);
      p2 = conv_get_p(min_index);
      q2 = min_index - p2 + 2;
      
      //subsequent fractions of the same block
      if((p1 == p2) && ((q1+1) == q2))
        {
          // do nothing.
        }
    
      //append another block of the same size
      else if((p1 == p2) && (q2==1))
        ++blocks_per_segment[num_segments-1];
    
      //start a new larger block
      else if((p1 < p2) && (q2==1))
        {
          if((num_segments + 1) < maxnum_segments)
            {
              block_sizes[num_segments] = p2 * latency;
              blocks_per_segment[num_segments] = 1;
              ++num_segments;
            }
          //else
            //fprintf(stderr, "ran out of space, increase maxnum_segments in Convolution.c");
        }
        
      p1 = p2; q1 = q2;
    }
  
  free(cost_matrix);
  
  Convolution* self = conv_new_multiple_partitioning(ir, ir_len, num_segments, block_sizes, blocks_per_segment, output_callback, callback_self);
  if(self != NULL)
    {
      self->cost = min_cost;
    }
  return self;
}

/*--------------------------------------------------------------------*/
int conv_get_latency(Convolution* self)
{
  return self->filter_segments[0]->block_length;
}

/*--------------------------------------------------------------------*/
int conv_get_num_segments(Convolution* self)
{
  return self->num_filter_segments;
}

/*--------------------------------------------------------------------*/
int conv_get_ir_length(Convolution* self)
{
  int length = 0;
  int i;
  for(i=0; i<self->num_filter_segments; i++)
    length += self->filter_segments[i]->filter_length;
    
  return length;
}

/*--------------------------------------------------------------------*/
void conv_get_block_sizes(Convolution* self, int* block_sizes, int* blocks_per_segment /*user makes sure theses arrays are of length conv_get_num_segments()*/)
{
  int i;
  for(i=0; i<self->num_filter_segments; i++)
    {
      block_sizes[i] = self->filter_segments[i]->block_length;
      blocks_per_segment[i] = self->filter_segments[i]->num_blocks;
    }
}

/*--------------------------------------------------------------------*/
void conv_print_partitioning(Convolution* self)
{
  int i;
  int ir_length = conv_get_ir_length(self);
  int latency = conv_get_latency(self);
  int num_segments = conv_get_num_segments(self);
  int T = ir_length / latency;
  int block_sizes[num_segments];
  int blocks_per_segment[num_segments];
  conv_get_block_sizes(self, block_sizes, blocks_per_segment);
  
  
  fprintf(stderr, "ir_size: %i x %i = %i, partition: [", latency, T, ir_length);
  for(i=0; i<num_segments; i++)
    {
      fprintf(stderr, "%i x %i", blocks_per_segment[i], block_sizes[i]);
      if(i < num_segments-1)
        fprintf(stderr, " --> ");
    }
  fprintf(stderr, "], cost: %i\r\n", (int)self->cost);
}

/*--------------------------------------------------------------------*/
Convolution* conv_destroy(Convolution* self)
{
  if(self != NULL)
    {
      if(self->filter_segments != NULL)
        {
          int i;
          for(i=0; i<self->num_filter_segments; i++)
            self->filter_segments[i] = conv_destroy_uniform_paritioning(self->filter_segments[i]);
          free(self->filter_segments);
        }
      
      if(self->output_buffer != NULL)
        free(self->output_buffer);
      
      free(self);
    }
    
  return (Convolution*) NULL;
}

/*--------------------------------------------------------------------*/
void conv_process(Convolution* self, dft_sample_t* input, int num_samples)
{
  int i;
  int smallest_block_size = self->filter_segments[0]->block_length;
  
  while(num_samples > 0)
    {
      int samples_to_process = CONV_MIN(num_samples, smallest_block_size);
      for(i=0; i<self->num_filter_segments; i++)
        conv_process_uniform_partitioning(self->filter_segments[i], input, samples_to_process);
      
      input += samples_to_process;
      num_samples -= samples_to_process;
    }
}
