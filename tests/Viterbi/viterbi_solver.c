//gcc *.c ../../lib/*.c

#include "../../lib/Convolution.h"
//#include "../../lib/DFT.h"
//#include "Matrix.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#define BUFFER_SIZE (1<<5)

void find_optimum(int latency, int ir_size);

void dummy_callback(void* SELF, dft_sample_t* output, int num_samples){}

int main(void)
{
  Convolution* conv;
  
  int latency = 256;
  int maxlen  = 131072;
  int len = latency;
  
  dft_sample_t* dummy_ir = (dft_sample_t*)calloc(maxlen, sizeof(*dummy_ir));
  if(dummy_ir == NULL)
    {perror("unable to create buffer"); return 0;}
  
  while(len <= maxlen)
    {
      conv = conv_new_optimum_partitioning (dummy_ir, len, latency, dummy_callback, NULL);
      conv_print_partitioning(conv);
      len += latency;
      conv = conv_destroy(conv);
    }
  free(dummy_ir);
/*
  int latency = 127;
  int maxlat  = 1024;
  int len = 10 * 44100;
  
  while(latency <= maxlat)
    {
      find_optimum(latency, len);
      latency *= 2;
    }
*/
}

