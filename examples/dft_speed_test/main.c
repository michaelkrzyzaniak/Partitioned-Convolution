//gcc *.c ../common/*.c ../../lib/*.c -O2

#include "../../lib/DFT.h"
#include "../common/Timestamp.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define NUM_TESTS   1000
#define INCLUDE_FORWARD_TRANSFORMS
#define INCLUDE_INVERSE_TRANSFORMS
#define SAMPLE_RATE 44100

void fill_buffers(dft_sample_t *r1, dft_sample_t *r2, dft_sample_t *i1, dft_sample_t *i2, int buff_size)
{
  //double phase_1 = 0;
  //double phase_2 = 0;
  int i;
  
  for(i=0; i<buff_size; i++)
    {
      //real_1[i] = sin(phase_1);
      //phase_1 += 1000 * 2 * M_PI / (double)SAMPLE_RATE;

      //real_2[i] = sin(phase_2);
      //phase_2 += 2000 * 2 * M_PI / (double)SAMPLE_RATE;
      
      r1[i] = 2 * (random() / (double)RAND_MAX) - 1;
      r2[i] = 2 * (random() / (double)RAND_MAX) - 1;
      i1[i] = 0;
      i2[i] = 0;
    }
}

int main()
{
  timestamp_microsecs_t start, end;
  int time_1, time_2, time_3, time_4, time_5;
  unsigned test, buff_size;
  
  fprintf(stderr, "N\t2x dft_complex_forward_dft\t2x dft_real_forward_dft\t2x rdft_real_forward_dft\t1x dft_2_real_forward_dfts\t1x rdft_2_real_forward_dfts\r\n");
  
  for(buff_size=8; buff_size<=4096; buff_size*=2)
    {
      dft_sample_t *real_1 = calloc(buff_size, sizeof(dft_sample_t));
      dft_sample_t *imag_1 = calloc(buff_size, sizeof(dft_sample_t));
      dft_sample_t *real_2 = calloc(buff_size, sizeof(dft_sample_t));
      dft_sample_t *imag_2 = calloc(buff_size, sizeof(dft_sample_t));
  
      start = timestamp_get_current_time();
      for(test=0; test<NUM_TESTS; test++)
        {
          fill_buffers(real_1, real_2, imag_1, imag_2, buff_size);
#ifdef INCLUDE_FORWARD_TRANSFORMS
          dft_complex_forward_dft(real_1, imag_1, buff_size);
          dft_complex_forward_dft(real_2, imag_2, buff_size);
#endif
#ifdef INCLUDE_INVERSE_TRANSFORMS
          dft_complex_inverse_dft(real_1, imag_1, buff_size);
          dft_complex_inverse_dft(real_2, imag_2, buff_size);
#endif
        }
      end = timestamp_get_current_time();
      time_1 = (unsigned)(end - start);
    
      start = timestamp_get_current_time();
      for(test=0; test<NUM_TESTS; test++)
        {
          fill_buffers(real_1, real_2, imag_1, imag_2, buff_size);
#ifdef INCLUDE_FORWARD_TRANSFORMS
          dft_real_forward_dft(real_1, imag_1, buff_size);
          dft_real_forward_dft(real_2, imag_2, buff_size);
#endif
#ifdef INCLUDE_INVERSE_TRANSFORMS
          dft_real_inverse_dft(real_1, imag_1, buff_size);
          dft_real_inverse_dft(real_2, imag_2, buff_size);
#endif
        }
      end = timestamp_get_current_time();
      time_2 = (unsigned)(end - start);
    
      start = timestamp_get_current_time();
      for(test=0; test<NUM_TESTS; test++)
        {
          fill_buffers(real_1, real_2, imag_1, imag_2, buff_size);
#ifdef INCLUDE_FORWARD_TRANSFORMS
          rdft_real_forward_dft(real_1, buff_size);
          rdft_real_forward_dft(real_2, buff_size);
#endif
#ifdef INCLUDE_INVERSE_TRANSFORMS
          rdft_real_inverse_dft(real_1, buff_size);
          rdft_real_inverse_dft(real_2, buff_size);
#endif
        }
      end = timestamp_get_current_time();
      time_3 = (unsigned)(end - start);
    
      start = timestamp_get_current_time();
      for(test=0; test<NUM_TESTS; test++)
        {
          fill_buffers(real_1, real_2, imag_1, imag_2, buff_size);
#ifdef INCLUDE_FORWARD_TRANSFORMS
          dft_2_real_forward_dfts(real_1, real_2, imag_1, imag_2, buff_size);
#endif
#ifdef INCLUDE_INVERSE_TRANSFORMS
          dft_2_real_inverse_dfts(real_1, real_2, imag_1, imag_2, buff_size);
#endif
        }
      end = timestamp_get_current_time();
      time_4 = (unsigned)(end - start);
    
    
      start = timestamp_get_current_time();
      for(test=0; test<NUM_TESTS; test++)
        {
          fill_buffers(real_1, real_2, imag_1, imag_2, buff_size);
#ifdef INCLUDE_FORWARD_TRANSFORMS
          rdft_2_real_forward_dfts(real_1, real_2, buff_size);
#endif
#ifdef INCLUDE_INVERSE_TRANSFORMS
          rdft_2_real_inverse_dfts(real_1, real_2, buff_size);
#endif
        }
      end = timestamp_get_current_time();
      time_5 = (unsigned)(end - start);
      
      /*
      fprintf(stderr, "%i trials of length %i\r\n\n", NUM_TESTS, buff_size);
      fprintf(stderr, "2x dft_complex_forward_dft: %u us (baseline)\r\n\n", time_1);
      fprintf(stderr, "2x dft_real_forward_dft: %u us (%.2fx improvement)\r\n\n", time_2, ((double)time_1/(double)time_2));
      fprintf(stderr, "2x rdft_real_forward_dft: %u us (%.2fx improvement)\r\n\n", time_3, ((double)time_1/(double)time_3));
      fprintf(stderr, "1x dft_2_real_forward_dfts: %u us (%.2fx improvement)\r\n\n", time_4, ((double)time_1/(double)time_4));
      fprintf(stderr, "1x rdft_2_real_forward_dfts: %u us (%.2fx improvement)\r\n\n", time_5, ((double)time_1/(double)time_5));
      */
      
      fprintf(stderr, "%i\t%i\t%i\t%i\t%i\t%i\r\n", buff_size, time_1, time_2, time_3, time_4, time_5);

      free(real_1);
      free(real_2);
      free(imag_1);
      free(imag_2);
    }

  return 0;
}
