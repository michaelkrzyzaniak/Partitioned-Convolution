#include "DFT.h"
#include "Timestamp.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define BUFFER_SIZE (1<<5)

void print_binary(int x, int n)
{
  while(n > 1)
    {
       n >>= 1;
      fprintf(stderr, "%i", (x & n) ? 1 : 0);
    }
}

int main(void)
{
  dft_sample_t buffer[BUFFER_SIZE];
  
  int i;
  
  for(i=0; i<BUFFER_SIZE; i++)
    buffer[i] = i;
    
  rdft_bit_reverse_indices(buffer, BUFFER_SIZE);
  
  //for(i=0; i<BUFFER_SIZE; i++)
    //fprintf(stderr, "[%i] -> [%i]\r\n", i, (int)buffer[i]);
  

  //for(i=2; i<BUFFER_SIZE; i+=2)
    //fprintf(stderr, "([%i], [%i]) -> ([%i], [%i])\r\n", i, (int)buffer[BUFFER_SIZE-(int)buffer[i]], (int)buffer[i], BUFFER_SIZE-(int)buffer[i]);
    
  for(i=2; i<BUFFER_SIZE; i+=2)
    {
      int x1 = i;
      int x2 = (int)buffer[BUFFER_SIZE-(int)buffer[i]];
      fprintf(stderr, "[%i], [%i] ==> [", x1, x2);
      print_binary(x1, BUFFER_SIZE);
      fprintf(stderr, "], [");
      print_binary(x2, BUFFER_SIZE);
      fprintf(stderr, "]\r\n");
    }
}

/*
[2], [3] ==> [00010], [00011]	1 bit
[4], [7] ==> [00100], [00111]	2 bits
[6], [5] ==> [00110], [00101]	2 bits
[8], [15] ==> [01000], [01111]	3
[10], [13] ==> [01010], [01101]	3
[12], [11] ==> [01100], [01011]	3
[14], [9] ==> [01110], [01001]	3
[16], [31] ==> [10000], [11111]	4
[18], [29] ==> [10010], [11101]	4
[20], [27] ==> [10100], [11011]	4
[22], [25] ==> [10110], [11001]	4
[24], [23] ==> [11000], [10111]	4
[26], [21] ==> [11010], [10101]	4
[28], [19] ==> [11100], [10011]	4
[30], [17] ==> [11110], [10001]	4
*/
