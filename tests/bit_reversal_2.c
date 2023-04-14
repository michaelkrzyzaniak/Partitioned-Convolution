#include <stdio.h>

int main()
{
  int N = 16;
  int N_over_2 = N >> 1;
  int N_Minus_1 = N - 1;
  unsigned n, bit, rev;
  unsigned n_reversed = N_over_2;
  
  for(n=1; n<N_Minus_1; n++)
    {
      fprintf(stderr, "n: %i, reverse: %i\r\n", n, n_reversed);
      
      bit = ~n & (n + 1);
      rev = N_over_2 / bit;
      n_reversed ^= (N - 1) & ~(rev - 1);
    }
}
