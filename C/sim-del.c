#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "mt.h"

#define ITER 10000000

int main(int argc, char **argv) {

   char **ignore;

   unsigned int K    = atoi(argv[1]);
   unsigned int D    = atoi(argv[2]);
   double       prob = strtod(argv[3], ignore);
   double       del  = strtod(argv[4], ignore);

   if (K == 0 || D == 0 || prob == 0 || del == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   // Set the random seed.
   seedMT(123);

   const unsigned long int p = (prob * 4294967295);
   const unsigned long int d = (del * 4294967295);

   long int total = ITER;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      // Initialize stack.
      int stack = 0;
      for (int i = 0 ; i < K ; i++) {
         if (stack > 0) {
            if (randomMT() < p) {
               // Substitution. Reset the stack.
               stack = 0;
            }
            else if (randomMT() < d) {
               // Deletion. Set the stack to 1.
               stack = 1;
            }
            else {
               // No error. Increase the stack
               stack++;
               if (stack >= D) {
                  // We have a seed.
                  total--;
                  break;
               }
            }
         }
         else {
            // The stack is 0.
            stack = randomMT() > p;
         }
      }

   }

   fprintf(stdout, "%d\t%ld\n", K, total);

}
