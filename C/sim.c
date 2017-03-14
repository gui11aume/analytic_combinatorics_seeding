#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#define ITER 10000000

int main(int argc, char **argv) {

   char **ignore;

   unsigned int K = atoi(argv[1]);
   unsigned int D = atoi(argv[2]);
   double       p = strtod(argv[3], ignore);

   if (K == 0 || D == 0 || p == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   const double pmax = p* RAND_MAX;

   // Set the random seed.
   srand(123);

   long int total = ITER;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      // Initialize stack.
      int stack = 0;
      for (int i = 0 ; i < K ; i++) {
         if (drand48() < p) {
            // Error. Reset the stack.
            stack = 0;
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

   }

   fprintf(stdout, "%ld\n", total);

}
