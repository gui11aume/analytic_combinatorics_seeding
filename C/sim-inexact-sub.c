#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "mt.h"

int main(int argc, char **argv) {

   char **ignore;

   unsigned int      KK    = atoi(argv[1]);
   unsigned int      D     = atoi(argv[2]);
   double            prob  = strtod(argv[3], ignore);
   long unsigned int ITER  = (long unsigned int) atoi(argv[4]);

   if (KK == 0 || D == 0 || prob == 0 || ITER == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   // Set the random seed.
   seedMT(123);

   const unsigned int long p = (prob * 4294967295);

   long int total = ITER;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      // Initialize stacks. It is important to set the old stack to
      // report a seed in case the first mutation is at position D-1.
      int old_stack = 1;
      int new_stack = 0;
      for (int i = 0 ; i < KK ; i++) {
         if (randomMT() < p) {
            // Substitution. Shift stacks and reset the new one.
            old_stack = new_stack + 1;
            new_stack = 0;
         }
         else {
            // No substitution. Increase the new stack
            new_stack++;
            if (old_stack + new_stack >= D) {
               // We have a seed.
               total--;
               break;
            }
         }
      }

   }

   fprintf(stdout, "%d\t%ld\n", KK, total);

}
