#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "mt.h"

#define ITER 10000000

int main(int argc, char **argv) {

   char **ignore;

   unsigned int K     = atoi(argv[1]);
   unsigned int D     = atoi(argv[2]);
   double       prob  = strtod(argv[3], ignore);
   double       kappa = strtod(argv[4], ignore);

   if (K == 0 || D == 0 || prob == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   const unsigned long int p = (prob * 4294967295);
   const unsigned long int kp = (kappa * 4294967295);
   const unsigned long int ak = ((1-kappa/3) * 4294967295);

   // Set the random seed.
   seedMT(123);

   long int total = 0; //ITER;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      int false_found = 0;
      int true_found = 0;

      // Initialize stacks.
      int stack_true = 0;
      int stack_false = 0;

      for (int i = 0 ; i < K ; i++) {
         if (randomMT() < p) {
            // Error. Reset the stack.
            stack_true = 0;
            if (randomMT() < ak) {
               // Error and the duplicat has a different nucleotide.
               // Reset the stack of the false positive.
               stack_false = 0;
            }
            else {
               // Error and duplicate has the same nucleotide.
               // Increase the stack of the false positive.
               stack_false++;
               if (stack_false >= D) false_found = 1;
            }
         }
         else {
            // No error. Increase the stack
            stack_true++;
            if (randomMT() < kp) {
               // No error and duplicate has a different nucleotide.
               // Reset the stack of the false positive.
               stack_false = 0;
            }
            else {
               // No error and duplicate has the same nucleotide.
               // Increase the stack of the false positive.
               stack_false++;
               if (stack_false >= D) false_found = 1;
            }
            if (stack_true >= D) {
               // We have a seed.
//               total--;
               true_found = 1;
               break;
            }
         }
      }

      if (false_found && !true_found) total++;

   }

   fprintf(stdout, "%d\t%ld\n", K, total);

}
