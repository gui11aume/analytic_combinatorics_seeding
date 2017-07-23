#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "mt.h"

int main(int argc, char **argv) {

   char **ignore;

   unsigned int K     = atoi(argv[1]);
   unsigned int D     = atoi(argv[2]);
   double       prob  = strtod(argv[3], ignore);
   double       kappa = strtod(argv[4], ignore);
   long int     ITER  = atoi(argv[5]);

   if (K == 0 || D == 0 || prob == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   const unsigned long int p = (prob * 4294967295);
   const unsigned long int kp = (kappa * 4294967295);
   const unsigned long int ak = ((1-kappa/3) * 4294967295);

   // Set the random seed.
   seedMT(123);

   long int total = ITER;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      // Initialize stacks.
      int stack_true = 0;
      int stack_false = 0;
      int false_positive = 0;

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
               if (stack_false >= D) {
                  false_positive = 1;
               }
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
               if (stack_false >= D) {
                  false_positive = 1;
               }
            }
            if (stack_true >= D) {
               // We have a true positive.
               false_positive = 0;
               break;
            }
         }
      }

      if (!false_positive) total--;

   }

   fprintf(stdout, "%d\t%ld\n", K, total);

}
