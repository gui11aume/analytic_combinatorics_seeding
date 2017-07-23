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

   if (K == 0 || D == 0 || prob == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   const unsigned long int p = (prob * 4294967295);

   // Set the random seed.
   seedMT(123);

   double toterr = 0;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      int redo = 1;
      int nerr;

      // Accept-reject loop.
      while (redo) {

         // Redo if there is a seed.
         redo = 0;

         // Initialize number of errors.
         nerr = 0;
         // Initialize stack.
         int stack = 0;

         for (int i = 0 ; i < K ; i++) {
            if (randomMT() < p) {
               // Error. Reset the stack.
               stack = 0;
               nerr++;
            }
            else {
               // No error. Increase the stack
               stack++;
               if (stack >= D) {
                  // We have a seed. Redo.
                  redo = 1;
                  break;
               }
            }
         }

      }

      toterr += nerr;

   }

   fprintf(stdout, "%d\t%.6f\n", K, toterr / ITER);

}
