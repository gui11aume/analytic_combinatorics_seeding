#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "mt.h"

#define ITER 100000

int main(int argc, char **argv) {

   char **ignore;

   const unsigned int K    = atoi(argv[1]);
   const unsigned int D    = atoi(argv[2]);
   const double       prob = strtod(argv[3], ignore);
   const double       del  = strtod(argv[4], ignore);
   const double       rrob = strtod(argv[5], ignore);
   const double       arob = strtod(argv[6], ignore);

   if (K == 0 || D == 0 || prob == 0 || del == 0 ||
         rrob == 0 || arob == 0) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }
   if (prob + rrob > 1) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   // Set the random seed.
   seedMT(123);

   const unsigned long int p = (prob * 4294967295);
   const unsigned long int d = (del  * 4294967295);
   const unsigned long int r = (rrob * 4294967295);
   const unsigned long int a = (arob * 4294967295);

   double toterr = 0;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {

      int redo = 1;
      int nerr;

      // Accept-reject loop.
      while (redo) {

         redo = 0;       // Redo if there is a seed.
         nerr = 0;       // Initialize number of errors.
         int stack = 0;  // Initialize stack.
         int insertion_mode = 0;

         for (int i = 0 ; i < K ; i++) {
            if (insertion_mode) {
               if (randomMT() > a) {
                  insertion_mode = 0;
                  // Going out of insertion mode. Get a mismatch with
                  // probability p / (1-r).
                  stack = randomMT() > p / (1-rrob);
               }
               else {
                  nerr++;
               }
            }
            else {
               // Not in insertion mode.
               double otherprob = randomMT();
               if (otherprob < p) {
                  // Substitution. Reset the stack.
                  stack = 0;
               }
               else if (otherprob < p+r) {
                  // Insertion. Reset the stack and enter insertion mode.
                  insertion_mode = 1;
                  stack = 0;
                  nerr++;
               }
               // Neither substitution nor insertion.
               else {
                  // If the stack was empty, just set it to 1,
                  // otherwise check for deletions.
                  if (stack == 0) {
                     stack = 1;
                  }
                  else {
                     stack = randomMT() < d ? 1 : stack + 1;
                  }
                  if (stack >= D) {
                     // We have a seed.
                     redo = 1;
                     break;
                  }
               }
            }
         }
      }

      // Generated the read without redo.
      // Add up deletions and move to the next read.
      toterr += nerr;

   }

   fprintf(stdout, "%d\t%.6f\n", K, toterr / ITER);

}
