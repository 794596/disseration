#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

int writeResults(long length, double results, char *file) {
  FILE *fp = fopen(file, "a");
   fprintf(fp, "%ld, %f\n", length, results);

  fclose(fp);
  return 0;
}
