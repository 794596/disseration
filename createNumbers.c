#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include "local.h"

/**
 * Function randoms.
 * generates a random number as an int and returns it as a string
 *
 */
char randoms()
{
  int place = rand()%10;
  char got[2];
  got[0] = place+'0';
  return *got;
}

/**
 * Function create_number.
 * Opens file an inputs a string of integers and then closes.
 *
 * @param length number of integers inserted onto the file.
 * @param file the file that is being opened and the content added.
 */
void create_number(long length, char *file)
{
  FILE *fptr;
  fptr = fopen(file , "w");
  if (fptr != NULL) {
    printf("File created successfully!\n");
  }
  else {
    printf("Failed to create the file.\n");
  }
  for (int i = 0; i < length; ++i){
    fputc(randoms(), fptr);
  }
  fclose(fptr);
}

void main(int argc, char *argv[])
{
  create_number(1000000, "number1.txt");
  create_number(1000000, "number2.txt");
}
