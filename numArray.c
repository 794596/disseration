#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include "numArray.h"

/**
 * Function file_to_array.
 * Opens external file and converts the output into an array
 *
 * @param length size of the array that is being created.
 * @param file the file that is being opened and the content extracted.
 */
char *file_to_array(long long length, char *file)
{
  char number[length];
  char *str;
  str = (char *) malloc(sizeof(number));
  //opens file
  FILE *fp = fopen(file, "rb");
  // prints  if the file fails to opens
  if (fp == NULL) printf("Failed to open the file.\n");
  //gets string from file and adds it to number
  fgets(number, length, fp);
  //closes file
  fclose(fp);

  //copies string from number to str
  strcpy(str, number);

  return str;
  // frees memory for str
  free(str);
}
