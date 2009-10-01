#include <stdio.h>

#include "binaryFileIO.h"


// TODO: remove those functions : unused anywhere !!!!

void write_vector(char *filename, void *data, int typesize, long nn) {

  FILE *fp;

  fp = fopen(filename,"w");
  fwrite(data,typesize, nn, fp);
  fclose(fp);

}


void read_vector(char *filename, void *data, int typesize, long nn) {

  FILE *fp;
  size_t result;

  fp = fopen(filename,"r");
  result = fread(data,typesize, nn, fp);
  fclose(fp);

}
