
#ifndef BINARYFILEIO_H_
#define BINARYFILEIO_H_

//#include <stdio.h>

void write_vector(char *filename, void *data, int typesize, long nn);
void read_vector(char *filename, void *data, int typesize, long nn);


#endif /* BINARYFILEIO_H_ */
