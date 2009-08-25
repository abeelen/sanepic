

#ifndef DATAIO_H_
#define DATAIO_H_

#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int read_data_std(string fname, int frame, int fs, int ns, void* data, string field, char type); // on garde
// int read_data(string fname, int frame, int fs, int ns, void* data, string field, char type); // on garde

#endif /* DATAIO_H_ */
