#ifndef BOLOIO_H_
#define BOLOIO_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

void read_bolofile(string fname, std::vector<string>& bolos);
void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets);

#endif /* BOLOIO_H_ */
