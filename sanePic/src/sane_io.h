#ifndef SANE_IO
#define SANE_IO

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <list>

using namespace std;

extern "C" {
#include <fitsio.h>
}

void print_fits_error(int status);
void write_fits(string fname, double pixsize, long nx, long ny, double *tancoord, double *tanpix, int coordsyst, char dtype, void *data);

void write_psd_tofits(string fname, long nx, long ny, char dtype, void * data);

void write_vector(char *filename, void *data, int typesize, long nn);
void read_vector(char *filename, void *data, int typesize, long nn);
void read_bolofile(string fname, list<string>& bolos);
void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets);

#endif
