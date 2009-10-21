
#ifndef IMAGEIO_H_
#define IMAGEIO_H_

#include <iostream>
#include <string>

using namespace std;

extern "C" {
#include <fitsio.h>
#include "wcslib/wcs.h"
}

void print_fits_error(int status);
void write_fits(string fname, double pixsize, long nx, long ny, double *tancoord, double *tanpix, int coordsyst, char dtype, void *data);
void write_fits_wcs(string fname, struct wcsprm * wcs, unsigned long NAXIS1, unsigned long NAXIS2, char dtype, void *data);
void read_fits_signal(string fname, double *S, long long* indpix, long &NAXIS1, long &NAXIS2, long long npix);

#endif /* IMAGEIO_H_ */
