
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
//void write_fits(string fname, double pixsize, long nx, long ny, double *tancoord, double *tanpix, int coordsyst, char dtype, void *data);
//void write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2, char dtype, void *data);
void write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data, string table_name ,bool fits_already_exist);
int read_mask_wcs(string fname, string extname, /*char dtype,*/ struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2,  short *& data);

//TODO : Rewrite this one...
void read_fits_signal(string fname, double *S, long long* indpix, long &NAXIS1, long &NAXIS2, long long npix);

void save_MapHeader(string outdir, struct wcsprm * wcs, long NAXIS1, long NAXIS2);
void read_MapHeader(string outdir, struct wcsprm *& wcs, long *NAXIS1, long *NAXIS2);
void print_MapHeader(struct wcsprm * wcs);


#endif /* IMAGEIO_H_ */
