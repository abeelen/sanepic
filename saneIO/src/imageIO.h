
#ifndef IMAGEIO_H_
#define IMAGEIO_H_

#include <iostream>
#include <string>
#include <vector>

#include "struct_definition.h"

extern "C" {
#include <fitsio.h>
#include "wcslib/wcs.h"
}

using namespace std;



void print_fits_error(int status);
//void write_fits(string fname, double pixsize, long nx, long ny, double *tancoord, double *tanpix, int coordsyst, char dtype, void *data);
//void write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2, char dtype, void *data);
void write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data, string table_name ,bool fits_already_exist);
void write_fits_hitory(std::string fnaivname,long NAXIS1, long NAXIS2, struct param_process proc_param, struct param_positions pos_param, std::vector<double> fcut, struct detectors det, struct samples samples_struct, long ncomp=-1);
int read_mask_wcs(string fname, string extname, /*char dtype,*/ struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2,  short *& data);

void read_fits_signal(string fname, double *S, long long* indpix, long &NAXIS1, long &NAXIS2, long long npix);

void save_MapHeader(string outdir, struct wcsprm * wcs, long NAXIS1, long NAXIS2);
void read_MapHeader(string outdir, struct wcsprm *& wcs, long *NAXIS1, long *NAXIS2);
void print_MapHeader(struct wcsprm * wcs);


#endif /* IMAGEIO_H_ */

