/*
 * inline_IO.h
 *
 *  Created on: 23 juin 2009
 *      Author: matthieu
 */

#ifndef INLINE_IO_H_
#define INLINE_IO_H_


#include <cstdlib>
#include <string>
#include <cmath>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>

#include <fftw3.h>

extern "C" {
#include "nrutil.h"
}

/*
#ifndef INLINE
# if __GNUC__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif*/
using namespace std;

// sanePos functions
void write_info_pointing(int nn, string outdir, string termin, int coordsyst, double *tanpix, double *tancoord);

void write_samptopix(long ns, long *samptopix, string termin, string outdir, int idet, long iframe);

void write_indpix(long ind_size, int npix, long *indpix, string termin, string outdir, int flagon);

//sanePre functions

void read_info_pointing(int &nn, string outdir, string termin, int &coordsyst2, double *&tanpix, double *&tancoord);

void read_indpix(long ind_size, int &npix, long *&indpix, string termin, string outdir, int &flagon);


void write_PNd(double *PNd, int npix, string termin, string outdir);

void read_PNd(double *&PNdtot, int &npix, string termin, string outdir);

void read_samptopix(long ns, long *&samptopix, string termin, string outdir, int idet, long iframe);

void write_fdata(long ns, fftw_complex *fdata, string termin, string outdir, int idet, long iframe);


void read_noise_file(long &nbins, double *&ell, double *&SpN, double **&SpN_all, string nameSpfile, long ndet);

void read_fdata(long ns, fftw_complex *&fdata, string prefixe, string termin, string outdir, int idet, long iframe);

void write_fPs(long ns, fftw_complex *fdata, string termin, string outdir, long idet, long iframe);

void write_info_for_second_part(string outdir, string termin, int nn, int npix,
		double pixdeg, double *tancoord, double* tanpix, int coordsyst, bool flagon, long* indpix);


#endif
