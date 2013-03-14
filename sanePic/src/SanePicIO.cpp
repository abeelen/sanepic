#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_math.h>

#include "ImageIO.h"
#include "TemporaryIO.h"
#include "MPIConfiguration.h"
#include "StructDefinition.h"
#include "SanePicIO.h"
#include "InputFileIO.h"



using namespace std;



int writeMapsToFits(string fname, double *S, double *Mptot, long long addnpix,
		long NAXIS1, long NAXIS2, long long *indpix, long long *indpsrc,
		long long npixsrc, int factdupl, long ntotscan,
		struct wcsprm *wcs, char * subheader, int nsubkeys, bool extend){

	double *map1d_d;
	long mi;

	// Build Image
	map1d_d= new double [NAXIS1*NAXIS2];

	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d_d[mi] = S[indpix[mi]];
			} else {
				if (addnpix){
					double b = 0.;
					for (long iframe = 0;iframe<ntotscan;iframe++){
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1) && (indpix[ll] != -1)){
							map1d_d[mi] += S[indpix[ll]]/gsl_pow_2(Mptot[indpix[ll]]);
							b += 1./gsl_pow_2(Mptot[indpix[ll]]);
						}
					}
					if (b >= 1) map1d_d[mi] /= b; else map1d_d[mi] = NAN;
				} else
					map1d_d[mi] = NAN;
			}
		}
	}


	if(write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d_d, (char *)"Image", extend, subheader, nsubkeys)){
		cerr << "EE - Error Writing Image map... \n";
		return EXIT_FAILURE;
	}

	// Build Error
	fill(map1d_d,map1d_d+NAXIS1*NAXIS2,0.0);

	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d_d[mi] = Mptot[indpix[mi]];
			} else {
				if (addnpix){

					for (long iframe = 0;iframe<ntotscan;iframe++){

						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1)  && (indpix[ll] != -1))
							map1d_d[mi] += 1./gsl_pow_2(Mptot[indpix[ll]]) ;
					}
					if (map1d_d[mi] != 0)
						map1d_d[mi] = 1./sqrt(map1d_d[mi]); else map1d_d[mi] = NAN;
				} else
					map1d_d[mi] = NAN;
			}
		}
	}

	if(write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d_d, (char *)"Error",true, subheader, nsubkeys)){
		cerr << "EE - Error Writing Error map... \n";
	}

	delete [] map1d_d;


	return EXIT_SUCCESS;

}

int exportExtraToFits(string fname, struct param_common dir, struct samples samples_struct,
		struct param_saneProc Proc_param, struct param_sanePos Pos_param, struct param_sanePS PS_param,
		struct param_sanePic Pic_param, struct param_saneInv Inv_param){

	if (Pos_param.maskfile != "")
		if(copy_fits_mask(fname, dir.input_dir + Pos_param.maskfile))
			cerr << "WW - No mask will be included in the file : " << fname << endl;

	if(write_fits_inifile(fname,  dir, Proc_param, Pos_param , PS_param, Pic_param, Inv_param)) // write sanePre parameters in naive Map fits file header
		cerr << "WW - No history will be included in the file : " << fname << endl;
	if(write_fits_inputfile(fname, samples_struct))
		cerr << "WW - No input files will be included in the file : " << fname << endl;

	// end of write map function
	return EXIT_SUCCESS;

}

