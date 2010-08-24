#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <gsl/gsl_math.h>

#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "write_maps_to_disk.h"

using namespace std;


int write_maps_to_disk(double *S, long NAXIS1, long NAXIS2, string outdir, long long *indpix, long long *indpsrc,
		double *Mptot, long long addnpix, long long npixsrc, int factdupl, long ntotscan,
		struct param_process proc_param, struct param_positions pos_param, struct detectors det,
		struct samples samples_struct, std::vector<double> fcut, struct wcsprm *wcs, string maskfile){



	ostringstream temp_stream;
	double *map1d;
	string fname;
	long mi;

	map1d= new double [NAXIS1*NAXIS2];


	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d[mi] = S[indpix[mi]];
			} else {
				if (addnpix){
					double b = 0.;
					for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1) && (indpix[ll] != -1)){
							map1d[mi] += S[indpix[ll]]/gsl_pow_2(Mptot[indpix[ll]]);
							b += 1./gsl_pow_2(Mptot[indpix[ll]]);
						}
					}
					if (b >= 1) map1d[mi] /= b; else map1d[mi] = NAN;
				} else
					map1d[mi] = NAN;
			}
		}
	}

	//	fname = '!' + outdir + "optimMap_flux.fits";
	fname = outdir + "optimMap_sanePic.fits";
	if(write_fits_wcs("!" + fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Image", 0)){
		cerr << "Error Writing map : EXITING ... \n";
		return 1;
	}


	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d[mi] = Mptot[indpix[mi]];
			} else {
				if (addnpix){

					for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){

						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1)  && (indpix[ll] != -1))
							map1d[mi] += 1./gsl_pow_2(Mptot[indpix[ll]]) ;
					}
					if (map1d[mi] != 0)
						map1d[mi] = 1./sqrt(map1d[mi]); else map1d[mi] = NAN;
				} else
					map1d[mi] = NAN;
			}
		}
	}


	//	fname = '!' + outdir + "optimMap_noisevar.fits";
	if(write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Error",1)){
		cerr << "Error Writing map : EXITING ... \n";
	}
	//
	//	if (addnpix){
	//		// initialize the container
	//		for (long jj=0; jj<NAXIS2 ; jj++){
	//			for (long ii=0; ii<NAXIS1 ; ii++){
	//				mi = jj*NAXIS1 + ii;
	//				map1d[mi] = 0.0;
	//			}
	//		}
	//		// loop thru frame to coadd all pixels
	//		for (long ii=0; ii<NAXIS1; ii++) {
	//			for (long jj=0; jj<NAXIS2; jj++)	 {
	//				mi = jj*NAXIS1 + ii;
	//				double b = 0.;
	//				for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
	//					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
	//					if ((indpsrc[mi] != -1) && (indpix[ll] != -1)){
	//						map1d[mi] += S[indpix[ll]]/gsl_pow_2(Mptot[indpix[ll]]);
	//						b += 1./gsl_pow_2(Mptot[indpix[ll]]);
	//					}
	//				}
	//				if (b >= 1)
	//					map1d[mi] /= b;
	//			}
	//		}
	//		// replace the non observed pixels by NAN
	//		for (long ii=0; ii<NAXIS1; ii++) {
	//			for (long jj=0; jj<NAXIS2; jj++) {
	//				mi = jj*NAXIS1 + ii;
	//				if (map1d[mi] == 0.0)
	//					map1d[mi] = NAN;
	//			}
	//		}
	//
	//		write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, "CCR Image", 1);
	//
	//		for (long ii=0; ii<NAXIS1; ii++) {
	//			for (long jj=0; jj<NAXIS2; jj++) {
	//				mi = jj*NAXIS1 + ii;
	//				for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
	//
	//					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
	//					if ((indpsrc[mi] != -1)  && (indpix[ll] != -1))
	//						map1d[mi] += 1./gsl_pow_2(Mptot[indpix[ll]]) ;
	//				}
	//				if (map1d[mi] != 0)
	//					map1d[mi] = 1./sqrt(map1d[mi]);
	//			}
	//		}
	//		write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"CCR Error",1);
	//	}

	if (maskfile != "")
		if(write_fits_mask(fname, maskfile))
			cerr << "WARNING ! No mask will be included in the file : " << fname << endl;

	if(write_fits_hitory(fname , NAXIS1, NAXIS2, outdir, proc_param, pos_param , fcut, det, samples_struct))
		cerr << "WARNING ! No history will be included in the file : " << fname << endl;

	// clean
	delete [] map1d;

	// end of write map function

	return 0;

}
