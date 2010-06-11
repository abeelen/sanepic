


#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "write_maps_to_disk.h"

using namespace std;


void write_maps_to_disk(double *S, long NAXIS1, long NAXIS2, string outdir, long long *indpix, long long *indpsrc,
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
				map1d[mi] = NAN;
			}
		}
	}

	//	fname = '!' + outdir + "optimMap_flux.fits";
	fname = outdir + "optimMap_sanePic.fits";
	write_fits_wcs("!" + fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Image", 0);


	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d[mi] = Mptot[indpix[mi]];
			} else {
				map1d[mi] = NAN;
			}
		}
	}


	//	fname = '!' + outdir + "optimMap_noisevar.fits";
	write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Error",1);

	if (addnpix){
		for (long iframe = 0;iframe<ntotscan;iframe++){
			for (long ii=0; ii<NAXIS1; ii++) {
				for (long jj=0; jj<NAXIS2; jj++) {
					mi = jj*NAXIS1 + ii;
					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
					if ((indpsrc[mi] != -1)  && (indpix[ll] != -1)){
						map1d[mi] = S[ll];
					} else {
						map1d[mi] = 0.0;
					}
				}
			}



			//			temp_stream << "!" + outdir + "optimMap_flux_fr" << iframe << ".fits";
			//			// Transform into string
			//			fname= temp_stream.str();
			//			// Clear ostringstream buffer
			//			temp_stream.str("");
			write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *) "Duplicated Image",1);

			for (long ii=0; ii<NAXIS1; ii++) {
				for (long jj=0; jj<NAXIS2; jj++) {
					mi = jj*NAXIS1 + ii;
					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
					if ((indpsrc[mi] != -1)  && (indpix[ll] != -1)){
						map1d[mi] = Mptot[indpix[ll]];
					} else {
						map1d[mi] = 0.0;
					}
				}
			}

			//fname = '!' + outdir + "optimMap_" + termin + "_noisevar_fr" + iframestr + ".fits";
			//			temp_stream << "!" + outdir + "optimMap_noisevar_fr" << iframe << ".fits";

			// Transform into string
			//			fname= temp_stream.str();
			//			// Clear ostringstream buffer
			//			temp_stream.str("");
			//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
			write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Duplicated Error",1);
		}
	}

	write_fits_hitory(fname , NAXIS1, NAXIS2, outdir, proc_param, pos_param , fcut, det, samples_struct);

	if (maskfile != "")
		write_fits_mask(fname, maskfile);

	// clean
	delete [] map1d;

	// end of write map function


}
