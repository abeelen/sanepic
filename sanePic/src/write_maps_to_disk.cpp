#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <gsl/gsl_math.h>

#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "write_maps_to_disk.h"
#include "inputFileIO.h"

using namespace std;


int write_maps_to_disk(double *S, long NAXIS1, long NAXIS2, long npix, struct param_common dir, long long *indpix, long long *indpsrc,
		double *Mptot, long long addnpix, long long npixsrc, int factdupl, long ntotscan,
		struct param_sanePre proc_param, struct param_sanePos pos_param,
		struct samples samples_struct, std::vector<double> fcut, struct wcsprm *wcs, string maskfile, struct param_sanePS structPS,
		struct param_sanePic sanePic_struct, struct param_saneInv saneInv_struct, char * subheader, int nsubkeys,
		std::vector<std::vector<std::string> > bolo_vect){



	ostringstream temp_stream;
	double *map1d_d;
	long *map1d_l;
	string fname;
	long mi;

	string outdir = dir.output_dir;
	map1d_d= new double [NAXIS1*NAXIS2];


	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d_d[mi] = S[indpix[mi]];
			} else {
				if (addnpix){
					double b = 0.;
					for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
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

	fname = outdir + sanePic_struct.map_prefix + "_sanePic.fits";
	if(write_fits_wcs("!" + fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d_d, (char *)"Image", 0, subheader, nsubkeys)){
		cerr << "Error Writing map : EXITING ... \n";
		return 1;
	}

	fill(map1d_d,map1d_d+NAXIS1*NAXIS2,0.0);

	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d_d[mi] = Mptot[indpix[mi]];
			} else {
				if (addnpix){

					for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){

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

	if(write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d_d, (char *)"Error",1, subheader, nsubkeys)){
		cerr << "Error Writing map : EXITING ... \n";
	}

	delete [] map1d_d;
	map1d_l= new long [NAXIS1*NAXIS2];

	//compute hits map
	long *hits;
	long long *samptopix;
	hits  = new long[npix];
	fill(hits,hits+npix,0.0);

	for (long iframe=0;iframe<samples_struct.ntotscan;iframe++){
		long ns = samples_struct.nsamples[iframe];
		samptopix=new long long [ns];

		std::vector<string> det_vect = bolo_vect[iframe];
		long ndet = (long)det_vect.size();

		for (long idet1=0;idet1<ndet;idet1++){

			string field1 = det_vect[idet1];

			if(read_samptopix(samples_struct.dirfile_pointer, ns, samptopix, samples_struct.basevect[iframe], field1))
				return 1;
			//compute hit counts
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}
		delete [] samptopix;
	}


	for (long jj=0; jj<NAXIS2; jj++) {
		for (long ii=0; ii<NAXIS1; ii++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d_l[mi] = hits[indpix[mi]];
			} else {
				map1d_l[mi] = 0;
			}
		}
	}

	if (addnpix){
		for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++){
			for (long jj=0; jj<NAXIS2; jj++) {
				for (long ii=0; ii<NAXIS1; ii++) {
					mi = jj*NAXIS1 + ii;
					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
					if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
						map1d_l[mi] += hits[indpix[ll]];
				}
			}
		}
	}

	if(	write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'l', (void *)map1d_l,"Coverage",1, subheader, nsubkeys)){
		cerr << "Error Writing coverage map  ... \n";
	}

	delete [] map1d_l;


	if (maskfile != "")
		if(write_fits_mask(fname, dir.input_dir + maskfile))
			cerr << "WARNING ! No mask will be included in the file : " << fname << endl;

	if(write_fits_history2(fname, NAXIS1, NAXIS2, dir, proc_param, pos_param , samples_struct.fcut, samples_struct, structPS, sanePic_struct, saneInv_struct)) // write sanePre parameters in naive Map fits file header
		cerr << "WARNING ! No history will be included in the file : " << fname << endl;

	// clean
	delete [] hits;


	// end of write map function
	return 0;

}
