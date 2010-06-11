/*
 * sanePre_preprocess.cpp
 *
 *  Created on: 25 juin 2009
 *      Author: matthieu
 */

#include "sanePos_preprocess.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "sanePos_map_making.h"
#include "dataIO.h"
#include "inline_IO2.h"




using namespace std;

void computePixelIndex(string outdir, std::vector<string> bolonames,
		struct samples samples_struct, struct param_process proc_param, struct param_positions pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout)
{

	/*!
	 * \fn Get coordinates of pixels that are seen
	 * Compute the position to pixel projetcion matrices :
	 * One binary file per bolometer and per scan
	 */

	// TODO : samptopix unsigned long
	long long  *samptopix;
	long ndet = bolonames.size();

	string field;
	string fits_file;

	long ns;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		fits_file=samples_struct.fits_table[iframe];
		ns = samples_struct.nsamples[iframe];

		double *ra, *dec, *phi, **offsets;
		//		short *flag;


		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		read_all_bolo_offsets_from_fits(fits_file, bolonames, offsets);
		//		cout << offsets[100][0] << " " << offsets[100][1] << endl;

		// read reference position
		long test_ns;
		read_ReferencePosition_from_fits(fits_file, ra, dec, phi, test_ns);
		//		cout << ra[100] << " " << dec[100] << " " << phi[100] <<  endl;

		if (test_ns != ns) {
			cerr << "Read position does not correspond to frame size : Check !" << endl;
			exit(-1);
		}

		// find the pointing solution at each time stamp for each detector
		struct celprm celestial;
		celini(&celestial);

		// TODO: use the PRJCODE read from the file...
		tanset(&celestial.prj);

		// Warning !!
		// We need to deproject the bolometer into the sky (detector parallelization) and then project into map pixel
		// as we need to process an entire frame to save data per frame, we cannot parallelize by detector,
		// since we would then need to save a 2*ndet*ns matrix

		// First compute the sin and cos of the phi (will be reused ns times..)
		double *cosphi, *sinphi;
		cosphi = new double[ns];
		sinphi = new double[ns];

		for (long ii=0; ii<ns; ii++){
			cosphi[ii] = cos(phi[ii]/180.0*M_PI);
			sinphi[ii] = sin(phi[ii]/180.0*M_PI);
		}

		delete [] phi;

		for (long idet = 0; idet<ndet; idet++){

			string field = bolonames[idet];

			double offxx, offyy, lon, lat, ra_deg, dec_deg;
			int status;

			long long *xx, *yy;
			double *world, *imgcrd, *pixcrd;
			double *theta, *phi;
			int    *wcsstatus;
			theta  = new double[ns];
			phi    = new double[ns];
			world  = new double[2*ns];
			imgcrd = new double[2*ns];
			pixcrd = new double[2*ns];
			xx     = new long long[ns];
			yy     = new long long[ns];
			wcsstatus = new int[ns];

			// First deproject the bolometer position ....
			for (long ii=0; ii <ns; ii++){

				celestial.ref[0] =  ra[ii]*15.0;
				celestial.ref[1] =  dec[ii];
				celset(&celestial);

				//TODO : check this -1 factor... just a stupid convention...
				offxx = (cosphi[ii] * offsets[idet][0]
				                                    - sinphi[ii] * offsets[idet][1])*-1;
				offyy =  sinphi[ii] * offsets[idet][0]
				                                    + cosphi[ii] * offsets[idet][1];

				if (celx2s(&celestial, 1, 0, 0, 0, &offxx, &offyy, &lon, &lat, &ra_deg, &dec_deg, &status) == 1) {
					printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
					continue;
				}

				world[2*ii]   = ra_deg;
				world[2*ii+1] = dec_deg;
			}

			// ... and reproject it back onto the map
			if ((status = wcss2p(wcs, ns, 2, world, phi, theta, imgcrd, pixcrd, wcsstatus))) {
				printf("   wcss2p(1) ERROR %2d \n", status);
				continue;
			}
			// Transform pixel coordinates to pixel index
			for (long ii = 0; ii < ns; ii++){
				xx[ii]  = (long long) (pixcrd[2*ii]-0.5);
				yy[ii]  = (long long) (pixcrd[2*ii+1]-0.5);
			}

			//						for (unsigned long ii = 0; ii < 20; ii++){
			//							cout << world[2*ii] << " " << world[2*ii+1] << " : ";
			//							cout << int(pixcrd[2*ii]-0.5)<< " " << int(pixcrd[2*ii+1]-0.5) << endl;
			//						}

			delete [] world;
			delete [] imgcrd;
			delete [] pixcrd;
			delete [] theta;
			delete [] phi;
			delete [] wcsstatus;



			// Combine position and bolo flags
			// and check

			int *bolo_flag=NULL;

			long test_ns;
			read_flag_from_fits(fits_file, field, bolo_flag, test_ns);
			if (test_ns != ns) {
				cerr << "Read flags does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}



			for (long ii=0; ii<ns; ii++){

				// TODO : Update this to read the flag corresponding to the channel...
				//				if (flpoint[ii] == 1) bolo_flag[ii] = 1;
				if (bolo_flag[ii] != 0) bolo_flag[ii] = 1;

				if ((xx[ii] < 0)   || (yy[ii] < 0  ))          bolo_flag[ii] = 2;
				if ((xx[ii] >=  NAXIS1) || (yy[ii] >= NAXIS2)) bolo_flag[ii] = 2;
				if (proc_param.NOFILLGAP && (bolo_flag[ii] == 1 ))        bolo_flag[ii] = 2;

				if ((ii < proc_param.napod) || (ii >= ns-proc_param.napod)) bolo_flag[ii] = 3;

			}

			//			cout << "NOFILLGAP : " << NOFILLGAP << endl;
			//			for (unsigned long ii=0; ii<ns; ii++)
			//			cout << ii << " " << (NOFILLGAP && (bolo_flag[ii] == 1) << endl;
			//				if (bolo_flag[ii]==1)
			//					cout << ii << " bolo flag" << endl;

			// Compute sample pixel index according to its flag
			samptopix = new long long[ns];

			for (long ii=0 ; ii<ns; ii++){

				// image1        + crossing constrain removal + flagged pixel
				// NAXIS1*NAXIS2 + addnpix*nframe             + 1

				long long ll=0;
				switch (bolo_flag[ii]) {
				case 0:			// sample is not rejected

					ll = NAXIS1*yy[ii]+xx[ii];
					if (mask[ll] != 0)								//if pixel is in the box constraint removal mask
						ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[ll];
					break;

				case 1:												// sample is flagged

					if (pos_param.flgdupl)								// if flagged pixels are in a duplicated map
						ll = NAXIS1*yy[ii]+xx[ii] + NAXIS1*NAXIS2;	// index in the second map...
					else											// else every flagged sample is projected to the same pixel (outside the map)
						ll = factdupl*NAXIS1*NAXIS2 + addnpix + 2;

					break;

				case 2:												// sample is rejected
					ll = factdupl*NAXIS1*NAXIS2 + addnpix + 1;
					pixout = 1; 									// pixel is out of map
					printf("[%2.2i] %s PIXEL OUT, ii = %ld, xx = %lld, yy = %lld\n",rank, field.c_str(), ii,xx[ii],yy[ii]);
					getchar();
					break;
				case 3:												// apodized data -> flag
					ll = factdupl*NAXIS1*NAXIS2 + addnpix + 2;
					flagon = 1;
					break;
				default:
					break;
				}

				pixon[ll] += 1;
				samptopix[ii] = ll;

			}

			write_samptopix(ns, samptopix,  outdir, iframe, bolonames[idet]);

			delete [] bolo_flag;
			delete [] samptopix;
			delete [] xx;
			delete [] yy;

		}//end of idet loop
		delete [] ra;
		delete [] dec;
		//		delete [] flag;
		free_dmatrix(offsets,(long)0,ndet-1,(long)0,2-1);
		delete [] cosphi;
		delete [] sinphi;

	} // end of iframe loop
}

void computePixelIndex_HIPE(string outdir, std::vector<string> bolonames,
		struct samples samples_struct, struct param_process proc_param, struct param_positions pos_param,long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout)
{

	/*!
	 * \fn Get coordinates of pixels that are seen
	 * Compute the position to pixel projetcion matrices :
	 * One binary file per bolometer and per scan
	 */

	// TODO : samptopix unsigned long
	long long  *samptopix;
	long ndet = bolonames.size();

	string field;

	string fits_file;
	long ns, test_ns;
	int status;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		fits_file=samples_struct.fits_table[iframe];
		ns = samples_struct.nsamples[iframe];

		double *ra, *dec;
		int *flag=NULL;

		for (long idet = 0; idet<ndet; idet++){
			string field = bolonames[idet];

			double *world, *imgcrd, *pixcrd;
			double *theta, *phi;
			long long    *xx, *yy;
			int    *wcsstatus;

			theta  = new double[ns];
			phi    = new double[ns];
			world  = new double[2*ns];
			imgcrd = new double[2*ns];
			pixcrd = new double[2*ns];
			xx     = new long long[ns];
			yy     = new long long[ns];
			wcsstatus = new int[ns];

			read_ra_dec_from_fits(fits_file, field, ra, dec, test_ns);
			if (test_ns != ns) {
				cerr << "Read ra does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}
			//			read_dec_from_fits(fits_file, field, dec, test_ns);
			//			if (test_ns != ns) {
			//				cerr << "Read dec does not correspond to frame size : Check !!" << endl;
			//				exit(-1);
			//			}
			read_flag_from_fits(fits_file, field, flag, test_ns);
			if (test_ns != ns) {
				cerr << "Read flag does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}
			for (long ii=0; ii <ns; ii++){
				world[2*ii]   = ra[ii];
				world[2*ii+1] = dec[ii];
			}
			//			for (long ii=0; ii<ns; ii++){
			//				cout << ii << " " << world[2*ii] << " " << world[2*ii+1] << endl;
			//			}
			// ... and reproject it back onto the map
			//			cout << "before reproject" << endl;

			if ((status = wcss2p(wcs, ns, 2, world, phi, theta, imgcrd, pixcrd, wcsstatus))) {
				printf("   wcss2p(1) ERROR %2d \n", status);
				for (long ii=0; ii<ns; ii++){
					printf("   wcsstatus (%li): %2d\n",ii,wcsstatus[ii]);
					// Something went wrong in the transformation, flag the data point
					if ( wcsstatus[ii] == 1 )
						flag[ii] = -1;
				}
				continue;
			}
			//			cout << "after reproject" << endl;

			// Transform pixel coordinates to pixel index
			for (long ii = 0; ii < ns; ii++){
				xx[ii]  = (long long) (pixcrd[2*ii]-0.5);
				yy[ii]  = (long long) (pixcrd[2*ii+1]-0.5);
			}

			//			for (unsigned long ii = 0; ii < 20; ii++){
			//				cout << world[2*ii] << " " << world[2*ii+1] << " : ";
			//				cout << int(pixcrd[2*ii]-0.5)<< " " << int(pixcrd[2*ii+1]-0.5) << endl;
			//			}

			delete [] world;
			delete [] imgcrd;
			delete [] pixcrd;
			delete [] theta;
			delete [] phi;
			delete [] wcsstatus;



			// Combine position and bolo flags
			// and check

			int *bolo_flag=NULL;

			read_flag_from_fits(fits_file, field, bolo_flag, test_ns);

			if (test_ns != ns) {
				cerr << "Read flags does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}


			for (long ii=0; ii<ns; ii++){

				//TODO: Bad Fix, check flag values and handle properly
				if (flag[ii] != 0)                               bolo_flag[ii] = 1;
				if (flag[ii] == -1)                              bolo_flag[ii] = 2;

				if ((xx[ii] < 0)   || (yy[ii] < 0  ))            bolo_flag[ii] = 2;
				if ((xx[ii] >=  NAXIS1) || (yy[ii] >=  NAXIS2))  bolo_flag[ii] = 2;
				if (proc_param.NOFILLGAP && (bolo_flag[ii] == 1 ))      bolo_flag[ii] = 2;

				if ((ii < proc_param.napod) || (ii >= ns-proc_param.napod))    bolo_flag[ii] = 3;

			}

			//			cout << "NOFILLGAP : " << NOFILLGAP << endl;
			//			for (unsigned long ii=0; ii<ns; ii++)
			//			cout << ii << " " << (NOFILLGAP && (bolo_flag[ii] == 1) << endl;
			//				if (bolo_flag[ii]==1)
			//					cout << ii << " bolo flag" << endl;

			// Compute sample pixel index according to its flag

			samptopix = new long long[ns];

			for (long ii=0 ; ii<ns; ii++){

				// image1        + crossing constrain removal + flagged pixel + apodized data
				// NAXIS1*NAXIS2 + addnpix*nframe             + 1             + 1

				long long ll=0;
				switch (bolo_flag[ii]) {
				case 0:			// sample is not rejected

					ll = NAXIS1*yy[ii]+xx[ii];
					if (mask[ll] != 0)								//if pixel is in the box constraint removal mask
						ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[ll];
					break;

				case 1:												// sample is flagged

					if (pos_param.flgdupl)							// if flagged pixels are in a duplicated map
						ll = NAXIS1*yy[ii]+xx[ii] + NAXIS1*NAXIS2;	// index in the second map...
					else											// else every flagged sample is projected to the same pixel (outside the map)
						ll = factdupl*NAXIS1*NAXIS2 + addnpix + 1;

					break;

				case 2:												// sample is rejected
					ll = factdupl*NAXIS1*NAXIS2 + addnpix + 1;
					pixout = 1; 									// pixel is out of map
					printf("[%2.2i] %s PIXEL OUT, ii = %ld, xx = %lld, yy = %lld\n",rank, field.c_str(), ii,xx[ii],yy[ii]);
					break;
				case 3:												// apodized data -> flag
					ll = factdupl*NAXIS1*NAXIS2 + addnpix + 2;
					flagon = 1;
					break;
				default:
					break;
				}

				pixon[ll] += 1;
				samptopix[ii] = ll;

			}

			write_samptopix(ns, samptopix,  outdir, iframe, bolonames[idet]);

			delete [] bolo_flag;
			delete [] samptopix;
			delete [] xx;
			delete [] yy;

		}
	}

}
