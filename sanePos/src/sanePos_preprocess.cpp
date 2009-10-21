/*
 * sanePre_preprocess.cpp
 *
 *  Created on: 25 juin 2009
 *      Author: matthieu
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <cmath>

#include <fcntl.h>
#include <unistd.h>

#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>


#include "sanePos_map_making.h"
#include "positionsIO.h"
#include "dataIO.h"
#include "blastSpecific.h"
#include "inline_IO2.h"
#include "positionsIO.h"
#include "mpi_architecture_builder.h"

#include "sanePos_preprocess.h"

//// TODO: Will not be used soon, replaced by compute MapMinima
//
//void find_coordinates_in_map(long ndet,std::vector<string> bolonames, string *fits_table,/*,string bextension,
//		string fextension,*//*string file_offsets,foffset *foffsets,float *scoffsets,	double *offsets,*/long iframe_min, long iframe_max,
//		/*long *fframes,*/long *nsamples,string dirfile,/*string ra_field,string dec_field,string phi_field, string scerr_field,
//		string flpoint_field,int nfoff,*/ double pixdeg, int *&xx, int *&yy,int nn, double *&coordscorner, double *tancoord,
//		double *tanpix, bool bfixc, double radius, /*double *offmap,*/ double *srccoord,char type,double *&ra,double *&dec,
//		double *&phi, short *&flpoint,double &ra_min,double &ra_max,double &dec_min,double &dec_max,bool default_projection){
//
//
//
//	double *offsets, *froffsets;
//	offsets = new double[2]; //
//	froffsets = new double[2]; //
//
//	//short *flpoint2;
//	//double *ra2, *dec2, *phi2;
//	//double *signal, *signal2;
//
//	//signal = new double[29341];
//	//signal2 = new double[29341];
//
//	//flpoint2 = new short[29341];
//	//ra2 = new double[29341];
//	//dec2 = new double[29341];
//	//phi2 = new double[29341];
//
//	string field;
//	//string bolofield;
//	//string flagfield;
//	string fits_file;
//	long ns;
//	//double *ra,*dec,*phi,*scerr, *flpoint;
//
//	// ndet = number of channels
//	for (int idet=0;idet<ndet;idet++){
//		//cout << idet << endl;
//
//		// field = actual boloname
//		field = bolonames[idet];
//
//		// bolofield = boloname + bextension
//		//bolofield = field+bextension;
//
//		//if (cextension != "NOCALP")
//		//	calfield  = field+cextension;
//		//if (fextension != "NOFLAG")
//		//flagfield = field+fextension;
//
//		//read bolometer offsets
//		//read_bolo_offsets(field,file_offsets,scoffsets,offsets);
//		//cout << "before " << offsets[0] << " " << offsets[1] << endl;
//
//		fits_file=fits_table[0];
//
////		cout << fits_file << " " << field << endl;
//
//		read_bolo_offsets_from_fits(fits_file, field, offsets);
//
//
////		cout << "after  " << offsets[0] << " " << offsets[1] << endl;
////		exit(0);
//
//		///////////////// Procedure de test pour read_data_from_fits /////////////
//
//		//read_signal_from_fits(fits_file, signal, field);
//
//		//cout << "signal : " << signal[0] << " " << signal[1] << " " << signal[2] << " " << signal[3] << endl;
//
//		/*ff=0;
//		ns=29340;
//		read_data_std(dirfile, ff, 0, ns, signal2, field+"_data", 'd');*/
//
//		//cout << "signal2 : " << signal2[0] << " " << signal2[1] << " " << signal2[2] << " " << signal2[3] << endl;
//
//		/*for(int ii=0; ii<ns;ii++)
//			if(signal[ii]!=signal2[ii]){
//				cout << signal[ii] << " != " << signal2[ii] << endl;
//				getchar();
//			}*/
//
//
//		///////////////// Procedure de test pour read_data_from_fits /////////////
//
//
//		long ns2;
//
//		// for each scan
//		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
//
//			//	cout << "[" << rank << "] " << idet << "/" << ndet << " " << iframe << "/" << ntotscan << endl;
//
//			// read pointing files
//			ns = nsamples[iframe];
//			//ff = fframes[iframe];
//
//			//cout << "iframe " << iframe << endl;
//
//			// lis les données RA/DEC, Phi, et scerr (sky coordinates err, BLASPEC) pour ff et ns
//
//			//			read_data_std(dirfile, ff, 0, ns, ra,   ra_field,  type); // type = 'd' 64-bit double
//			//			cout << "before " << ra[0] << " " << ra[1] << endl;
//			//
//			//			read_data_std(dirfile, ff, 0, ns, dec,  dec_field, type);
//			//			cout << "before dec " << dec[0] << " " << dec[1] << endl;
//			//
//			//			read_data_std(dirfile, ff, 0, ns, phi,  phi_field, type);
//			//			cout << "before phi " << phi[0] << " " << phi[1] << endl;
//
//
//			fits_file=fits_table[iframe];
//
//			read_position_from_fits(fits_file, ra, dec, phi, flpoint, 0, NULL, ns2, field);
//			/*cout << "after " << ra2[0] << " " << ra2[1] << endl << endl;
//			cout << "after dec " << dec2[0] << " " << dec2[1] << endl;
//			cout << "after phi " << phi2[0] << " " << phi2[1] << endl;
//
//			cout << "flpoint : " << flpoint2[0] <<  flpoint2[1] << flpoint2[2] << flpoint2[3] << endl;
//			cout << "ns : " << ns << " " << "ns 2 : " << ns2 << endl;*/
//
//			//getchar();
//
//			//read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, type);
//
//			//read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c'); // flpoint = donnee vs time, take the data or not
//
//			/*for (long ii=0;ii<ns;ii++)
//				if (isnan(ra[ii]) || isnan(dec[ii]) || isnan(phi[ii]))
//					flpoint[ii] = 1; // sample is flagged, don't take this sample*/
//
//
//
//			// find offset based on frame range
////			correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets); // TODO : A virer
//
//			if (default_projection){
//				// spheric coords to squared coords (tangent plane)
//				sph_coord_to_sqrmap(pixdeg, ra, dec, phi, offsets, ns, xx, yy, &nn, coordscorner,
//						tancoord, tanpix, bfixc, radius, /*offmap,*/ srccoord,0);
//			}else{
//				//compute ramin/max, etc... with WCS
//			}
//
//			// store ra/dec min and max (of the map) and for this processor
//			if (coordscorner[0] < ra_min) ra_min = coordscorner[0];
//			if (coordscorner[1] > ra_max) ra_max = coordscorner[1];
//			if (coordscorner[2] < dec_min) dec_min = coordscorner[2];
//			if (coordscorner[3] > dec_max) dec_max = coordscorner[3];
//
//		}
//
//	} //// end of idet loop
//
//	delete [] offsets;
//	delete [] froffsets;
//
//	//delete [] flpoint2;
//}
//
////TODO : remove : not used anymore
//void compute_seen_pixels_coordinates(long ntotscan,string outdir, std::vector<string> bolonames,  string *fits_table,/*,string bextension, string fextension, string termin,*/
//		/*string file_offsets,foffset *foffsets,float *scoffsets,*/ long iframe_min, long iframe_max,/*long *fframes,*/
//		long *nsamples,string dirfile,/*string ra_field,string dec_field,string phi_field, string scerr_field,
//		string flpoint_field,int nfoff*/ double pixdeg, int *&xx, int *&yy, unsigned short *&mask,int &nn, double *&coordscorner, double *tancoord,
//		double *tanpix, bool bfixc, double radius, /*double *offmap,*/ double *srccoord, char type, double *&ra,double *&dec,
//		double *&phi, short *&flpoint,int shift_data_to_point,double &ra_min,double &ra_max,double &dec_min,double &dec_max,short* &flag,
//		long napod, double errarcsec, bool NOFILLGAP,bool flgdupl, int factdupl,long addnpix,unsigned char *&rejectsamp, long *&samptopix, long *&pixon, int rank,
//		long *indpsrc, long npixsrc, int &flagon, bool &pixout){
//
//	/*!
//	 * \fn Get coordinates of pixels that are seen
//	 * Compute the position to pixel projetcion matrices :
//	 * One binary file per bolometer and per scan
//	 */
//
//	unsigned long ndet = bolonames.size();
//
//	string field;
//
//	string fits_file;
//	long ns, ns2;
//
//	int ll=0;
//
//	double * offsets;
//	offsets = new double[2]; //
//
//	/// loop again on detectors
//	//for (idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){
//	for (unsigned int idet=0;idet<ndet;idet++){
//
//
//		field = bolonames[idet];
//
//		// read bolometer offsets
//		//read_bolo_offsets(field,file_offsets,scoffsets,offsets);
//
//		fits_file=fits_table[0];
//
//		// TODO: Bolo offsets change per frames, so switch the loop between bolo and frame
//		read_bolo_offsets_from_fits(fits_file, field, offsets);
//
//		//loop to get coordinates of pixels that are seen
//		for (unsigned long iframe=0;iframe<ntotscan;iframe++){
//			ns = nsamples[iframe];
//			//ff = fframes[iframe];
//
//			fits_file=fits_table[iframe];
//
//			read_position_from_fits(fits_file, ra, dec, phi, flpoint, 1, flag, ns2, field);
//
//			for (long ii=0;ii<ns;ii++)
//				if (isnan(ra[ii]) || isnan(dec[ii]))
//					flpoint[ii] = 1;
//
//			if(ns!=ns2)
//				cerr << "Error. Fits file has a wrong number of samples\n";
//
//			// return coordscorner, nn ,xx, yy, tancoord, tanpix
//			sph_coord_to_sqrmap(pixdeg, ra, dec, phi, offsets, ns, xx, yy, &nn, coordscorner,
//					tancoord, tanpix, 1, radius, /*offmap,*/ srccoord,1);
//
////			for (long ii=0; ii<20; ii++){
////				cout << xx[ii] << " " << yy[ii] << endl;
////
////			}
//
//
//			// create flag field -----> conditions to flag data,
//			// flag
//			// scerr = pointing error : BLASPEC
//			// flpoint = do we consider those data at time t
//			// napod = number of samples to apodize
//			// errarcsec critere de flag, default = 15.0, BLAST specific : pointing error threshold
//			// NOFILLGAP = fill the gap ? default = yes
//			// rejectsamples: rejected samples array, rejectsample value =0,1,2 or 3
//			flag_conditions(flag,/*scerr,*/flpoint,ns,napod,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);
//			//flag_conditions(flag,scerr,flpoint,ns,napod,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);
//
//			// returns rejectsamp
//
//
//
//			for (long ii=0;ii<ns;ii++){
//				//sample is rejected
//				if (rejectsamp[ii] == 2){
//					pixout = 1; // pixel is out of map
//					pixon[factdupl*nn*nn + addnpix] += 1; // add one to a pixel outside the map
//					samptopix[ii] = factdupl*nn*nn + addnpix; // the ii sample corresponds to a pixel outside the map
//
//
//					printf("[%2.2i] PIXEL OUT, ii = %ld, xx = %i, yy = %i\n",rank, ii,xx[ii],yy[ii]);
//
//				}
//				// sample is not rejected
//				if (rejectsamp[ii] == 0) {
//					// if not in the  box constraint removal mask (defined by user)
//					if (mask[yy[ii]*nn + xx[ii]] == 1){
//						ll = yy[ii]*nn + xx[ii]; // ll=map indice
//					} else {//if pixel is in the mask
//						ll = indpsrc[yy[ii]*nn + xx[ii]]+factdupl*nn*nn+iframe*npixsrc;
//					}
//					pixon[ll] += 1;
//					samptopix[ii] = ll; // the ii sample ii coresponds to pixel ll
//
//					//printf("ll=%ld\n",ll);
//
//				}
//				// sample is flagged or scerr > scerr_threshold or flpoint =0
//				if (rejectsamp[ii] == 1){
//					if (flgdupl){ // if flagged pixels are in a duplicated map
//						ll = yy[ii]*nn + xx[ii];
//						pixon[nn*nn+ll] += 1;
//						samptopix[ii] = nn*nn+ll;
//					} else { // else every flagged sample is projected to the same pixel (outside the map)
//						pixon[nn*nn+1 + addnpix] += 1;
//						samptopix[ii] = nn*nn+1 + addnpix;
//					}
//				}
//				// pixel dans la zone d'apodisation //a suppr
//				if (rejectsamp[ii] == 3){
//					pixon[factdupl*nn*nn+1 + addnpix] += 1;
//					flagon = 1;
//					samptopix[ii] = factdupl*nn*nn+1 + addnpix;
//				}
//			}
//
//
//
//			//printf("pixon[nn*nn] = %d/n",pixon[nn*nn]);
//
//
//			//if (rank == 0){
//			write_samptopix(ns, samptopix, /*termin, */ outdir, idet, iframe);
//			//write_samptopix(ns, samptopix, termin, outdir, idet, index_table[iframe]);
//			//}
//
//		} // end of iframe loop
//
//	}//end of idet loop
//
//
//	delete [] offsets;
//
//}
//


void computePixelIndex(long ntotscan,string outdir, std::vector<string> bolonames,
		string *fits_table, long iframe_min, long iframe_max, long *nsamples,
		struct wcsprm & wcs, long NAXIS1, long NAXIS2,
		unsigned short *&mask,
		long napod, bool NOFILLGAP,bool flgdupl, int factdupl,
		long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout){

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
		fits_file=fits_table[iframe];
		ns = nsamples[iframe];

		double *ra, *dec, *phi, **offsets;
		short *flpoint;


		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		read_all_bolo_offsets_from_fits(fits_file, bolonames, offsets);

		// read reference position
		long test_ns;
		read_ReferencePosition_from_fits(fits_file, ra, dec, phi, flpoint, test_ns);

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

				celestial.ref[0] =  ra[ii]*15.;
				celestial.ref[1] =  dec[ii];
				celset(&celestial);

				//TODO : check this -1 factor... just a stupid convention...
				offxx = (cosphi[ii] * offsets[idet][0]
				                                - sinphi[ii] * offsets[idet][1])*-1;
				offyy =  sinphi[ii] * offsets[idet][0]
				                                + cosphi[ii] * offsets[idet][1];

				if (celx2s(&celestial, 1, 1, 0, 1, &offxx, &offyy, &lon, &lat, &ra_deg, &dec_deg, &status) == 1) {
					printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
					continue;
				}

				world[2*ii]   = ra_deg;
				world[2*ii+1] = dec_deg;
			}

			// ... and reproject it back onto the map
			if ((status = wcss2p(&wcs, ns, 2, world, phi, theta, imgcrd, pixcrd, wcsstatus))) {
				printf("   wcss2p(1) ERROR %2d \n", status);
				continue;
			}
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

			short *bolo_flag;

			long test_ns;
			read_flag_from_fits(fits_file, field, bolo_flag, test_ns);

			if (test_ns != ns) {
				cerr << "Read flags does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}


			for (long ii=0; ii<ns; ii++){

				if (flpoint[ii] == 1) bolo_flag[ii] = 1;

				if ((xx[ii] < 0)   || (yy[ii] < 0  ))          bolo_flag[ii] = 2;
				if ((xx[ii] >=  NAXIS1) || (yy[ii] >= NAXIS2)) bolo_flag[ii] = 2;
				if (NOFILLGAP && (bolo_flag[ii] == 1 ))        bolo_flag[ii] = 2;

				if ((ii < napod) || (ii >= ns-napod)) bolo_flag[ii] = 3;

			}

//			cout << "NOFILLGAP : " << NOFILLGAP << endl;
//			for (unsigned long ii=0; ii<ns; ii++)
//			cout << ii << " " << (NOFILLGAP && (bolo_flag[ii] == 1) << endl;
//				if (bolo_flag[ii]==1)
//					cout << ii << " bolo flag" << endl;

			// Compute sample pixel index according to its flag

			samptopix = new long long[ns];

			for (long ii=0 ; ii<ns; ii++){

				// image1        + flag pixel + out pixel + crossing constrain removal
				// NAXIS1*NAXIS2 + 1          + 1         + addnpix*nframe

				long long ll;
				switch (bolo_flag[ii]) {
				case 0:			// sample is not rejected

					ll = NAXIS1*yy[ii]+xx[ii];

					if (mask[ll] != 1) //if pixel is in the box constraint removal mask
						ll = factdupl*NAXIS1*NAXIS2 + 2 + iframe*npixsrc + indpsrc[ll];

					pixon[ll] += 1;
					samptopix[ii] = ll;

					break;

				case 1:	   // sample is flagged

					if (flgdupl){ // if flagged pixels are in a duplicated map
						ll = NAXIS1*yy[ii]+xx[ii]; // index in the second map...
						pixon[NAXIS1*NAXIS2 + ll] += 1;
						samptopix[ii] = NAXIS1*NAXIS2 + ll;
					} else { // else every flagged sample is projected to the same pixel (outside the map)
						ll = factdupl*NAXIS1*NAXIS2 + 1;
						pixon[ll] += 1;
						samptopix[ii] = ll;
					}

					break;

				case 2: // sample is rejected
					ll = factdupl*NAXIS1*NAXIS2 + 2;
					pixout = 1; // pixel is out of map
					pixon[ll] += 1;
					samptopix[ii] = ll;

					printf("[%2.2i] PIXEL OUT, ii = %ld, xx = %ld, yy = %ld\n",rank, ii,xx[ii],yy[ii]);

					break;
				case 3:	// apodized data -> flag
					ll = factdupl*NAXIS1*NAXIS2 + 1;
					flagon = 1;
					pixon[ll] += 1;
					samptopix[ii] = ll;
					break;
				default:
					break;
				}
			}

		write_samptopix(ns, samptopix,  outdir, idet, iframe);

		}//end of idet loop
	} // end of iframe loop
}

void computePixelIndex_HIPE(long ntotscan,string outdir, std::vector<string> bolonames,
		string *fits_table, long iframe_min, long iframe_max, long *nsamples,
		struct wcsprm & wcs, long NAXIS1, long NAXIS2,
		unsigned short *&mask,
		long napod, bool NOFILLGAP,bool flgdupl, int factdupl,
		long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout){

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
		fits_file=fits_table[iframe];
		ns = nsamples[iframe];

		double *ra, *dec;
		short *flag;

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

			read_ra_from_fits(fits_file, field, ra, test_ns);
			if (test_ns != ns) {
				cerr << "Read ra does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}
			read_dec_from_fits(fits_file, field, dec, test_ns);
			if (test_ns != ns) {
				cerr << "Read dec does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}
			read_flag_from_fits(fits_file, field, flag, test_ns);
			if (test_ns != ns) {
				cerr << "Read dec does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}
			for (long ii=0; ii <ns; ii++){
				//TODO: Very bad fix.... prevent from mapping flagged data to a seperate map
				if (flag[ii] == 0) {
					world[2*ii]   = ra[ii];
					world[2*ii+1] = dec[ii];
				} else {
					world[2*ii]   = wcs.crval[0];
					world[2*ii+1] = wcs.crval[1];
				}
			}
//			for (long ii=0; ii<ns; ii++){
//				cout << ii << " " << world[2*ii] << " " << world[2*ii+1] << endl;
//			}
			// ... and reproject it back onto the map
//			cout << "before reproject" << endl;

			if ((status = wcss2p(&wcs, ns, 2, world, phi, theta, imgcrd, pixcrd, wcsstatus))) {
				printf("   wcss2p(1) ERROR %2d \n", status);
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

			short *bolo_flag;

			read_flag_from_fits(fits_file, field, bolo_flag, test_ns);

			if (test_ns != ns) {
				cerr << "Read flags does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}


			for (long ii=0; ii<ns; ii++){

				//TODO: Bad Fix, check flag values and handle properly
				if (flag[ii] != 0) bolo_flag[ii] = 1;

				if ((xx[ii] < 0)   || (yy[ii] < 0  ))         bolo_flag[ii] = 2;
				if ((xx[ii] >=  NAXIS1) || (yy[ii] >=  NAXIS2)) bolo_flag[ii] = 2;
				if (NOFILLGAP && (bolo_flag[ii] == 1 ))       bolo_flag[ii] = 2;

				if ((ii < napod) || (ii >= ns-napod)) bolo_flag[ii] = 3;

			}

//			cout << "NOFILLGAP : " << NOFILLGAP << endl;
//			for (unsigned long ii=0; ii<ns; ii++)
//			cout << ii << " " << (NOFILLGAP && (bolo_flag[ii] == 1) << endl;
//				if (bolo_flag[ii]==1)
//					cout << ii << " bolo flag" << endl;

			// Compute sample pixel index according to its flag

			samptopix = new long long[ns];

			for (long ii=0 ; ii<ns; ii++){

				// image1        + flag pixel + out pixel + crossing constrain removal
				// NAXIS1*NAXIS2 + 1          + 1         + addnpix*nframe

				long long ll;
				switch (bolo_flag[ii]) {
				case 0:			// sample is not rejected

					ll = NAXIS1*yy[ii]+xx[ii];

					if (mask[ll] != 1) //if pixel is in the box constraint removal mask
						ll = factdupl*NAXIS1*NAXIS2 + 2 + iframe*npixsrc + indpsrc[ll];

					pixon[ll] += 1;
					samptopix[ii] = ll;

					break;

				case 1:	   // sample is flagged

					if (flgdupl){ // if flagged pixels are in a duplicated map
						ll = NAXIS1*yy[ii]+xx[ii]; // index in the second map...
						pixon[NAXIS1*NAXIS2 + ll] += 1;
						samptopix[ii] = NAXIS1*NAXIS2 + ll;
					} else { // else every flagged sample is projected to the same pixel (outside the map)
						ll = factdupl*NAXIS1*NAXIS2 + 1;
						pixon[ll] += 1;
						samptopix[ii] = ll;
					}

					break;

				case 2: // sample is rejected
					ll = factdupl*NAXIS1*NAXIS2 + 2;
					pixout = 1; // pixel is out of map
					pixon[ll] += 1;
					samptopix[ii] = ll;

					printf("[%2.2i] PIXEL OUT, ii = %ld, xx = %ld, yy = %ld\n",rank, ii,xx[ii],yy[ii]);

					break;
				case 3:	// apodized data -> flag
					ll = factdupl*NAXIS1*NAXIS2 + 1;
					flagon = 1;
					pixon[ll] += 1;
					samptopix[ii] = ll;
					break;
				default:
					break;
				}
			}

		write_samptopix(ns, samptopix,  outdir, idet, iframe);

		}//end of idet loop
	} // end of iframe loop
}



//TODO: remove : not used
/*
int map_offsets(string file_frame_offsets,long ntotscan, float *&scoffsets, foffset *&foffsets,long *fframes, int rank){

	int ff, nfoff;

	if (file_frame_offsets != "NOOFFS") // file_frame_offset in ini file, containing offsets for different frame ranges
		foffsets = read_mapoffsets(file_frame_offsets, scoffsets, &nfoff); //-> return foffsets, scoffsets, nfoff
	else {
		printf("[%2.2i] No offsets between visits\n", rank);
		nfoff = 2;
		ff = fframes[0];
		for (int ii=0;ii<ntotscan;ii++) if (fframes[ii] > ff) ff = fframes[ii];
		foffsets = new foffset [2];
		(foffsets[0]).frame = 0;
		(foffsets[0]).pitch = 0.0;
		(foffsets[0]).yaw   = 0.0;
		(foffsets[1]).frame = ff+1; // ff = last frame number here
		(foffsets[1]).pitch = 0.0;
		(foffsets[1]).yaw   = 0.0;
		for (int ii=0;ii<6;ii++)
			scoffsets[ii] = 0.0;
	}
	return nfoff;
}
*/
