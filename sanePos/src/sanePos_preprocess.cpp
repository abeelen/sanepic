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


void find_coordinates_in_map(long ndet,std::vector<string> bolonames,string bextension,
		string fextension,string file_offsets,foffset *foffsets,float *scoffsets,	/*double *offsets,*/long iframe_min, long iframe_max,
		long *fframes,long *nsamples,string dirfile,string ra_field,string dec_field,string phi_field, string scerr_field,
		string flpoint_field,int nfoff,double pixdeg, int *&xx, int *&yy,int nn, double *&coordscorner, double *tancoord,
		double *tanpix, bool bfixc, double radius, double *offmap, double *srccoord,char type,double *&ra,double *&dec,
		double *&phi,double *&scerr, short *&flpoint,double &ra_min,double &ra_max,double &dec_min,double &dec_max,bool default_projection){

	double * offsets,*froffsets;
	offsets = new double[2]; //
	froffsets = new double[2]; //

	string field;
	string bolofield;
	string flagfield;
	long ns,ff;
	//Sdouble *ra,*dec,*phi,*scerr, *flpoint;

	// ndet = number of channels
	for (int idet=0;idet<ndet;idet++){
		//cout << idet << endl;
		// field = actual boloname
		field = bolonames[idet];

		// bolofield = boloname + bextension
		bolofield = field+bextension;

		//if (cextension != "NOCALP")
		//	calfield  = field+cextension;
		if (fextension != "NOFLAG")
			flagfield = field+fextension;

		//read bolometer offsets
		read_bolo_offsets(field,file_offsets,scoffsets,offsets);
		cout << "before " << offsets[0] << " " << offsets[1] << endl;
		read_bolo_offsets_from_fits("data_00.fits", field, offsets);
		cout << "after  " << offsets[0] << " " << offsets[1] << endl;

		long ns2;

		// for each scan
		for (long iframe=iframe_min;iframe<iframe_max;iframe++){

			//	cout << "[" << rank << "] " << idet << "/" << ndet << " " << iframe << "/" << ntotscan << endl;

			// read pointing files
			ns = nsamples[iframe];
			ff = fframes[iframe];

			cout << "Length  : " << readFitsLength("data_00.fits") << endl;


			cout << "iframe " << iframe << endl;
			// lis les données RA/DEC, Phi, et scerr (sky coordinates err, BLASPEC) pour ff et ns
//			read_data_std(dirfile, ff, 0, ns, ra,   ra_field,  type); // type = 'd' 64-bit double
//			cout << "before " << ra[0] << " " << ra[1] << endl;
//
//			read_data_std(dirfile, ff, 0, ns, dec,  dec_field, type);
//			cout << "before dec " << dec[0] << " " << dec[1] << endl;
//
//			read_data_std(dirfile, ff, 0, ns, phi,  phi_field, type);
//			cout << "before phi " << phi[0] << " " << phi[1] << endl;

			short *flpoint2;
			flpoint2 = new short[ns];

			read_position_from_fits("data_00.fits", ra, dec, phi, flpoint2, ns2);
			cout << "after " << ra[0] << " " << ra[1] << endl << endl;
			cout << "after dec " << dec[0] << " " << dec[1] << endl;
			cout << "after phi " << phi[0] << " " << phi[1] << endl;
			cout << "flpoint : " << flpoint2[0] <<  flpoint2[1] << flpoint2[2] << flpoint2[3] << endl;
			cout << "ns : " << ns << " " << "ns 2 : " << ns2 << endl;

			getchar();
			//read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, type); // TODO : ca n'existe plus

			//read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c'); // flpoint = donnee vs time, take the data or not

			for (long ii=0;ii<ns;ii++)
				if (isnan(ra[ii]) || isnan(dec[ii]) || isnan(phi[ii]))
					flpoint[ii] = 1; // sample is flagged, don't take this sample



			// find offset based on frame range
			//correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets); // TODO : A virer

			if (default_projection){
				// spheric coords to squared coords (tangent plane)
				sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
						tancoord, tanpix, bfixc, radius, offmap, srccoord,0);
			}else{
				//compute ramin/max, etc... with WCS
			}

			// store ra/dec min and max (of the map) and for this processor
			if (coordscorner[0] < ra_min) ra_min = coordscorner[0];
			if (coordscorner[1] > ra_max) ra_max = coordscorner[1];
			if (coordscorner[2] < dec_min) dec_min = coordscorner[2];
			if (coordscorner[3] > dec_max) dec_max = coordscorner[3];

		}

	} //// end of idet loop

	delete [] offsets;
	delete [] froffsets;
}



long Compute_indpsrc_addnpix(int nn, long ntotscan, std::vector<long> xxi, std::vector<long> xxf, std::vector<long> yyi,
		std::vector<long> yyf, long* &indpsrc,long &npixsrc,short *&mask){

	long addnpix=0;

	//************************************* Deal with masking the point sources
	// define the mask



	for (long ii=0;ii<nn*nn;ii++){
		mask[ii] = 1;
		indpsrc[ii] = -1;
	}


	if (xxi.size() != 0){
		for (long ib = 0;ib < (long)xxi.size(); ib++){ // to avoid warning, mat-27/05
			// for each box crossing constraint removal
			for (long ii=xxi[ib];ii<xxf[ib];ii++)
				for (long ll=yyi[ib];ll<yyf[ib];ll++){
					mask[ll*nn + ii] = 0;  // mask is initialised to 0
					indpsrc[ll*nn + ii] = npixsrc;
					npixsrc += 1;
				}
		}
	}





	/*for (long ii=0;ii<nn*nn;ii++){
		if (mask[ii] == 0){
			indpsrc[ii] = npixsrc;
			npixsrc += 1;
		} else {
			indpsrc[ii] = -1;
		}
	}*/
	addnpix = ntotscan*npixsrc; // addnpix = number of pix to add in pixon = number of scans * number of pix in box crossing constraint removal

	//debug
	cout  << "addnpix : " << addnpix << endl;



	return addnpix;
}

void compute_seen_pixels_coordinates(int ndet,long ntotscan,string outdir, std::vector<string> bolonames,string bextension, string fextension,string termin,
		string file_offsets,foffset *foffsets,float *scoffsets,long iframe_min, long iframe_max,long *fframes,
		long *nsamples,string dirfile,string ra_field,string dec_field,string phi_field, string scerr_field,
		string flpoint_field,int nfoff,double pixdeg, int *&xx, int *&yy, short *&mask,int &nn, double *&coordscorner, double *tancoord,
		double *tanpix, bool bfixc, double radius, double *offmap, double *srccoord, char type, double *&ra,double *&dec,
		double *&phi,double *&scerr, short *&flpoint,int shift_data_to_point,double &ra_min,double &ra_max,double &dec_min,double &dec_max,short* &flag,
		long napod, double errarcsec, bool NOFILLGAP,bool flgdupl, int factdupl,long addnpix,short *&rejectsamp, long *&samptopix, long *&pixon, int rank,
		long *indpsrc, long npixsrc, int &flagon, bool &pixout){


	string field;
	string bolofield;
	string flagfield;
	long ns,ff, ns2;
	short *flpoint2;


	int ll=0;

	double * offsets,*froffsets;
	offsets = new double[2]; //
	froffsets = new double[2]; //

	/// loop again on detectors
	//for (idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){
	for (int idet=0;idet<ndet;idet++){


		field = bolonames[idet];
		//printf("%s\n",field.c_str());
		bolofield = field+bextension;
		//calfield  = field+cextension;
		flagfield = field+fextension;


		// read bolometer offsets
		read_bolo_offsets(field,file_offsets,scoffsets,offsets);

		read_bolo_offsets_from_fits("data_00.fits", field, offsets);


		//loop to get coordinates of pixels that are seen
		for (long iframe=0;iframe<ntotscan;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
//
//			read_data_std(dirfile, ff, 0, ns, ra,   ra_field,  type);
//			read_data_std(dirfile, ff, 0, ns, dec,  dec_field, type);
//			read_data_std(dirfile, ff, 0, ns, phi,  phi_field, type);
//			//read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, type);
//			read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
//			for (long ii=0;ii<ns;ii++)
//				if (isnan(ra[ii]) || isnan(dec[ii]))
//					flpoint[ii] = 1;
//
//			flpoint2 = new short[ns];

			// TODO : flag should be read here, all time, and define as short int
			read_position_from_fits("data_00.fits", ra, dec, phi, flpoint, ns2);
//
//			if (fextension != "NOFLAG"){
//				read_data_std(dirfile, ff, shift_data_to_point, ns, flag, flagfield,  'c');
//				read_data_from_fits("data_00.fits", ra, dec, phi, flpoint2, 1, flag, ns2, field);
//			} else {
//				read_data_from_fits("data_00.fits", ra, dec, phi, flpoint2, 0, NULL, ns2, field);
//				for (long ii=0;ii<ns;ii++)
//					flag[ii] = 0;
//			}



			// find offset based on frame range
			//correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets);



			// return coordscorner, nn ,xx, yy, tancoord, tanpix
			sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
					tancoord, tanpix, 1, radius, offmap, srccoord,1);




			// create flag field -----> conditions to flag data,
			// flag
			// scerr = pointing error : BLASPEC
			// flpoint = do we consider those data at time t
			// napod = number of samples to apodize
			// errarcsec critere de flag, default = 15.0, BLAST specific : pointing error threshold
			// NOFILLGAP = fill the gap ? default = yes
			// rejectsamples: rejected samples array, rejectsample value =0,1,2 or 3
			flag_conditions(flag,scerr,flpoint,ns,napod,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);
			//flag_conditions(flag,scerr,flpoint,ns,napod,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);

			// returns rejectsamp



			for (long ii=0;ii<ns;ii++){
				//sample is rejected
				if (rejectsamp[ii] == 2){
					pixout = 1; // pixel is out of map
					pixon[factdupl*nn*nn + addnpix] += 1; // add one to a pixel outside the map
					samptopix[ii] = factdupl*nn*nn + addnpix; // the ii sample corresponds to a pixel outside the map


					printf("[%2.2i] PIXEL OUT, ii = %ld, xx = %i, yy = %i\n",rank, ii,xx[ii],yy[ii]);

				}
				// sample is not rejected
				if (rejectsamp[ii] == 0) {
					// if not in the  box constraint removal mask (defined by user)
					if (mask[yy[ii]*nn + xx[ii]] == 1){
						ll = yy[ii]*nn + xx[ii]; // ll=map indice
					} else {//if pixel is in the mask
						ll = indpsrc[yy[ii]*nn + xx[ii]]+factdupl*nn*nn+iframe*npixsrc;
					}
					pixon[ll] += 1;
					samptopix[ii] = ll; // the ii sample ii coresponds to pixel ll

					//printf("ll=%ld\n",ll);

				}
				// sample is flagged or scerr > scerr_threshold or flpoint =0
				if (rejectsamp[ii] == 1){
					if (flgdupl){ // if flagged pixels are in a duplicated map
						ll = yy[ii]*nn + xx[ii];
						pixon[nn*nn+ll] += 1;
						samptopix[ii] = nn*nn+ll;
					} else { // else every flagged sample is projected to the same pixel (outside the map)
						pixon[nn*nn+1 + addnpix] += 1;
						samptopix[ii] = nn*nn+1 + addnpix;
					}
				}
				// pixel dans la zone d'apodisation //a suppr
				if (rejectsamp[ii] == 3){
					pixon[factdupl*nn*nn+1 + addnpix] += 1;
					flagon = 1;
					samptopix[ii] = factdupl*nn*nn+1 + addnpix;
				}
			}



			//printf("pixon[nn*nn] = %d/n",pixon[nn*nn]);


			//if (rank == 0){
			write_samptopix(ns, samptopix, termin, outdir, idet, iframe);
			//}

		} // end of iframe loop

	}//end of idet loop


	delete [] offsets;
	delete [] froffsets;
}


int compute_indpix(long *&indpix,int factdupl,int nn,long addnpix,long *pixon){

	int npix;
	int ll=0;


	// initialize indpix
	fill(indpix, indpix+(factdupl*nn*nn+2 + addnpix),-1);



	for(int ii=0;ii<factdupl*nn*nn+2 + addnpix;ii++){
		if (pixon[ii] != 0){
			indpix[ii] = ll; // pixel indice : 0 -> number of seen pixels
			ll++;
		}
	}
	npix = ll;  // npix = number of filled pixelsx;

	return npix;
}


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
