
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "SanePosPreProcess.h"
#include "SanePosMapMaking.h"
#include "DataIO.h"
#include "TemporaryIO.h"
#include "InputFileIO.h"

extern "C" {
#include "getdata.h"
}



using namespace std;



int modify_mask_flag_in_dirfile(struct samples samples_struct, std::vector<std::vector<std::string> > bolo_list, long long *indpsrc,
		long NAXIS1, long NAXIS2){


	long long * samptopix;
	int *mask;

	for(long iframe=samples_struct.iframe_min; iframe< samples_struct.iframe_max; iframe++){
		// set dirfile name and binary name
		string filename = samples_struct.basevect[iframe];
		long ns = samples_struct.nsamples[iframe];

		mask  = new int[ns];
		samptopix = new long long[ns];

		string output_read = "";
		std::vector<string> bolonames = bolo_list[iframe];
		//		if(read_channel_list(output_read, samples_struct.bolovect[iframe], bolonames)){
		//			cout << output_read << endl;
		//			return 1;
		//		}

		long ndet = (long)bolonames.size();

		for(long idet=0; idet < ndet; idet ++){


			string field = bolonames[idet];
			string flag_outfile = "flag_" + filename + "_" + field;

			if(readFlagFromDirfile(samples_struct.dirfile_pointers[iframe], filename, field, mask, ns))
				return 1;

			if(readSampToPix(samples_struct.dirfile_pointers[iframe], filename, field, samptopix, ns))
				return 1;

			// modify if in indpsrc !
			for(long ii=0; ii< ns ; ii++)
				if(indpsrc[samptopix[ii]]!=-1)
					mask[ii]=-1;

			//put back in dirfile
			int n_write = gd_putdata(samples_struct.dirfile_pointers[iframe], (char*) flag_outfile.c_str(), 0, 0, 0, ns, GD_INT32,
					mask);
			if (gd_error(samples_struct.dirfile_pointers[iframe]) != 0) {
				cout << "error putdata in modify_mask_flag_in_dirfile : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			gd_flush(samples_struct.dirfile_pointers[iframe], NULL);
		}

		delete [] mask;
		delete [] samptopix;

	}

	return 0;
}

// deprecated
int computePixelIndex(string dirfile, struct samples samples_struct, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout, std::vector<std::vector<std::string> > bolo_vect)
{

	long long  *samptopix;

	string field;
	string fits_file, base_file;

	long ns;

	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){
		// for each scan
		fits_file=samples_struct.fitsvect[iframe];
		base_file = samples_struct.basevect[iframe];
		ns = samples_struct.nsamples[iframe];

		std::vector<string> det_vect = bolo_vect[iframe];

		long ndet = (long)det_vect.size();


		double *lon, *lat, *phi, **offsets;


		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		if(read_all_bolo_offsets_from_fits(fits_file, det_vect, offsets))
			return 1;

		// read reference position
		long test_ns;
		if(read_ReferencePosition_from_fits(fits_file, lon, lat, phi, test_ns))
			return 1;

		if (test_ns != ns) {
			cerr << "Read position does not correspond to frame size : Check !" << endl;
			return 1;
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

		// First compute the sin and cos of the phi (will be reused ndet times..)
		double *cosphi, *sinphi;
		cosphi = new double[ns];
		sinphi = new double[ns];

		for (long ii=0; ii<ns; ii++){
			cosphi[ii] = cos(phi[ii]/180.0*M_PI);
			sinphi[ii] = sin(phi[ii]/180.0*M_PI);
		}

		delete [] phi;

		for (long idet = 0; idet<ndet; idet++){

			string field = det_vect[idet];

			double offxx, offyy, dummy1, dummy2, lon_deg, lat_deg;
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

				celestial.ref[0] =  lon[ii];
				celestial.ref[1] =  lat[ii];
				celset(&celestial);

				//TODO : check this -1 factor... just a stupid convention...
				offxx = (cosphi[ii] * offsets[idet][0]
				                                    - sinphi[ii] * offsets[idet][1])*-1;
				offyy =  sinphi[ii] * offsets[idet][0]
				                                    + cosphi[ii] * offsets[idet][1];

				if (celx2s(&celestial, 1, 0, 0, 0, &offxx, &offyy, &dummy1, &dummy2, &lon_deg, &lat_deg, &status) == 1) {
					printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
					continue;
				}

				world[2*ii]   = lon_deg;
				world[2*ii+1] = lat_deg;
			}


			// ... and reproject it back onto the map
			if ((status = wcss2p(wcs, ns, 2, world, phi, theta, imgcrd, pixcrd, wcsstatus))) {
				printf("   wcss2p ERROR %2d \n", status);
				for (long ii = 0; ii < ns ; ii++){
					cout << ii << " " << world[2*ii] << " " << world[2*ii+1] << " : " << phi[ii] << " " << theta[ii] << " : " << pixcrd[2*ii] << " " << pixcrd[2*ii+1] << " : " << wcsstatus[ii] << endl;

				}
				exit(1);
			}
			// Transform pixel coordinates to pixel index
			for (long ii = 0; ii < ns; ii++){
				xx[ii]  = (long long) (pixcrd[2*ii]-0.5);
				yy[ii]  = (long long) (pixcrd[2*ii+1]-0.5);
			}

			delete [] world;
			delete [] imgcrd;
			delete [] pixcrd;
			delete [] theta;
			delete [] phi;
			delete [] wcsstatus;

			int *bolo_flag;

			bolo_flag = new int[ns];

			if(readFlagFromDirfile(samples_struct.dirfile_pointers[iframe], base_file, field, bolo_flag, ns))
				return 1;



			for (long ii=0; ii<ns; ii++){

				if (bolo_flag[ii] != 0) bolo_flag[ii] = 1;

				if ((xx[ii] < 0)   || (yy[ii] < 0  ))          bolo_flag[ii] = 2;
				if ((xx[ii] >=  NAXIS1) || (yy[ii] >= NAXIS2)) bolo_flag[ii] = 2;
				if ( (! proc_param.fill_gap ) && (bolo_flag[ii] == 1 ))        bolo_flag[ii] = 2;

				if ((ii < proc_param.napod) || (ii >= ns-proc_param.napod)) bolo_flag[ii] = 3;

			}

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

			if(writeSampToPix(samples_struct.dirfile_pointers[iframe], base_file, det_vect[idet], ns, samptopix))
				return 1;

			delete [] bolo_flag;
			delete [] samptopix;
			delete [] xx;
			delete [] yy;

		}//end of idet loop
		delete [] lon;
		delete [] lat;
		free_dmatrix(offsets,(long)0,ndet-1,(long)0,2-1);
		delete [] cosphi;
		delete [] sinphi;

	} // end of iframe loop

	return 0;
}

//int computeMaskIntersection(samples_struct, pos_param, iframe_min, iframe_max, rank,
//				wcs, NAXIS1, NAXIS2, mask, indpsrc, npixsrc, bolo_list, count_pixsrc)
//{
//
//
//}

int computePixelIndex_HIPE(struct samples samples_struct, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *indpsrc, long long npixsrc,
		long long *&pixon, int &flagon, int &pixout,
		int sub_rank, int sub_size)
{

	long long  *samptopix;

	string field;

	string fits_file, base_file;
	long ns;
	int status;

	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){
		// for each scan
		fits_file=samples_struct.fitsvect[iframe];
		base_file= samples_struct.basevect[iframe];
		ns = samples_struct.nsamples[iframe];

		std::vector<string> det_vect = samples_struct.bolo_list[iframe];

		long ndet = (long)det_vect.size();

		double *lon, *lat;
		int *flag, *bolo_flag;

		lon  = new double[ns];
		lat  = new double[ns];
		flag = new int[ns];
		bolo_flag = new int[ns];

		for (long idet = floor(sub_rank*ndet*1.0/sub_size); idet<floor((sub_rank+1)*ndet*1.0/sub_size); idet++){

			string field = det_vect[idet];

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


			if(readLonFromDirfile(samples_struct.dirfile_pointers[iframe], base_file, field, lon, ns))
				return 1;
			if(readLatFromDirfile(samples_struct.dirfile_pointers[iframe], base_file, field, lat, ns))
				return 1;

			if(readFlagFromDirfile(samples_struct.dirfile_pointers[iframe], base_file, field, flag, ns))
				return 1;
			for (long ii=0; ii <ns; ii++){
				world[2*ii]   = lon[ii];
				world[2*ii+1] = lat[ii];
			}

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


			delete [] world;
			delete [] imgcrd;
			delete [] pixcrd;
			delete [] theta;
			delete [] phi;
			delete [] wcsstatus;


			// Combine position and bolo flags
			// and check

			if(readFlagFromDirfile(samples_struct.dirfile_pointers[iframe], base_file, field, bolo_flag, ns))
				return 1;

			for (long ii=0; ii<ns; ii++){

				//TODO: Bad Fix, check flag values and handle properly
				if (flag[ii] != 0)                               bolo_flag[ii] = 1;
				if (flag[ii] == -1)                              bolo_flag[ii] = 2;

				if ((xx[ii] < 0)   || (yy[ii] < 0  ))            bolo_flag[ii] = 2;
				if ((xx[ii] >=  NAXIS1) || (yy[ii] >=  NAXIS2))  bolo_flag[ii] = 2;
				if ( (! proc_param.fill_gap ) && (bolo_flag[ii] == 1 ))        bolo_flag[ii] = 2;

				if ((ii < proc_param.napod) || (ii >= ns-proc_param.napod))    bolo_flag[ii] = 3;

			}

			// Compute sample pixel index according to its flag

			samptopix = new long long[ns];

			for (long ii=0 ; ii<ns; ii++){

				// image1        + (if 2nd image)  + crossing constrain removal + flagged pixel + apodized data
				// NAXIS1*NAXIS2 + NAXIS1*NAXIS2   + npixsrc*nframe             + 1             + 1
				//                                   addnpix


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
					printf("II - %s PIXEL OUT, ii = %ld, xx = %lld, yy = %lld\n",field.c_str(), ii,xx[ii],yy[ii]);
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

			if(writeSampToPix(samples_struct.dirfile_pointers[iframe], base_file, det_vect[idet], ns, samptopix))
				return 1;

			delete [] samptopix;
			delete [] xx;
			delete [] yy;


		}
		delete [] lon;
		delete [] lat;
		delete [] flag;
		delete [] bolo_flag;

	}

	return 0;

}
