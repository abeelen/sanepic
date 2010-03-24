

#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}



using namespace std;

void file_compatibility_verification(struct samples samples_struct){

	long numscan = samples_struct.ntotscan;
	fitsfile * fptr1;
	fitsfile * fptr2;
	int status = 0;

	for(long ii = 1; ii< numscan; ii++ ){

		double *time_one, *time_two;
		struct detectors det_one, det_two;
		// TODO : Do we have to test that offsets are the same ??
		string file1, file2;

		file1 = samples_struct.fitsvect[ii-1];
		file2 = samples_struct.fitsvect[ii];

		// open both files
		if (fits_open_file(&fptr1, file1.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		if (fits_open_file(&fptr2, file2.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		long ns1, ns2;
		ns1=samples_struct.nsamples[ii-1];
		ns2=samples_struct.nsamples[ii];

		time_one=new double [samples_struct.nsamples[ii-1]];
		time_two=new double [samples_struct.nsamples[ii]];

		read_time_from_fits(file1, time_one, ns1);
		read_time_from_fits(file2, time_two, ns2);

		if(time_one[ns1-1]>=time_two[0]){
			cout << "There is a problem between " << file1 << " and " << file2 << " time tables.\n";
			cout << "Please order them correctly in order to have a time ordered list of input files. Exiting...\n";
			exit(EXIT_FAILURE);
		}

		read_bolo_list(file1, det_one);
		read_bolo_list(file2, det_two);

		if(det_one.ndet!=det_two.ndet){
			cout << "Error ! The channel lists in " << file1 << " and " << file2 << " are not the same. Exiting... \n";
			exit(EXIT_FAILURE);
		}

		for(long jj=0; jj<det_one.ndet; jj++)
			if(det_one.boloname[jj]!=det_two.boloname[jj]){
				cout << "Error ! The channel lists in " << file1 << " and " << file2 << " are not the same.\n";
				cout << det_one.boloname[jj] << " != " << det_two.boloname[jj] << ". Exiting...\n";
				exit(EXIT_FAILURE);
			}


		delete [] time_one;
		delete [] time_two;
		det_one.boloname.clear();
		det_two.boloname.clear();
		det_one.ndet=0;
		det_two.ndet=0;

		if (fits_close_file(fptr1, &status))
			fits_report_error(stderr, status);
		if (fits_close_file(fptr2, &status))
			fits_report_error(stderr, status);


	}

}

void copy_ref_pos(fitsfile *outfptr, struct samples samples_struct, long ns_final){


	fitsfile * fptr;
	long indice_debut=0;
	long indice_fin=0;
	double *RA_bis, *DEC_bis, *PHI_bis;

	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];
	PHI_bis = new double [ns_final];

	int status=0;

	for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){

		long ns_temp=0;


		string fname=samples_struct.fitsvect[iframe];
		indice_debut = indice_fin;


		// open each fits file and copy each table in the final file
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		if(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status))
			fits_report_error(stderr, status);

		if(iframe==0){
			if(fits_copy_header(fptr, outfptr, &status))
				fits_report_error(stderr, status);
			if(fits_movnam_hdu(outfptr, BINARY_TBL, (char*) "reference position", NULL, &status))
				fits_report_error(stderr, status);
			if(fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_final, (char*)"Number of rows", &status))
				fits_report_error(stderr, status);
		}

		double *RA, *DEC, *PHI;

		read_ReferencePosition_from_pointer(fptr, RA, DEC, PHI, ns_temp);


		for(long ii = 0; ii< ns_temp; ii++){
			RA_bis[indice_debut+ii]=RA[ii]*15.0;
			DEC_bis[indice_debut+ii]=DEC[ii];
			PHI_bis[indice_debut+ii]=PHI[ii];
		}


		indice_fin = indice_debut + samples_struct.nsamples[iframe];

		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		delete [] RA;
		delete [] DEC;
		delete [] PHI;
	}

	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, RA_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_final, DEC_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_final, PHI_bis, &status);

	delete [] RA_bis;
	delete [] DEC_bis;
	delete [] PHI_bis;

}

void copy_offsets(fitsfile * fptr, fitsfile *outfptr){

	int status=0;

	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);

	for(int col=1;col<4;col++)
		fits_copy_col(fptr, outfptr,  col, col,	0, &status);

}


void copy_channels(fitsfile * fptr, fitsfile *outfptr){

	int status=0;


	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);


	fits_copy_col(fptr, outfptr,  1, 1,	0, &status);

}

void copy_time(fitsfile *outfptr, struct samples samples_struct, long ns_final){


	fitsfile * fptr;
	long indice_debut=0;
	long indice_fin=0;
	double *time_bis;

	time_bis = new double [ns_final];

	int status=0;

	for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){

		long ns_temp=samples_struct.nsamples[iframe];


		string fname=samples_struct.fitsvect[iframe];
		indice_debut = indice_fin;


		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);


		fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status);
		if(iframe==0){
			fits_copy_header(fptr, outfptr, &status);
			fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);
		}

		double *time;
		read_time_from_fits(fname, time, ns_temp);

		for(long ii = 0; ii< ns_temp; ii++){
			time_bis[indice_debut + ii]=time[ii];
		}
		indice_fin = indice_debut + samples_struct.nsamples[iframe];
		delete [] time;
	}
	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status);


	delete [] time_bis;

}

void copy_signal(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final){

	fitsfile * fptr;
	long indice_debut=0;
	long indice_fin=0;
	double  *signal_bis;
	signal_bis = new double [ns_final];

	for(long jj=0;jj<det.ndet;jj++){
		int status=0;

		for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){

			long ns_temp=0;


			string fname=samples_struct.fitsvect[iframe];
			//			cout << fname << endl;
			indice_debut = indice_fin;


			// open each fits file and copy each table in the final file
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);

			fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status);

			if((jj==0)&&(iframe==0)){
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);
			}

			double *signal;

			read_signal_from_fits(fname, det.boloname[jj], signal, ns_temp);
			for(long ii = 0; ii< ns_temp; ii++)
				signal_bis[indice_debut+ii]=signal[ii];


			indice_fin = indice_debut + samples_struct.nsamples[iframe];
			delete [] signal;

			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

		}
		string field= det.boloname[jj];
		string fname=samples_struct.fitsvect[0];
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		long rowIndex = find_channel_index(fptr, field.c_str());
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);
		long fpixel[2]={1,rowIndex};
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, signal_bis, &status);
		indice_debut=0;
		indice_fin=0;
	}
	delete [] signal_bis;

}

void copy_mask(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final){

	fitsfile * fptr;
	long indice_debut=0;
	long indice_fin=0;
	double  *mask_bis;
	mask_bis = new double [ns_final];

	for(long jj=0;jj<det.ndet;jj++){
		int status=0;

		for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){

			long ns_temp=0;
			string fname=samples_struct.fitsvect[iframe];
			indice_debut = indice_fin;


			// open each fits file and copy each table in the final file
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);

			fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status);

			if((jj==0)&&(iframe==0)){
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);
			}

			double *mask;

			read_signal_from_fits(fname, det.boloname[jj], mask, ns_temp);
			for(long ii = 0; ii< ns_temp; ii++)
				mask_bis[indice_debut+ii]=mask[ii];


			indice_fin = indice_debut + samples_struct.nsamples[iframe];
			delete [] mask;

			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

		}

		string field= det.boloname[jj];
		string fname=samples_struct.fitsvect[0];
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		long rowIndex = find_channel_index(fptr, field.c_str());

		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		long fpixel[2]={1,rowIndex};
		fits_write_pix(outfptr, TINT, fpixel, ns_final, mask_bis, &status);
		indice_debut=0;
		indice_fin=0;
	}
	delete [] mask_bis;
}


void copy_RA_DEC(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final){

	fitsfile * fptr;
	long indice_debut=0;
	long indice_fin=0;
	double *RA_bis;
	double *DEC_bis;
	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];

	for(long jj=0;jj<det.ndet;jj++){
		int status=0;

		for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){

			long ns_temp=0;


			string fname=samples_struct.fitsvect[iframe];
			indice_debut = indice_fin;


			// open each fits file and copy each table in the final file
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);


			if((jj==0)&&(iframe==0)){
				fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status);
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

				fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status);
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);
			}

			double *RA, *DEC;
			read_ra_dec_from_fits(fname, det.boloname[jj], RA, DEC, ns_temp);

			for(long ii = 0; ii< ns_temp; ii++){
				RA_bis[indice_debut + ii]=RA[ii];
				DEC_bis[indice_debut + ii]=DEC[ii];
			}



			indice_fin = indice_debut + samples_struct.nsamples[iframe];

			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

			delete [] RA;
			delete [] DEC;
		}

		string field= det.boloname[jj];
		string fname=samples_struct.fitsvect[0];
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		long rowIndex = find_channel_index(fptr, field.c_str());
		long fpixel[2]={1,rowIndex};
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "ra", NULL, &status);
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, RA_bis, &status);
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "dec", NULL, &status);
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, DEC_bis, &status);
		indice_debut=0;
		indice_fin=0;
	}
	delete [] RA_bis;
	delete [] DEC_bis;


}






