

#include "covMatrixIO.h"
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



//void read_Split_file(string fname, std::vector< long > &cut_sample, struct samples sample_struct){
//
//	std::ifstream file;
//	string s, line;
//	int nb_elem=0;
//
//	file.open(fname.c_str(), ios::in);
//	if(!file.is_open()){
//		cerr << "File [" << fname << "] Invalid." << endl;
//		exit(-1);
//	}
//
//	nb_elem=0;
//
//
//	cout << " while \n";
//	while(file >> s){
//		size_t found;
//		s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
//		found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
//		if (found == 0) continue;
//		cut_sample.push_back(atol(s.c_str()));
//		nb_elem++;
//		cout << s << endl;
//	}
//
//
//	if(cut_sample[0]<0){
//		cout << "Warning : in " << fname << " you must provide Positives cut limits ! Exiting\n";
//		exit(EXIT_FAILURE);
//	}
//
//	for(int ii=1;ii<nb_elem;ii++)
//		if((cut_sample[ii]<0)||((cut_sample[ii]-cut_sample[ii-1])<0)){
//			cout << "Warning : in " << fname << " you must provide a crescent order of Positives cut limits !\nExiting\n";
//			exit(EXIT_FAILURE);
//		}
//
//	file.close();
//}


void copy_ref_pos(fitsfile * fptr, fitsfile *outfptr, string name, long ns_final){

	long ns_temp;
	double *RA,*DEC,*PHI;
	double *RA_bis, *DEC_bis, *PHI_bis;
	int status;

	read_ReferencePosition_from_fits(name, RA, DEC, PHI, ns_temp);


	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];
	PHI_bis = new double [ns_final];
	//	time_bis = new double [ns_final];

	for(long ii = 0; ii< ns_final; ii++){
		RA_bis[ii]=RA[ii]*15.0;
		DEC_bis[ii]=DEC[ii];
		PHI_bis[ii]=PHI[ii];
		//		time_bis[ii]=time[ii];
	}



	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);


	// insert column
	//fits_insert_col(fptr, 1, TDOUBLE, , &status);
	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, RA_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_final, DEC_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_final, PHI_bis, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_final, (char*)"Number of rows", &status);


	delete [] RA;
	delete [] DEC;
	delete [] PHI;
	delete [] RA_bis;
	delete [] DEC_bis;
	delete [] PHI_bis;

}

void copy_offsets(fitsfile * fptr, fitsfile *outfptr){

	int status;

	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);

	for(int col=1;col<4;col++)
		fits_copy_col(fptr, outfptr,  col, col,	0, &status);

}


void copy_channels(fitsfile * fptr, fitsfile *outfptr){

	int status;


	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);


	fits_copy_col(fptr, outfptr,  1, 1,	0, &status);

}

void copy_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final){

	int status;
	double *time_bis;

	time_bis = new double [ns_final];

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);

	for(long ii = 0; ii< ns_final; ii++){
		time_bis[ii]=time[ii];
	}

	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	delete [] time_bis;

}

void copy_signal(fitsfile * fptr, fitsfile *outfptr, string name, long ns_final, struct detectors det){

	int status;
	long ns_temp;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	double *signal, *signal_bis;
	signal_bis = new double [ns_final];
	for(long jj=0;jj<det.ndet;jj++){

		read_signal_from_fits(name, det.boloname[jj], signal, ns_temp);
		for(long ii = 0; ii< ns_final; ii++)
			signal_bis[ii]=signal[ii];

		long fpixel[2]={1,jj+1};
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, signal_bis, &status);
	}

	delete [] signal_bis;
	delete [] signal;

}

void copy_mask(fitsfile * fptr, fitsfile *outfptr,  string name, long ns_final, struct detectors det){

	int status;
	long ns_temp;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);


	int *mask, *mask_bis;
	mask_bis = new int [ns_final];
	for(long jj=0;jj<det.ndet;jj++){

		read_flag_from_fits(name, det.boloname[jj], mask, ns_temp);
		for(long ii = 0; ii< ns_final; ii++)
			mask_bis[ii]=mask[ii];

		long fpixel[2]={1,jj+1};
		fits_write_pix(outfptr, TINT, fpixel, ns_final, mask_bis, &status);
	}

	delete [] mask_bis;
	delete [] mask;
}

