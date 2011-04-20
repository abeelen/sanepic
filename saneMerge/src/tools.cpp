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
#include <fitsio.h>
}



using namespace std;

void file_compatibility_verification(struct samples samples_struct)
/*! this function tests that the detectors lists are the same in whole fits files, also tests that the files have a crescent time reference */
{

	long numscan = samples_struct.ntotscan;
	fitsfile * fptr1; // 2 fits file pointers
	fitsfile * fptr2;
	int status = 0; // fits error status

	for(long ii = 1; ii< numscan; ii++ ){ // for each scan

		double *time_one, *time_two; // 2 fits input time tables
		std::vector<string> det_one;
		std::vector<string> det_two; // 2 detectors lists
		long ndet1, ndet2;

		// TODO : Do we have to test that offsets are the same ??

		string file1, file2; // 2 filenames

		file1 = samples_struct.fitsvect[ii-1]; // file2 will be checked in regards to file1
		file2 = samples_struct.fitsvect[ii];

		// open both files
		if (fits_open_file(&fptr1, file1.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		if (fits_open_file(&fptr2, file2.c_str(), READONLY, &status))
			fits_report_error(stderr, status);


		long ns1, ns2; // get each number of samples
		ns1=samples_struct.nsamples[ii-1];
		ns2=samples_struct.nsamples[ii];

		read_time_from_fits(file1, time_one, ns1); // read both time tables
		read_time_from_fits(file2, time_two, ns2);

		if(time_one[ns1-1]>=time_two[0]){ // check time is crescent between 1st and 2nd file
			cout << "There is a problem between " << file1 << " and " << file2 << " time tables.\n";
			cout << "Please order them correctly in order to have a time ordered list of input files. Exiting...\n";
			exit(EXIT_FAILURE);
		}

		read_bolo_list(file1, det_one, ndet1); // read both bolometers lists and store in detectors structure
		read_bolo_list(file2, det_two, ndet2);

		if(ndet1!=ndet2){ // check number of detector is the same
			cout << "Error ! The channel lists in " << file1 << " and " << file2 << " are not the same. Exiting... \n";
			exit(EXIT_FAILURE);
		}

		for(long jj=0; jj<ndet2; jj++) // check each detector from first file is present in second one
			if(det_one[jj]!=det_two[jj]){
				cout << "Error ! The channel lists in " << file1 << " and " << file2 << " are not the same.\n";
				cout << det_one[jj] << " != " << det_two[jj] << ". Exiting...\n";
				exit(EXIT_FAILURE);
			}


		// ckean up
		delete [] time_one;
		delete [] time_two;
		det_one.clear();
		det_two.clear();
		ndet1=0;
		ndet2=0;

		// close both fits file
		if (fits_close_file(fptr1, &status))
			fits_report_error(stderr, status);
		if (fits_close_file(fptr2, &status))
			fits_report_error(stderr, status);


	} // ii loop : 1 -> numscan

}

void copy_ref_pos(fitsfile *outfptr, struct samples samples_struct, long ns_final)
/*! copy reference position tables from each file to output file */
{


	fitsfile * fptr; // fits file pointer
	long indice_debut=0; // output last written sample indice
	double *RA_bis, *DEC_bis, *PHI_bis; // RA DEC and PHI output tables

	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];
	PHI_bis = new double [ns_final];

	int status=0; // fits error status

	for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){ // for each scan

		long ns_temp=0;


		string fname=samples_struct.fitsvect[iframe];
		//		indice_debut = indice_fin;


		// open each fits file and copy each table in the final file
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		if(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status)) // move HDU pointer to desired table
			fits_report_error(stderr, status);

		if(iframe==0){ // if this is the first file of the list
			if(fits_copy_header(fptr, outfptr, &status))
				fits_report_error(stderr, status);
			if(fits_movnam_hdu(outfptr, BINARY_TBL, (char*) "reference position", NULL, &status)) // move HDU pointer to desired table
				fits_report_error(stderr, status);
			if(fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_final, (char*)"Number of rows", &status)) // update table key
				fits_report_error(stderr, status);
		}

		double *RA, *DEC, *PHI; // input RA DEC PHI tables

		read_ReferencePosition_from_fits(fname, RA, DEC, PHI, ns_temp); // read RA DEC and PHI table from input file

		for(long ii = 0; ii< ns_temp; ii++){ // add RA DEC and PHI values to the ouput table
			RA_bis[indice_debut+ii]=RA[ii]*15.0;
			DEC_bis[indice_debut+ii]=DEC[ii];
			PHI_bis[indice_debut+ii]=PHI[ii];
		}

		// update last written sample indice
		indice_debut+=samples_struct.nsamples[iframe];

		if (fits_close_file(fptr, &status)) // close input fits file
			fits_report_error(stderr, status);

		//clean up input tables
		delete [] RA;
		delete [] DEC;
		delete [] PHI;

	} // end of iframe loop : 1 -> ntotscan

	// write output tables RA DEC and PHI in output file
	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, RA_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_final, DEC_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_final, PHI_bis, &status);

	// clean up
	delete [] RA_bis;
	delete [] DEC_bis;
	delete [] PHI_bis;

}

void copy_time(fitsfile *outfptr, struct samples samples_struct, long ns_final)
/*! copy time tables from each file to output file */
{


	fitsfile * fptr; // fits file pointer
	long indice_debut=0; // output last written sample indice
	double *time_bis; // output file table

	time_bis = new double [ns_final];

	int status=0; // fits error status

	for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){ // for each scan

		long ns_temp=samples_struct.nsamples[iframe]; // get number of samples

		string fname=samples_struct.fitsvect[iframe]; // get filename

		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open input file
			fits_report_error(stderr, status);


		fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status); // move HDU pointer to time table

		if(iframe==0){ // if first scan ...
			fits_copy_header(fptr, outfptr, &status); // copy talb header
			fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status); // update number of rows
		}

		double *time;
		read_time_from_fits(fname, time, ns_temp); // read time from input file

		for(long ii = 0; ii< ns_temp; ii++){ // add in the output time table
			time_bis[indice_debut + ii]=time[ii];
		}

		indice_debut += samples_struct.nsamples[iframe]; // update last written sample pointer

		// clean up
		delete [] time;

		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

	}
	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status); // write time table in output file


	// clean up
	delete [] time_bis;

}

void copy_signal(fitsfile *outfptr, struct samples samples_struct, std::vector<std::string> det, long ndet, long ns_final)
/*! copy signal tables from each file to output file */
{

	fitsfile * fptr; // fits pointer
	long indice_debut=0; // output last written sample indice

	double  *signal_bis; // output signal table
	signal_bis = new double [ns_final];

	for(long jj=0;jj<ndet;jj++){ // for each detector
		int status=0; // fits error status

		for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){ // for each scan

			long ns_temp=0; // input number of samples


			string fname=samples_struct.fitsvect[iframe]; // input filename

			// open fits file
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);

			// move HDU pointer to signal table
			fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status);

			if((jj==0)&&(iframe==0)){ // copy header only once !
				fits_copy_header(fptr, outfptr, &status); // copy table header
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status); // and update it
			}

			double *signal;

			read_signal_from_fits(fname, det[jj], signal, ns_temp); // read input signal table
			for(long ii = 0; ii< ns_temp; ii++)
				signal_bis[indice_debut+ii]=signal[ii]; // add to signal output table


			indice_debut += samples_struct.nsamples[iframe]; // update last written sample pointer value

			delete [] signal;

			// close file
			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

		} // end of iframe loop

		string field= det[jj]; // actual detector name

		string fname=samples_struct.fitsvect[0]; // first input file name

		// open first input file
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		// get this detector index in first file channel list
		long rowIndex = find_channel_index(fptr, field.c_str());

		// close first file
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		// store this index
		long fpixel[2]={1,rowIndex};

		// write the row in the ouput table
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, signal_bis, &status);
		indice_debut=0;
	} // end of jj loop : 0 -> ndet

	// clean up
	delete [] signal_bis;

}

void copy_mask(fitsfile *outfptr, struct samples samples_struct, std::vector<std::string> det, long ndet, long ns_final)
/*! copy flag tables from each file to output file */
{

	fitsfile * fptr; // fits file pointer
	long indice_debut=0; // last written sample pointer

	double  *mask_bis; // output flag table
	mask_bis = new double [ns_final];

	for(long jj=0;jj<ndet;jj++){ // for each detector in detector list
		int status=0; // fits error status

		for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){ // for each scan

			long ns_temp=0; // input number of samples
			string fname=samples_struct.fitsvect[iframe]; // input file name


			// open fits file a
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);

			// move HDU pointer to mask
			fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status);

			if((jj==0)&&(iframe==0)){ // copy header only once
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);
			}

			double *mask;

			read_signal_from_fits(fname, det[jj], mask, ns_temp); // read input file flag table

			for(long ii = 0; ii< ns_temp; ii++) // add flag values to ouput flag table
				mask_bis[indice_debut+ii]=mask[ii];


			indice_debut += samples_struct.nsamples[iframe]; // update last written sample index

			delete [] mask;

			// close input fits file
			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

		}

		string field= det[jj]; // actual bolometer's name
		string fname=samples_struct.fitsvect[0]; // first input file name

		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open first input file
			fits_report_error(stderr, status);

		long rowIndex = find_channel_index(fptr, field.c_str()); // find actual bolometer index in detector input list

		// close first file
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		// store bolometer index
		long fpixel[2]={1,rowIndex};
		// write row in output file
		fits_write_pix(outfptr, TINT, fpixel, ns_final, mask_bis, &status);

		// reset written sample for next bolometer
		indice_debut=0;

	} // end of bolo loop

	// clean up
	delete [] mask_bis;
}


void copy_RA_DEC(fitsfile *outfptr, struct samples samples_struct, std::vector<std::string> det, long ndet, long ns_final)
/*! copy RA and DEC tables (HIPE format only) from each file to output file */
{

	fitsfile * fptr; // fits file pointer
	long indice_debut=0;

	double *RA_bis; // output RA and DEC tables
	double *DEC_bis;
	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];

	for(long jj=0;jj<ndet;jj++){ // for each detector
		int status=0;

		for(long iframe=0; iframe<samples_struct.ntotscan;iframe++){ // for each scan

			long ns_temp=0; // actual fits file number of sample


			string fname=samples_struct.fitsvect[iframe]; // actual input file name


			// open fits file
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);


			if((jj==0)&&(iframe==0)){ // copy headers once
				fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status);
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

				fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status);
				fits_copy_header(fptr, outfptr, &status);
				fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);
			}

			double *RA, *DEC;
			read_ra_dec_from_fits(fname, det[jj], RA, DEC, ns_temp); // read RA and DEC tables

			for(long ii = 0; ii< ns_temp; ii++){ // add RA and DEC values to output tables
				RA_bis[indice_debut + ii]=RA[ii];
				DEC_bis[indice_debut + ii]=DEC[ii];
			}



			indice_debut += samples_struct.nsamples[iframe]; // update last written index

			if (fits_close_file(fptr, &status)) //close file
				fits_report_error(stderr, status);

			//clean up
			delete [] RA;
			delete [] DEC;
		} // end of scan loop

		// for this bolometer, find index in first input detector list
		string field= det[jj];
		string fname=samples_struct.fitsvect[0];

		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open fits file
			fits_report_error(stderr, status);

		long rowIndex = find_channel_index(fptr, field.c_str()); // find detector index

		// store it
		long fpixel[2]={1,rowIndex};

		if (fits_close_file(fptr, &status)) // close input file
			fits_report_error(stderr, status);

		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "ra", NULL, &status); // move HDU pointer to RA table to write down the row
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, RA_bis, &status); // write the row
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "dec", NULL, &status); // move HDU pointer to DEC table to write down the row
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, DEC_bis, &status); // write the row
		indice_debut=0;
	} // end of bolo loop

	//clean up
	delete [] RA_bis;
	delete [] DEC_bis;


}






