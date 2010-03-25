
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneFix.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}

using namespace std;


int main(int argc, char *argv[]) {

	int parsed=1;
	int rank=0;

	struct samples samples_struct;
	struct common dir;
	struct detectors det;
	std::vector<long> indice;
	std::vector <long> add_sample;
	double fsamp;
	//	double *RA, *DEC, *PHI;

	string outname;

	if(rank==0)
		printf("\nBeginning of saneFix:\n\n");


	parsed=parse_saneFix_ini_file(argv[1],dir,
			samples_struct, rank);


	for(long ii=0; ii<samples_struct.ntotscan;ii++){
		int format_fits;
		fitsfile * fptr;
		fitsfile *outfptr;


		if((read_indices_file(samples_struct.fitsvect[ii],dir,  indice, fsamp))){
			cout << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run saneCheck on this file before\n";
			continue;
		}

		format_fits=test_format(samples_struct.fitsvect[ii]);

		double *time;
		time=new double[samples_struct.nsamples[ii]];
		long samples_to_add=how_many(samples_struct.fitsvect[ii], samples_struct.nsamples[ii] ,indice, time, fsamp, add_sample);
		cout << "samptoadd : " << samples_to_add << endl;
		cout << "total : " << samples_struct.nsamples[ii] + samples_to_add << endl;
		delete [] time;

		long ns_total = samples_struct.nsamples[ii] + samples_to_add;


		std::ostringstream oss;
		oss << samples_struct.fitsvect[ii];
		string filename = oss.str();
		string fname2 = Basename(filename) + "_fixed.fits";
		oss.str("");
		oss << "!" << dir.output_dir << fname2;
		string temp = oss.str();
		oss.str("");
		string fname=samples_struct.fitsvect[ii];

		int status=0;
		struct detectors det;

		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		if (fits_create_file(&outfptr, temp.c_str(), &status))
			fits_report_error(stderr, status);

		// read channels list
		read_bolo_list(samples_struct.fitsvect[ii], det);

		// Copy primary Header
		fits_copy_header(fptr, outfptr, &status);


		if(format_fits==1){ // HIPE format
			cout << "HIPE format found\n";

			// 1 signal
			fix_signal(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			// 2 RA 3 DEC
			//			fix_RA(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);
			fix_RA_DEC(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			//			// 3 DEC
			//			fix_DEC(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			// 4 mask
			fix_mask(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			// 5 time
			fix_time_table(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample, samples_struct.nsamples[ii], fsamp);

			// 6 channels
			copy_channels(fptr, outfptr);

			// 7 reference positions
			fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			// 8 offsets
			copy_offsets(fptr, outfptr);

		}else{
			cout << "SANEPIC format found\n";

			// 1 ref pos
			fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			// 2 offsets
			copy_offsets(fptr, outfptr);

			// 3 channels
			copy_channels(fptr, outfptr);

			// 4 time
			fix_time_table(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample, samples_struct.nsamples[ii], fsamp);

			// 5 signal
			fix_signal(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

			// 6 mask
			fix_mask(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, indice, add_sample);

		}


		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		if (fits_close_file(outfptr, &status))
			fits_report_error(stderr, status);


	}

	//cleaning
	delete [] samples_struct.nsamples;

	if(rank==0)
		cout << "\nEnd of saneFix\n";



}
