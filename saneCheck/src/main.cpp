
#include "covMatrixIO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneCheck.h"
#include "tools.h"

#include <iostream>
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

	// data parameters
	/*!
	 * -ndet = number of detectors to output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	//	long ndet, ndetOrig, nbins;
	//	int nbolos;
	int rank = 0;
	int parsed=1;
	int checked = 0;

	struct samples samples_struct;
	struct directories dir;
	struct detectors det;

	string outname;
	string temp;
	size_t found;

	printf("\nBeginning of saneCheck:\n\n");


	parsed=parse_saneCheck_ini_file(argv[1],dir,
			det, samples_struct, rank);

	//	outname = dir.outdir;

	//cout << samples_struct.fitsvect.size() << endl;

	readFrames(samples_struct.fitsvect, samples_struct.nsamples);

	//	cout << "ntotscan : " << samples_struct.ntotscan;

	std::vector<std::string> bolo_bad;
	std::vector<std::string> bolo_bad_80;

	for(int ii=0;ii<samples_struct.ntotscan;ii++){

		cout << endl << endl << "Checking : " << samples_struct.fitsvect[ii] << endl;

		read_bolo_list(samples_struct.fitsvect[ii],det);

		check_hdu(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],det);
		check_NaN(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],det);
		check_time_gaps(samples_struct.fitsvect[ii],samples_struct.nsamples[ii]);

		temp = samples_struct.fitsvect[ii];
		found=temp.find_last_of('/');
		outname = dir.outdir + temp.substr(found+1) + "_bolos_flag.txt";
		cout << "Writing informations in :\n" << outname << endl << endl;

		checked=check_flag(samples_struct.fitsvect[ii],det, samples_struct.nsamples[ii],outname, bolo_bad,bolo_bad_80);
	}



	// generating log files :
	outname = dir.outdir + "bolo_totally_flagged.txt";
	log_gen(bolo_bad,outname);

	outname = dir.outdir + "bolo_80_percent_flagged.txt";
	log_gen(bolo_bad_80,outname);


	//cleaning
	delete [] samples_struct.nsamples;

	cout << "\nEnd of saneCheck\n";


}
