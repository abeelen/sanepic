
#include "covMatrixIO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneCheck.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <algorithm>
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


	/*
	 * TODO : This should be organized as :
	 *
	 *  - parse the input ini file
	 *  - check for existence of directory/files pointed from the ini file
	 *  - for each file, check for:
	 *      - presence of the different hdus (only the position change between the two format)
	 *      - consistent sizes of the different hdu
	 *      - presence of non flagged NaN values
	 *
	 *      - bad channels [latter, noisier channels]
	 *
	 *      - check the fsamp from the ini file, should be the most frequent one
	 *      - check for time gaps > 1.9 * the most frequent time gap
	 */


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

	struct samples samples_struct;
	struct common dir;
	struct detectors det;

	string outname;
	string temp;
	size_t found;

	printf("\nBeginning of saneCheck:\n\n");


	parsed=parse_saneCheck_ini_file(argv[1],dir,
			det, samples_struct, rank);

	//	outname = dir.output_dir;

	//cout << samples_struct.fitsvect.size() << endl;

	readFrames(samples_struct.fitsvect, samples_struct.nsamples);

	//	cout << "ntotscan : " << samples_struct.ntotscan;

	struct detectors bolo_fits;
	std::vector<std::string> bolo_bad;
	std::vector<std::string> bolo_bad_80;
	//	std::vector<string>::iterator it;

	int mycount;

	//TODO: MPI this.....

	for(int ii=0;ii<samples_struct.ntotscan;ii++){

		cout << endl << endl << "Checking : " << samples_struct.fitsvect[ii] << endl;

		read_bolo_list(samples_struct.fitsvect[ii],bolo_fits);
		//		it=bolo_fits.boloname.begin();

		//		cout << bolo_fits.ndet << endl;
		//		for(int uu =0;uu<bolo_fits.ndet; uu++)
		//			cout << bolo_fits.boloname[uu] << " ";
		//

		// TODO: Make a function for this...

		for(int jj=0;jj< det.ndet; jj++){
			mycount = (int) count (bolo_fits.boloname.begin(), bolo_fits.boloname.end(), det.boloname[jj]);
			if(mycount==0)
				cout << "Warning ! The detector " << det.boloname[jj] << " is not referenced in the fits " << samples_struct.fitsvect[ii] << endl << endl;
		}

		// TODO: Make 3 function check_commonHDU, check_positionHDU, check_positionAltHDU (hipe format position)
		// TODO: Same for NaN
		check_hdu(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
		check_NaN(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);

		check_time_gaps(samples_struct.fitsvect[ii],samples_struct.nsamples[ii]);


		temp = samples_struct.fitsvect[ii];
		found=temp.find_last_of('/');
		outname = dir.output_dir + temp.substr(found+1) + "_bolos_flag.txt";

		check_flag(samples_struct.fitsvect[ii],bolo_fits, samples_struct.nsamples[ii],outname, bolo_bad,bolo_bad_80);

		cout << "Writing informations in :\n" << outname << endl << endl;



	}



	// generating log files :
	outname = dir.output_dir + "bolo_totally_flagged.txt";
	log_gen(bolo_bad,outname);

	outname = dir.output_dir + "bolo_80_percent_flagged.txt";
	log_gen(bolo_bad_80,outname);


	//cleaning
	delete [] samples_struct.nsamples;

	cout << "\nEnd of saneCheck\n";


}
