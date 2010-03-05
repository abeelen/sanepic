
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

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;


int main(int argc, char *argv[]) {


	int rank, size;

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
	cout << "Mpi is not used for this step" << endl;
#endif

	/*
	 * TODO : This should be organized as :
	 *
	 *  - parse the input ini file (OK)
	 *  - check for existence of directory/files pointed from the ini file (OK)
	 *  - for each file, check for:
	 *      - presence of the different hdus (only the position change between the two format) (OK)
	 *      - consistent sizes of the different hdu (OK)
	 *      - presence of non flagged NaN values (OK)
	 *
	 *      - bad channels [latter, noisier channels] (OK)
	 *
	 *      - check the fsamp from the ini file, should be the most frequent one (OK)
	 *      - check for time gaps > 1.9 * the most frequent time gap (OK)
	 */


	// data parameters
	/*!
	 * -ndet = number of detectors to output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	//	long ndet, ndetOrig, nbins;
	//	int nbolos;
	int parsed=1;

	struct samples samples_struct;
	struct common dir;
	struct detectors det;
	double fsamp;

	string outname;
	//	string temp;
	//	size_t found;

	if(rank==0)
		printf("\nBeginning of saneCheck:\n\n");


	parsed=parse_saneCheck_ini_file(argv[1],dir,
			det, samples_struct, fsamp, rank);

	if(parsed==-1){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}

	readFrames(samples_struct.fitsvect, samples_struct.nsamples);


	struct detectors bolo_fits;
	std::vector<std::string> bolo_bad;
	std::vector<std::string> bolo_bad_80;



	// TODO : sort the processors with mpi

	for(int ii=0;ii<samples_struct.ntotscan;ii++){

		int format_fits=0;

		cout << endl << endl << "[" << rank <<  "] Checking : " << samples_struct.fitsvect[ii] << endl << endl;

		read_bolo_list(samples_struct.fitsvect[ii],bolo_fits);

		format_fits=test_format(samples_struct.fitsvect[ii]); // format = 1 => HIPE, else Sanepic

		check_detector_is_in_fits(det,bolo_fits,samples_struct.fitsvect[ii]);




		cout << "\n[" << rank <<  "] Checking presence of common HDU and position HDU\n";
		check_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
		check_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, format_fits);

		if(format_fits==1){
			cout << "[" << rank <<  "] HIPE format found, Checking Alt position HDU presence\n";
			check_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
		}

		cout << "\n[" << rank <<  "] Checking NANs in common HDU and position HDU\n";
		check_NAN_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
		check_NAN_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);

		if(format_fits==1){
			cout << "[" << rank <<  "] HIPE format found, Checking NANs in Alt position HDU\n";
			check_NAN_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
		}

		cout << "\n[" << rank <<  "] Checking time gaps in time table\n";
		check_time_gaps(samples_struct.fitsvect[ii],samples_struct.nsamples[ii], fsamp, dir);



		//		temp = samples_struct.fitsvect[ii];
		//		found=temp.find_last_of('/');
		//		outname = dir.output_dir + temp.substr(found+1) + "_bolos_flag.txt";

		check_flag(samples_struct.fitsvect[ii],bolo_fits, samples_struct.nsamples[ii],outname, bolo_bad,bolo_bad_80);


	}

	//	outname = dir.tmp_dir +
	//	indice_files_generation(outname,indice);

	if(rank==0){
		cout << endl;

		// generating log files :
		if(bolo_bad.size()>0){
			outname = dir.output_dir + "bolo_totally_flagged.txt";
			cout << "Writing informations in :\n" << outname << endl << endl;
			log_gen(bolo_bad,outname);
		}else
			cout << "There are no bolometers fully flagged in the files\n\n";


		if(bolo_bad_80.size()>0){
			outname = dir.output_dir + "bolo_80_percent_flagged.txt";
			cout << "Writing informations in :\n" << outname << endl;
			log_gen(bolo_bad_80,outname);
		}else
			cout << "There are no bolometers more than 80% flagged in the files\n";


		cout << "\nEnd of saneCheck\n";
	}
	//cleaning
	delete [] samples_struct.nsamples;


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

}
