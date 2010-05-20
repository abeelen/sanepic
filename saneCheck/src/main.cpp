
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


	int rank, size; /* MPI processor rank and MPI total number of used processors */

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
	 *  This is organized as :
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
	 *      - generate bad detectors list and temporary files for saneFix
	 */


	int parsed=0; /* parser error status */

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct common dir;  /*! structure that contains output input temp directories */
	struct detectors det;  /*! A structure that contains everything about the detectors names and number */
	double fsamp; /*! sampling frequency */

	string outname; /*! Ouput log files name */

	if(rank==0)
		printf("\nBeginning of saneCheck:\n\n");

	if (argc<2) /* not enough argument */
		parsed=-1;
	else {
		parsed=parse_saneCheck_ini_file(argv[1],dir,
				det, samples_struct, fsamp, rank);
	}

	if(parsed==-1){ /* error during parsing phase */
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}

	struct detectors bolo_fits_0; /* bolometers list of the first fits file given as input */
	read_bolo_list(samples_struct.fitsvect[0],bolo_fits_0); /* Read the first fits file bolo table */

	long *bolo_bad_tot; /*! bad detectors full list => fully flag detectors */
	long *bolo_bad_80_tot; /*! valid worst detectors full list => more than 80% flag detectors */

	if(rank==0){ // only for MPI_reduce
		bolo_bad_tot= new long [bolo_fits_0.ndet];
		bolo_bad_80_tot= new long [bolo_fits_0.ndet];
		fill(bolo_bad_tot, bolo_bad_tot + bolo_fits_0.ndet ,0);
		fill(bolo_bad_80_tot, bolo_bad_80_tot + bolo_fits_0.ndet ,0);
	}


	for(int ii=0;ii<samples_struct.ntotscan;ii++){ /* for each input fits file */

		int do_it=who_do_it(size, rank, ii); /* which rank do the job ? */

		if(rank==do_it){ /* if this rank has to do the job ... */

			struct detectors bolo_fits;
			long *bolo_bad;
			long *bolo_bad_80;

			struct checkHDU Check_it;

			// initialize default values
			Check_it.checkDEC=1;
			Check_it.checkRA=1;
			Check_it.checkREFERENCEPOSITION=1;
			Check_it.checkOFFSETS=1;
			Check_it.checkSIGNAL=1;
			Check_it.checkFLAG=1;
			Check_it.checkTIME=1;

			int format_fits=0; /* fits format indicator : 1 HIPE, 2 Sanepic */

			cout << endl << endl << "[" << rank <<  "] Checking : " << samples_struct.fitsvect[ii] << endl << endl;

			format_fits=test_format(samples_struct.fitsvect[ii]); // format = 1 => HIPE, else Sanepic

			read_bolo_list(samples_struct.fitsvect[ii],bolo_fits); // read fits file detector list
			if(check_bolos(bolo_fits.boloname, bolo_fits_0.boloname)){ // compare to the first input fits file to ensure compatibility between scans
				cout << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run only together scans that correspond to the same field\n";
				continue;
			}

			check_detector_is_in_fits(det,bolo_fits,samples_struct.fitsvect[ii]); // check wether used detector user list is correct

			bolo_bad = new long[bolo_fits.ndet]; // this scan bad detectors list
			bolo_bad_80 = new long[bolo_fits.ndet]; // this scan valid worst detectors list
			fill(bolo_bad,bolo_bad+bolo_fits.ndet,0);
			fill(bolo_bad_80,bolo_bad_80+bolo_fits.ndet,0);

			cout << "\n[" << rank <<  "] Checking presence of common HDU and position HDU\n";
			check_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits,Check_it); // check presence of channels, time, signal and mask HDUs
			check_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, format_fits,Check_it); // check presence of reference positions and offsets HDUs

			if(format_fits==1){ // check RA/DEC table presence for HIPE format
				cout << "[" << rank <<  "] HIPE format found, Checking Alt position HDU presence\n";
				check_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits,Check_it);
			}

			cout << "\n[" << rank <<  "] Checking NANs in common HDU and position HDU\n"; // check non-flag NANs presence in whole tables
			check_NAN_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits,Check_it);
			check_NAN_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits,Check_it);

			if(format_fits==1){ // check NANs presence in RA/DEc tables for HIPE format
				cout << "[" << rank <<  "] HIPE format found, Checking NANs in Alt position HDU\n";
				check_NAN_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits,Check_it);
			}

			cout << "\n[" << rank <<  "] Checking time gaps in time table\n"; // check for time gaps in time table
			check_time_gaps(samples_struct.fitsvect[ii],samples_struct.nsamples[ii], fsamp, dir,Check_it);

			// Lookfor fully or more than 80% flagged detectors, also flag singletons
			check_flag(samples_struct.fitsvect[ii],bolo_fits, samples_struct.nsamples[ii],outname, bolo_bad,bolo_bad_80,Check_it);

			if(Check_it.checkFLAG){
#ifdef USE_MPI
				// inform processor 0 of bad or worst bolometer presence in fits files
				MPI_Reduce(bolo_bad,bolo_bad_tot,bolo_fits_0.ndet,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
				MPI_Reduce(bolo_bad_80,bolo_bad_80_tot,bolo_fits_0.ndet,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#else
				for(long kk = 0; kk< bolo_fits_0.ndet; kk++){ // sum up the bad bolometers status
					bolo_bad_tot[kk]+=bolo_bad[kk];
					bolo_bad_80_tot[kk]+=bolo_bad_80[kk];
				}
#endif
			}
			delete [] bolo_bad;
			delete [] bolo_bad_80;
		}
	}

	if(rank==0){
		cout << endl;

		// generating log files :
		outname = dir.output_dir + "bolo_totally_flagged.txt";
		cout << "Writing informations in :\n" << outname << endl << endl;
		log_gen(bolo_bad_tot,outname, bolo_fits_0); // generate bad detectors log file


		outname = dir.output_dir + "bolo_80_percent_flagged.txt";
		cout << "Writing informations in :\n" << outname << endl;
		log_gen(bolo_bad_80_tot, outname, bolo_fits_0); // generate valid worst detectors log file


		cout << "\nEnd of saneCheck\n";
	}

	//clean up
	delete [] samples_struct.nsamples;
	delete [] bolo_bad_tot;
	delete [] bolo_bad_80_tot;

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	return EXIT_SUCCESS;
}
