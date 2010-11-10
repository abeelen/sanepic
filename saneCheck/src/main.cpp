#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneCheck.h"
#include "parser_functions.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>


extern "C" {
#include "nrutil.h"
}

#ifdef PARA_FRAME
#include "mpi.h"
#endif

using namespace std;


int main(int argc, char *argv[]) {


	int rank, size; /* MPI processor rank and MPI total number of used processors */

#ifdef PARA_FRAME

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
	struct param_common dir;  /*! structure that contains output input temp directories */
	//	struct detectors det;  /*! A structure that contains everything about the detectors names and number */
	//	std::vector<detectors> detector_tab;
	//	std::vector<double> bolometer_gain;
	samples_struct.fits_table=NULL;

	struct param_sanePos pos_param;
	struct param_sanePre proc_param;
	struct param_sanePic sanePic_struct;
	struct param_saneInv saneInv_struct;
	struct param_sanePS structPS;
	std::vector<double> fcut;
	struct saneCheck check_struct;
	string outname; /*! Ouput log files name */
	string output = "";
	string bolo_gain_filename=""; // TODO change this

	if(rank==0)
		printf("\nBeginning of saneCheck:\n\n");

	if (argc<2) /* not enough argument */
		parsed=-1;
	else {
		parsed=parse_saneCheck_ini_file(argv[1], output, dir, samples_struct, pos_param, proc_param, fcut,
				structPS, saneInv_struct, sanePic_struct, check_struct, rank, size);

		if(rank==0)
			// print parser warning and/or errors
			cout << endl << output << endl;
	}

	if(parsed==-1){ /* error during parsing phase */
		if(rank==0)
			cout << "Error during parsing step. Exiting ...\n";
#ifdef PARA_FRAME
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}


	if(read_bolo_for_all_scans(dir, samples_struct, rank, size)){
#ifdef PARA_FRAME
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return(EX_IOERR);
	}

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, sanePic_struct);
		print_saneCheck_ini(check_struct);
	}


	std::vector<string> bolo_fits_0; /* bolometers list of the first fits file given as input */
	long ndet0;
	read_bolo_list(samples_struct.fitsvect[0],bolo_fits_0,ndet0); /* Read the first fits file bolo table */

	long *bolo_bad_tot = NULL; /*! bad detectors full list => fully flag detectors */
	long *bolo_bad_80_tot = NULL; /*! valid worst detectors full list => more than 80% flag detectors */

	if(rank==0){ // only for MPI_reduce
		bolo_bad_tot= new long [ndet0];
		bolo_bad_80_tot= new long [ndet0];
		fill(bolo_bad_tot, bolo_bad_tot + ndet0 ,0);
		fill(bolo_bad_80_tot, bolo_bad_80_tot + ndet0 ,0);
	}


	for(int ii=0;ii<samples_struct.ntotscan;ii++){ /* for each input fits file */

		int do_it=who_do_it(size, rank, ii); /* which rank do the job ? */

		if(rank==do_it){ /* if this rank has to do the job ... */

			std::vector<string> bolo_fits;
			long ndet_fits;
			long *bolo_bad;
			long *bolo_bad_80;
			double *percent_tab;
			int return_value=0;

			// initialize default values
			check_struct.Check_it.checkDEC=1;
			check_struct.Check_it.checkRA=1;
			check_struct.Check_it.checkREFERENCEPOSITION=1;
			check_struct.Check_it.checkOFFSETS=1;

			int format_fits=0; /* fits format indicator : 1 HIPE, 2 Sanepic */

			cout << endl << endl << "[" << rank <<  "] Checking : " << samples_struct.fitsvect[ii] << endl << endl;

			format_fits=test_format(samples_struct.fitsvect[ii]); // format = 1 => HIPE, else Sanepic
			if(format_fits==0){
				cerr << "input fits file format is undefined : " << samples_struct.fitsvect[ii] << " . Exiting...\n";
#ifdef PARA_FRAME
				MPI_Finalize();
#endif
				return EX_IOERR;
			}

			read_bolo_list(samples_struct.fitsvect[ii],bolo_fits,ndet_fits); // read fits file detector list
			if(check_bolos(bolo_fits, bolo_fits_0)){ // compare to the first input fits file to ensure compatibility between scans
				cout << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run only together scans that correspond to the same field\n";
				continue;
			}

			string output_read = "";
			std::vector<string> det_vect;
			if(read_channel_list(output_read, samples_struct.bolovect[ii], det_vect)){
				cout << output_read << endl;
#ifdef PARA_FRAME
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return EX_CONFIG;
			}
			long ndet_vect = (long)det_vect.size();

			check_detector_is_in_fits(det_vect, ndet_vect, bolo_fits,samples_struct.fitsvect[ii]); // check wether used detector user list is correct

			bolo_bad = new long[ndet_fits]; // this scan bad detectors list
			bolo_bad_80 = new long[ndet_fits]; // this scan valid worst detectors list
			percent_tab = new double[ndet_fits];
			fill(bolo_bad,bolo_bad+ndet_fits,0);
			fill(bolo_bad_80,bolo_bad_80+ndet_fits,0);
			fill(percent_tab,percent_tab+ndet_fits,0.0);

			cout << "\n[" << rank <<  "] Checking presence of common HDU and position HDU\n";
			return_value+=check_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],ndet_fits,check_struct.Check_it); // check presence of channels, time, signal and mask HDUs
			return_value+=check_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],ndet_fits, format_fits,check_struct.Check_it); // check presence of reference positions and offsets HDUs
			if(format_fits==1){ // check RA/DEC table presence for HIPE format
				cout << "[" << rank <<  "] HIPE format found, Checking Alt position HDU presence\n";
				return_value+=check_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],ndet_fits,check_struct.Check_it);
			}


			if(return_value<0){
#ifdef PARA_FRAME
				MPI_Finalize();
#endif
				return EX_SOFTWARE;
			}



			if(((format_fits==2)&&((!check_struct.Check_it.checkREFERENCEPOSITION)||(!check_struct.Check_it.checkOFFSETS))) ||
					((format_fits==1)&&((!check_struct.Check_it.checkRA)||(!check_struct.Check_it.checkDEC)))){
				cout << "NO POSITION TABLES ARE PRESENTS : EXITING ...\n";
#ifdef PARA_FRAME
				MPI_Finalize();
#endif
				return EX_IOERR;
			}

			if(check_struct.checkNAN){
				cout << "\n[" << rank <<  "] Checking NANs in common HDU and position HDU\n"; // check non-flag NANs presence in whole tables
				check_NAN_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, ndet_fits,check_struct.Check_it);
				check_NAN_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, ndet_fits,check_struct.Check_it);
			}

			if((format_fits==1)&&(check_struct.checkNAN)){ // check NANs presence in RA/DEc tables for HIPE format
				cout << "[" << rank <<  "] HIPE format found, Checking NANs in Alt position HDU\n";
				check_NAN_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, ndet_fits,check_struct.Check_it);
			}

			if(check_struct.checktime){
				cout << "\n[" << rank <<  "] Checking time gaps in time table\n"; // check for time gaps in time table
				check_time_gaps(samples_struct.fitsvect[ii],samples_struct.nsamples[ii], proc_param.fsamp, dir,check_struct.Check_it);
			}

			if(check_struct.checkGain){
				//				cout << "\n[" << rank <<  "] Computing and Checking bolometer gain correction in signal table\n"; // check for time gaps in time table
				//				check_bolo_gain(samples_struct.fitsvect[ii],samples_struct.nsamples[ii], bolo_gain_filename, det, check_struct.Check_it); // TODO : quelle bolo list on prend ?? selon les differents cas !!
				//				getchar();
			}

			if(check_struct.checkflag){
				// Lookfor fully or more than 80% flagged detectors, also flag singletons
				cout << "\n[" << rank <<  "] Checking flagged detectors\n"; // check for time gaps in time table
				check_flag(samples_struct.fitsvect[ii],bolo_fits,ndet_fits, samples_struct.nsamples[ii],outname, bolo_bad,bolo_bad_80,percent_tab,check_struct.Check_it);

				// generating log files :
				outname = dir.output_dir + "bolo_totally_flagged_" + FitsBasename(samples_struct.fitsvect[ii]) +".txt";
				cout << "Writing informations in :\n" << outname << endl << endl;
				log_gen(bolo_bad,outname, bolo_fits_0,ndet0); // generate bad detectors log file


				outname = dir.output_dir + "bolo_80_percent_flagged_" + FitsBasename(samples_struct.fitsvect[ii]) +".txt";
				cout << "Writing informations in :\n" << outname << endl;
				log_gen(bolo_bad_80, outname, bolo_fits_0, ndet0, percent_tab); // generate valid worst detectors log file
			}

#ifdef PARA_FRAME
			// inform processor 0 of bad or worst bolometer presence in fits files
			MPI_Reduce(bolo_bad,bolo_bad_tot,ndet0,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(bolo_bad_80,bolo_bad_80_tot,ndet0,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#else
			for(long kk = 0; kk< ndet0; kk++){ // sum up the bad bolometers status
				bolo_bad_tot[kk]+=bolo_bad[kk];
				bolo_bad_80_tot[kk]+=bolo_bad_80[kk];
			}
#endif
			delete [] bolo_bad;
			delete [] bolo_bad_80;
		}
	}

	if((rank==0)&&(check_struct.checkflag)){
		cout << endl;

		// generating log files :
		outname = dir.output_dir + "bolo_totally_flagged.txt";
		cout << "Writing informations in :\n" << outname << endl << endl;
		log_gen(bolo_bad_tot,outname, bolo_fits_0, ndet0); // generate bad detectors log file


		outname = dir.output_dir + "bolo_80_percent_flagged.txt";
		cout << "Writing informations in :\n" << outname << endl;
		log_gen(bolo_bad_80_tot, outname, bolo_fits_0, ndet0); // generate valid worst detectors log file


		cout << "\nEnd of saneCheck\n";
	}

	//clean up
	delete [] samples_struct.nsamples;
	delete [] bolo_bad_tot;
	delete [] bolo_bad_80_tot;

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	return EXIT_SUCCESS;
}
