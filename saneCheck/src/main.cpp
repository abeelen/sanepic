#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>

#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneCheck.h"
#include "parser_functions.h"
#include "tools.h"
#include "struct_definition.h"


extern "C" {
#include "nrutil.h"
}

#ifdef PARA_FRAME
#include "mpi.h"
#endif

using namespace std;

/*!
 *  This is organized as :
 *
 *  - parse the input ini file
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - for each file, check for:
 *      - presence of the different hdus (only the position change between the two format)
 *      - consistent sizes of the different hdu
 *
 *      - presence of non flagged NaN values
 *      - bad channels [latter, noisier channels]
 *      - check the fsamp from the ini file, should be the most frequent one
 *      - check for time gaps > 1.9 * the most frequent time gap
 *
 *      - generate bad detectors list and temporary files for saneFix
 */

int main(int argc, char *argv[]) {


	int rank, size; /* MPI processor rank and MPI total number of used processors */


	printf("\nBeginning of saneCheck:\n\n");

#ifdef PARA_FRAME

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
#endif

	//	int parsed=0; /* parser error status */

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_common dir;  // structure that contains output input temp directories
	//	std::vector<double> bolometer_gain;

	struct param_sanePos pos_param;
	struct param_sanePre proc_param;
	struct param_sanePic sanePic_struct;
	struct param_saneInv saneInv_struct;
	struct param_sanePS structPS;
	struct saneCheck check_struct;
	string outname; // Ouput log files name
	string output = "";
	string bolo_gain_filename="";

	std::vector<std::vector<std::string> > bolo_list; // this vector contains all bolonames for all the scans


//	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

		uint16_t parsed=0x0000; // parser error status

		uint16_t mask_saneCheck = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
				BOLOFILE_NOT_FOUND | FSAMP_WRONG_VALUE | FITS_FILELIST_NOT_FOUND; // 0x411f

		uint16_t compare_to_mask; // parser error status


		if (argc<2)/* not enough argument */
			compare_to_mask=0x001;
		else {
			/* parse ini file and fill structures */
			parsed=parse_saneCheck_ini_file(argv[1], output, dir, samples_struct, pos_param, proc_param,
					structPS, saneInv_struct, sanePic_struct, check_struct, rank, size);

			compare_to_mask = parsed & mask_saneCheck;
		}

		// print parser warning and/or errors
		cout << endl << output << endl;

		// in case there is a parsing error or the dirfile format file was not created correctly
		if(compare_to_mask>0x0000){


			switch (compare_to_mask){/* error during parsing phase */

			case 0x0001: printf("Please run %s using a correct *.ini file\n",argv[0]);
			break;

			default : printf("Wrong program options or argument. Exiting !\n");
			break;


			}

#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return EX_CONFIG;
		}
//	}

#ifdef PARA_FRAME
	MPI_Bcast(&check_struct.checkNAN, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&check_struct.checktime, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&check_struct.checkGain, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&check_struct.checkflag, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

#endif

	/* ------------------------------------- READ bolo list ----------------------------*/

	if(channel_list_to_vect_list(samples_struct, bolo_list, rank)){
		cout << "error in channel_list_to_vect_list" << endl;
		return EX_CONFIG;
	}

	/* ------------------------------------------------------------------------------------*/


	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, sanePic_struct, saneInv_struct);
		print_saneCheck_ini(check_struct);
	}


	std::vector<string> bolo_fits_0; // bolometers list of the first fits file given as input
	long ndet0;
	read_bolo_list(dir.data_dir + samples_struct.fitsvect[0],bolo_fits_0,ndet0); // Read the first fits file bolo table

	long *bolo_bad_tot = NULL; // bad detectors full list => fully flag detectors
	long *bolo_bad_80_tot = NULL; // valid worst detectors full list => more than 80% flag detectors


	// for user report at the end of the program
	int *format, *format_tot=NULL;
	long *n_hund, *n_hund_tot=NULL;
	long *n_heig, *n_heig_tot=NULL;
	long *nb_Nan, *nb_Nan_tot=NULL;
	double *sampling_freq, *sampling_freq_tot=NULL;
	long *t_gaps, *t_gaps_tot = NULL;

	t_gaps= new long[samples_struct.ntotscan];
	format= new int[samples_struct.ntotscan];
	n_hund= new long[samples_struct.ntotscan];
	n_heig= new long[samples_struct.ntotscan];
	nb_Nan= new long[samples_struct.ntotscan];
	sampling_freq= new double[samples_struct.ntotscan];

	fill(t_gaps, t_gaps + samples_struct.ntotscan, 0);
	fill(format, format + samples_struct.ntotscan, 0);
	fill(n_hund, n_hund + samples_struct.ntotscan, 0);
	fill(n_heig, n_heig + samples_struct.ntotscan, 0);
	fill(nb_Nan, nb_Nan + samples_struct.ntotscan, 0);
	fill(sampling_freq, sampling_freq + samples_struct.ntotscan, 0.0);



	if(rank==0){ // only for MPI_reduce

		bolo_bad_tot= new long [ndet0];
		bolo_bad_80_tot= new long [ndet0];
		fill(bolo_bad_tot, bolo_bad_tot + ndet0 ,0);
		fill(bolo_bad_80_tot, bolo_bad_80_tot + ndet0 ,0);

#ifdef PARA_FRAME
		t_gaps_tot= new long[samples_struct.ntotscan];
		format_tot= new int[samples_struct.ntotscan];
		n_hund_tot= new long[samples_struct.ntotscan];
		n_heig_tot= new long[samples_struct.ntotscan];
		nb_Nan_tot= new long[samples_struct.ntotscan];
		sampling_freq_tot= new double[samples_struct.ntotscan];

		fill(t_gaps_tot, t_gaps_tot + samples_struct.ntotscan, 0); //
		fill(format_tot, format_tot + samples_struct.ntotscan, 0); //
		fill(n_hund_tot, n_hund_tot + samples_struct.ntotscan, 0);
		fill(n_heig_tot, n_heig_tot + samples_struct.ntotscan, 0);
		fill(nb_Nan_tot, nb_Nan_tot + samples_struct.ntotscan, 0); //
		fill(sampling_freq_tot, sampling_freq_tot + samples_struct.ntotscan, 0.0); //
#endif

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
			std::vector<long> indice;
			double Populated_freq=0.0;

			long init_flag_num=0;
			long end_flag_num=0;

			// initialize default values
			check_struct.Check_it.checkDEC=1;
			check_struct.Check_it.checkRA=1;
			check_struct.Check_it.checkREFERENCEPOSITION=1;
			check_struct.Check_it.checkOFFSETS=1;

#ifdef DEBUG
			cout << endl << endl << "[" << rank <<  "] Checking : " << samples_struct.fitsvect[ii] << endl << endl;
#endif

			format[ii]=test_format(dir.data_dir + samples_struct.fitsvect[ii]); // format = 1 => HIPE, else Sanepic
			if(format[ii]==0){/* fits format indicator : 1 HIPE, 2 Sanepic */
				cerr << "input fits file format is undefined : " << samples_struct.fitsvect[ii] << " . Skipping...\n";
				continue;
			}

			read_bolo_list(dir.data_dir + samples_struct.fitsvect[ii], bolo_fits, ndet_fits); // read fits file detector list
			if(check_bolos(bolo_fits, bolo_fits_0)){ // compare to the first input fits file to ensure compatibility between scans
				cout << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run only together scans that correspond to the same field\n";
				continue;
			}

			std::vector<string> det_vect=bolo_list[ii];
			long ndet_vect = (long)det_vect.size();

			check_detector_is_in_fits(det_vect, ndet_vect, bolo_fits, dir.data_dir + samples_struct.fitsvect[ii]); // check wether used detector user list is correct

			bolo_bad = new long[ndet_fits]; // this scan bad detectors list
			bolo_bad_80 = new long[ndet_fits]; // this scan valid worst detectors list
			percent_tab = new double[ndet_fits];
			fill(bolo_bad,bolo_bad+ndet_fits,0);
			fill(bolo_bad_80,bolo_bad_80+ndet_fits,0);
			fill(percent_tab,percent_tab+ndet_fits,0.0);

#ifdef DEBUG
			cout << "\n[" << rank <<  "] Checking presence of common HDU and position HDU\n";
#endif

			return_value+=check_commonHDU(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii],ndet_fits,check_struct.Check_it); // check presence of channels, time, signal and mask HDUs
			return_value+=check_positionHDU(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii],ndet_fits, format[ii],check_struct.Check_it); // check presence of reference positions and offsets HDUs
			if(format[ii]==1){ // check RA/DEC table presence for HIPE format
#ifdef DEBUG
				cout << "[" << rank <<  "] HIPE format found, Checking Alt position HDU presence\n";
#endif
				return_value+=check_altpositionHDU(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii],ndet_fits,check_struct.Check_it);
			}


			if(return_value<0){
				cerr << "Some Mandatory HDUs are missing in : " << dir.data_dir + samples_struct.fitsvect[ii] << " . Skipping...\n";
				continue;
			}


			if(((format[ii]==2)&&((!check_struct.Check_it.checkREFERENCEPOSITION)||(!check_struct.Check_it.checkOFFSETS))) ||
					((format[ii]==1)&&((!check_struct.Check_it.checkRA)||(!check_struct.Check_it.checkDEC)))){
				cout << "NO POSITION TABLES in : " << dir.data_dir + samples_struct.fitsvect[ii] << " ... Skipping...\n";
			}

			if(check_struct.checkNAN){
#ifdef DEBUG
				cout << "\n[" << rank <<  "] Checking NANs in common HDU and position HDU\n"; // check non-flag NANs presence in whole tables
#endif
				nb_Nan[ii] += check_NAN_commonHDU(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, ndet_fits,check_struct.Check_it);
				nb_Nan[ii] += check_NAN_positionHDU(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, ndet_fits,check_struct.Check_it);
			}

			if((format[ii]==1)&&(check_struct.checkNAN)){ // check NANs presence in RA/DEc tables for HIPE format
#ifdef DEBUG
				cout << "[" << rank <<  "] HIPE format found, Checking NANs in Alt position HDU\n";
#endif
				nb_Nan[ii] += check_NAN_altpositionHDU(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, ndet_fits,check_struct.Check_it);
			}

			if(check_struct.checktime){
#ifdef DEBUG
				cout << "\n[" << rank <<  "] Checking time gaps in time table\n"; // check for time gaps in time table
#endif
				check_time_gaps(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii], proc_param.fsamp, indice, Populated_freq, check_struct.Check_it);
				t_gaps[ii] = (long)indice.size();
			}else{
				Populated_freq=proc_param.fsamp;
			}

			sampling_freq[ii]=Populated_freq;

			if(check_struct.checkGain){
				//				cout << "\n[" << rank <<  "] Computing and Checking bolometer gain correction in signal table\n"; // check for time gaps in time table
				//				check_bolo_gain(dir.data_dir + samples_struct.fitsvect[ii],samples_struct.nsamples[ii], bolo_gain_filename, det, check_struct.Check_it); // TODO : uncomment if needed
				//				getchar();
			}

			if(check_struct.checkflag){
				// Lookfor fully or more than 80% flagged detectors, also flag singletons
#ifdef DEBUG
				cout << "\n[" << rank <<  "] Checking flagged detectors\n"; // check for time gaps in time table
#endif
				check_flag(dir.data_dir + samples_struct.fitsvect[ii],bolo_fits,ndet_fits, samples_struct.nsamples[ii], bolo_bad, n_hund[ii],bolo_bad_80, n_heig[ii],percent_tab, init_flag_num, end_flag_num, check_struct.Check_it);

				if((init_flag_num>0) && (init_flag_num<samples_struct.nsamples[ii])){
#ifdef DEBUG
					cout << "\nInitial Flagged samples will be removed from " << dir.data_dir + samples_struct.fitsvect[ii] << ".\n Considering " << init_flag_num << " samples...\n\n";
#endif
				}else
					init_flag_num=0;

				if((end_flag_num>0) && (end_flag_num<samples_struct.nsamples[ii])){
#ifdef DEBUG
					cout << "Final Flagged samples will be removed from " << dir.data_dir + samples_struct.fitsvect[ii] << ".\n Considering " << end_flag_num << " samples...\n\n";
#endif
				}else
					end_flag_num=0;

				// generating log files :
				outname = dir.output_dir + "bolo_totally_flagged_" + FitsBasename(samples_struct.fitsvect[ii]) +".txt";
#ifdef DEBUG
				cout << "\nWriting informations in :\n" << outname << endl;
#endif
				log_gen(bolo_bad,outname, bolo_fits_0,ndet0); // generate bad detectors log file


				outname = dir.output_dir + "bolo_80_percent_flagged_" + FitsBasename(samples_struct.fitsvect[ii]) +".txt";
#ifdef DEBUG
				cout << "Writing informations in :\n" << outname << endl;
#endif
				log_gen(bolo_bad_80, outname, bolo_fits_0, ndet0, percent_tab); // generate valid worst detectors log file
			}

			if(print_to_bin_file(dir.tmp_dir, samples_struct.fitsvect[ii], init_flag_num, end_flag_num, Populated_freq, indice))
				return 1;

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
			delete [] percent_tab;

		}
	}

#ifdef PARA_FRAME
	MPI_Reduce(t_gaps, t_gaps_tot, samples_struct.ntotscan, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(format, format_tot, samples_struct.ntotscan, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(n_hund, n_hund_tot, samples_struct.ntotscan, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(n_heig, n_heig_tot, samples_struct.ntotscan, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nb_Nan, nb_Nan_tot, samples_struct.ntotscan, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(sampling_freq, sampling_freq_tot, samples_struct.ntotscan, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	t_gaps_tot=t_gaps;
	format_tot=format;
	n_hund_tot=n_hund;
	n_heig_tot=n_heig;
	nb_Nan_tot=nb_Nan;
	sampling_freq_tot=sampling_freq;
#endif

	if(rank==0){

		cout << endl;

		if(check_struct.checkflag){
			// generating log files :
			outname = dir.output_dir + "bolo_totally_flagged.txt";
#ifdef DEBUG
			cout << "\nWriting informations in :\n" << outname << endl;
#endif
			log_gen(bolo_bad_tot,outname, bolo_fits_0, ndet0); // generate bad detectors log file


			outname = dir.output_dir + "bolo_80_percent_flagged.txt";
#ifdef DEBUG
			cout << "\nWriting informations in :\n" << outname << endl;
#endif
			log_gen(bolo_bad_80_tot, outname, bolo_fits_0, ndet0); // generate valid worst detectors log file
		}

		/* -------------------------- Screen output report ------------------------------ */

		cout << "Report :" << endl << endl;

		for(int ii=0;ii<samples_struct.ntotscan;ii++){

			cout << samples_struct.fitsvect[ii] << " : " << endl << endl;

			if(format_tot[ii]==1)
				cout << "- Hipe Format" << endl;
			else
				cout << "- Sanepic Format" << endl;



			cout << "- " << (nb_Nan_tot[ii]>0 ? StringOf(nb_Nan_tot[ii]) : "No") << " unflagged NaNs" << endl;

			cout << "- Sampling Freq. : " << sampling_freq_tot[ii] << " Hz" << endl;

			cout << "- " << t_gaps_tot[ii] << " time gaps" << endl;

			cout << "- " << n_hund_tot[ii] << " 100% flagged bolometers" << endl;

			cout << "- " << n_heig_tot[ii] << " 80% flagged bolometers" << endl;

			cout << endl;

			/* -------------------------------------------------------------------------- */
		}

		cout << "\nPlease run saneFix\n";

	}

	//clean up
	//	delete [] samples_struct.nsamples;
	delete [] bolo_bad_tot;
	delete [] bolo_bad_80_tot;

	delete [] t_gaps;
	delete [] format;
	delete [] n_hund;
	delete [] n_heig;
	delete [] nb_Nan;
	delete [] sampling_freq;

#ifdef PARA_FRAME
	delete [] t_gaps_tot;
	delete [] format_tot;
	delete [] n_hund_tot;
	delete [] n_heig_tot;
	delete [] nb_Nan_tot;
	delete [] sampling_freq_tot;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	if(rank==0)
		cout << "\nEnd of saneCheck\n";

	return EXIT_SUCCESS;
}
