#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>

#include "InputFileIO.h"
#include "MPIConfiguration.h"
#include "DataIO.h"
#include "SaneCheckParse.h"
#include "ParserFunctions.h"
#include "SaneCheckTools.h"
#include "StructDefinition.h"
#include "ErrorCode.h"

extern "C" {
#include "nrutil.h"
}

//TODO: In the case of saneCheck, one may wish to use a Master/Slave application style
// See http://www.lam-mpi.org/tutorials/one-step/ezstart.php

#ifdef USE_MPI
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


	int      rank,      size; /* MPI processor rank and MPI total number of used processors */
	int  bolo_rank,  bolo_size; /* As for parallel scheme */
	int node_rank, node_size; /* On a node basis, same as *sub* but for frame scheme */

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Comm MPI_COMM_NODE, MPI_COMM_MASTER_NODE;
#else
	size = 1;
	rank = 0;
	bolo_size  = 1;
	bolo_rank  = 0;
	node_size = 1;
	node_rank = 0;
#endif

	if(rank==0)
		cout << endl << "saneCheck: checking fits files consistency" << endl;

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_common dir;  // structure that contains output input temp directories

	struct param_sanePos pos_param;
	struct param_saneProc proc_param;
	struct param_sanePic sanePic_struct;
	struct param_saneInv saneInv_struct;
	struct param_sanePS structPS;
	struct param_saneCheck check_struct;

	string outname; // Ouput log files name
	string parser_output = "";
	string bolo_gain_filename="";

	uint16_t parsed=0x0000; // parser error status

	uint16_t mask_saneCheck = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | FSAMP_PROBLEM | FITS_FILELIST_NOT_FOUND; // 0x411f

	uint16_t compare_to_mask; // parser error status

	if (argc<2) {/* not enough argument */
		if (rank == 0)
			cerr << "EE - Please run  " << StringOf(argv[0]) << " with a .ini file" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);
	} else {
		/* parse ini file and fill structures */
		parsed=parse_saneCheck_ini_file(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, sanePic_struct, check_struct, size, rank);

		compare_to_mask = parsed & mask_saneCheck;

		// print parser warning and/or errors
		if (rank == 0)
			cout << endl << parser_output << endl;

		// in case there is a parsing error or the dirfile format file was not created correctly
		if(compare_to_mask>0x0000){


			switch (compare_to_mask){/* error during parsing phase */

			case 0x0001:
				if (rank==0)
					cerr << " EE - Please run " << StringOf(argv[0]) << " using a correct *.ini file" << endl;
				break;

			default :
				if (rank==0)
					cerr << "EE - Wrong program options or argument. Exiting ! " <<  "("<< hex << compare_to_mask << ")" << endl;
				break;

			}

#ifdef USE_MPI
			MPI_Finalize();
#endif
			return EX_CONFIG;
		}
	}

	/* ------------------------------------------------------------------------------------*/
	// Start...

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, sanePic_struct, saneInv_struct);
		print_saneCheck_ini(check_struct);
	}

#ifdef USE_MPI

	// Extra BCast for saneCheck structure...
	MPI_Bcast(&check_struct.checkNAN,  1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&check_struct.checktime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&check_struct.checkGain, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&check_struct.checkflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if(configureMPI(dir.output_dir, samples_struct, rank, size,
			bolo_rank,  bolo_size, node_rank, node_size,
			MPI_COMM_NODE, MPI_COMM_MASTER_NODE)){
		if (rank==0)
			cerr << endl << endl << "Exiting..." << endl;

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_CONFIG;
	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	/* ------------------------------------------------------------------------------------*/
	parser_output.clear();

	// Create tmpDir if needed (sub_rank_dummy here to prevent concurrence in the frame case on same node
	if (init_tmpdir(parser_output, samples_struct, dir.tmp_dir, node_rank) ){
		cerr << parser_output;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
#endif

	// Read file size once for all
	if (rank == 0)
		readFramesFromFits(samples_struct);
#ifdef USE_MPI
	MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_WORLD);
#endif

	// for user report at the end of the program
	int *return_value, *return_value_tot=NULL;
	int *format, *format_tot=NULL;

	long *nFullyFlagged, *nFullyFlagged_tot =NULL;
	long *nMostlyFlagged, *nMostlyFlagged_tot=NULL;
	long *nNan, *nNan_tot=NULL;
	long *nTimeGaps, *nTimeGaps_tot = NULL;

	double *sampFreqs , *sampFreqs_tot=NULL;

	return_value    = new  int[samples_struct.ntotscan];

	format    = new  int[samples_struct.ntotscan];
	sampFreqs = new double[samples_struct.ntotscan];

	nTimeGaps      = new long[samples_struct.ntotscan];
	nFullyFlagged  = new long[samples_struct.ntotscan];
	nMostlyFlagged = new long[samples_struct.ntotscan];
	nNan           = new long[samples_struct.ntotscan];

	fill(return_value,   return_value + samples_struct.ntotscan, 0);
	fill(nTimeGaps,      nTimeGaps      + samples_struct.ntotscan, 0);
	fill(format,         format         + samples_struct.ntotscan, 0);
	fill(nFullyFlagged,  nFullyFlagged  + samples_struct.ntotscan, 0);
	fill(nMostlyFlagged, nMostlyFlagged + samples_struct.ntotscan, 0);
	fill(nNan,           nNan           + samples_struct.ntotscan, 0);

	fill(sampFreqs ,     sampFreqs      + samples_struct.ntotscan, 0.0);


	if (bolo_rank == 0) {

		for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){
			std::vector<string> bolo_fits;
			long ndet_fits;

			long *bolo_flag;

			double *percent_tab;
			std::vector<long> indice;
			double Populated_freq=0.0;

			long init_flag_num=0;
			long end_flag_num=0;

			// initialize default values
			check_struct.Check_it.checkLAT=1;
			check_struct.Check_it.checkLON=1;
			check_struct.Check_it.checkREFERENCEPOSITION=1;
			check_struct.Check_it.checkOFFSETS=1;

#ifdef DEBUG
			cout << endl << endl << "[" << rank <<  "] Checking : " << samples_struct.fitsvect[iframe] << endl << endl;
#endif

			format[iframe]=test_format(samples_struct.fitsvect[iframe]);
			// 0 = Unknown format,
			// 1 = RefPos & offsets format (sanepic),
			// 2 = lon/lat format (HIPE),
			// 3 = both sanepic & HIPE

			if(format[iframe]==0){
				cerr << "input fits file format is undefined : " << samples_struct.fitsvect[iframe] << " . Skipping...\n";
				continue;
			}

			read_bolo_list(samples_struct.fitsvect[iframe], bolo_fits, ndet_fits); // read fits file detector list

			std::vector<string> det_vect=samples_struct.bolo_list[iframe];
			long ndet_vect = (long)det_vect.size();

			check_detector_is_in_fits(det_vect, ndet_vect, bolo_fits, samples_struct.fitsvect[iframe]); // check wether used detector user list is correct



#ifdef DEBUG
			cout << "\n[" << rank <<  "] Checking presence of common HDU and position HDU\n";
#endif

			return_value[iframe]+=check_commonHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],ndet_fits,check_struct.Check_it); // check presence of channels, time, signal and mask HDUs


			switch(format[iframe]){
			case 1: // Check RefPos/offsets tables only
				return_value[iframe]+=check_positionHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],ndet_fits, format[iframe],check_struct.Check_it); // check presence of reference positions and offsets HDUs
				break;
			case 3: // Check RefPos/offsets tables
				return_value[iframe]+=check_positionHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],ndet_fits, format[iframe],check_struct.Check_it); // check presence of reference positions and offsets HDUs
				// and also ...
			case 2: // Check lon/lat tables
				return_value[iframe]+=check_altpositionHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],ndet_fits,check_struct.Check_it);
				break;
			}


			if(return_value[iframe]<0){
				cerr << "EE - Some Mandatory HDUs are missing in : " << samples_struct.fitsvect[iframe] << " . Skipping...\n";
				continue;
			}


			if(((format[iframe]==1 || format[iframe] == 3)&&((!check_struct.Check_it.checkREFERENCEPOSITION)||(!check_struct.Check_it.checkOFFSETS))) ||
					((format[iframe]==2 || format[iframe] == 3)&&((!check_struct.Check_it.checkLON)||(!check_struct.Check_it.checkLAT)))){
				cerr << "EE - NO POSITION TABLES in : " << samples_struct.fitsvect[iframe] << " ... Skipping...\n";
			}


			if(check_struct.checkNAN){
#ifdef DEBUG
				cout << "\n[" << rank <<  "] Checking NANs in common HDU and position HDU\n"; // check non-flag NANs presence in whole tables
#endif
				nNan[iframe] += check_NAN_commonHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],bolo_fits, ndet_fits,check_struct.Check_it);

				switch(format[iframe]){
				case 1:{ // Check RefPos/offsets tables only
					nNan[iframe] += check_NAN_positionHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],bolo_fits, ndet_fits,check_struct.Check_it);
					break;
				}
				case 3: // Check RefPos/offsets tables
					nNan[iframe] += check_NAN_positionHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],bolo_fits, ndet_fits,check_struct.Check_it);
					// and also ...
				case 2: // Check lon/lat tables
					nNan[iframe] += check_NAN_altpositionHDU(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe],bolo_fits, ndet_fits,check_struct.Check_it);
					break;
				}

			}

			if(check_struct.checktime){
#ifdef DEBUG
				cout << "\n[" << rank <<  "] Checking time gaps in time table\n"; // check for time gaps in time table
#endif
				check_time_gaps(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe], samples_struct.fsamp[iframe], indice, Populated_freq, check_struct.Check_it);
				nTimeGaps[iframe] = (long)indice.size();
			}else{
				Populated_freq=samples_struct.fsamp[iframe];
			}

			sampFreqs [iframe]=Populated_freq;

			if(check_struct.checkGain){
				//				cout << "\n[" << rank <<  "] Computing and Checking bolometer gain correction in signal table\n"; // check for time gaps in time table
				//				check_bolo_gain(samples_struct.fitsvect[iframe],samples_struct.nsamples[iframe], bolo_gain_filename, det, check_struct.Check_it); // TODO : uncomment if needed
				//				getchar();
			}

			if(check_struct.checkflag){
				percent_tab = new double[ndet_fits];
				fill(percent_tab,percent_tab+ndet_fits,0.0);

				// Look for fully or more than 80% flagged detectors, also flag singletons

				check_flag(samples_struct.fitsvect[iframe],bolo_fits,ndet_fits, samples_struct.nsamples[iframe], percent_tab, init_flag_num, end_flag_num);

				if((init_flag_num>0) && (init_flag_num<samples_struct.nsamples[iframe])){
#ifdef DEBUG
					cout << "\nInitial Flagged samples will be removed from " << samples_struct.fitsvect[iframe] << ".\n Considering " << init_flag_num << " samples...\n\n";
#endif
				}else
					init_flag_num=0;

				if((end_flag_num>0) && (end_flag_num<samples_struct.nsamples[iframe])){
#ifdef DEBUG
					cout << "Final Flagged samples will be removed from " << samples_struct.fitsvect[iframe] << ".\n Considering " << end_flag_num << " samples...\n\n";
#endif
				}else
					end_flag_num=0;

				// generating log files...
				bolo_flag   = new long[ndet_fits]; // flag for goodness
				fill(bolo_flag,bolo_flag+ndet_fits,0);

				// ... first for very wrong channels... (more than 99 % flagged)
				for (long ii=0; ii< ndet_fits; ii++)
					if (percent_tab[ii] >= 99){
						nFullyFlagged[iframe] +=1 ;
						bolo_flag[ii] = 1;
					} else {
						bolo_flag[ii] = 0;
					}
				outname = dir.output_dir + FitsBasename(samples_struct.fitsvect[iframe]) +"_fully_flagged.bolo";
				log_gen(bolo_flag,outname, bolo_fits, ndet_fits, percent_tab);

				// ... next for mostly wrong channels... (more than 80% flagged)
				for (long ii=0; ii< ndet_fits; ii++)
					if (percent_tab[ii] >= 80 && percent_tab[ii] < 99){
						nMostlyFlagged[iframe] += 1;
						bolo_flag[ii] = 1;
					} else {
						bolo_flag[ii] = 0;
					}

				outname = dir.output_dir + FitsBasename(samples_struct.fitsvect[iframe]) +"_mostly_flagged.bolo";
				log_gen(bolo_flag, outname, bolo_fits, ndet_fits, percent_tab);

				// ... next for good channels ... (less than 80% flagged)
				for (long ii=0; ii< ndet_fits; ii++)
					if (percent_tab[ii] >= 80)
						bolo_flag[ii] = 0;
					else
						bolo_flag[ii] = 1;
				outname = dir.output_dir + FitsBasename(samples_struct.fitsvect[iframe]) +"_good.bolo";
				log_gen(bolo_flag, outname, bolo_fits, ndet_fits);


				// ... next for all channels ...
				for (long ii=0; ii< ndet_fits; ii++)
					bolo_flag[ii] = 1;
				outname = dir.output_dir + FitsBasename(samples_struct.fitsvect[iframe]) +"_all.bolo";
				log_gen(bolo_flag,outname, bolo_fits, ndet_fits, percent_tab);


				delete [] bolo_flag;
				delete [] percent_tab;


			}

			if(print_to_bin_file(dir.tmp_dir, samples_struct.fitsvect[iframe], init_flag_num, end_flag_num, Populated_freq, indice))
				return 1;

		}
	}


#ifdef USE_MPI


	if(rank==0){ // only for MPI_reduce

		format_tot         = new int[samples_struct.ntotscan];
		return_value_tot   = new int[samples_struct.ntotscan];

		nTimeGaps_tot      = new long[samples_struct.ntotscan];
		nFullyFlagged_tot  = new long[samples_struct.ntotscan];
		nMostlyFlagged_tot = new long[samples_struct.ntotscan];
		nNan_tot           = new long[samples_struct.ntotscan];

		sampFreqs_tot      = new double[samples_struct.ntotscan];

		fill(format_tot,         format_tot         + samples_struct.ntotscan, 0);
		fill(return_value_tot,   return_value_tot   + samples_struct.ntotscan, 0);
		fill(nTimeGaps_tot,      nTimeGaps_tot      + samples_struct.ntotscan, 0);
		fill(nFullyFlagged_tot , nFullyFlagged_tot  + samples_struct.ntotscan, 0);
		fill(nMostlyFlagged_tot, nMostlyFlagged_tot + samples_struct.ntotscan, 0);
		fill(nNan_tot,           nNan_tot           + samples_struct.ntotscan, 0);

		fill(sampFreqs_tot,      sampFreqs_tot      + samples_struct.ntotscan, 0.0);
	}

	// Eventually use the custom MPI_Node_Reduce...

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(format,         format_tot,         samples_struct.ntotscan, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(return_value,   return_value_tot,   samples_struct.ntotscan, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(nTimeGaps,      nTimeGaps_tot,      samples_struct.ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nFullyFlagged,  nFullyFlagged_tot,  samples_struct.ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nMostlyFlagged, nMostlyFlagged_tot, samples_struct.ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nNan,           nNan_tot,           samples_struct.ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(sampFreqs,      sampFreqs_tot,      samples_struct.ntotscan, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#else
	format_tot         = format;
	return_value_tot   = return_value;

	nTimeGaps_tot      = nTimeGaps;
	nFullyFlagged_tot  = nFullyFlagged;
	nMostlyFlagged_tot = nMostlyFlagged;
	nNan_tot           = nNan;

	sampFreqs_tot      = sampFreqs;
#endif

	if(rank==0){

		cout << endl;
		/* -------------------------- Screen parser_output report ------------------------------ */

		cout << "Report :" << endl << endl;

		for(int ii=0;ii<samples_struct.ntotscan;ii++){

			cout << samples_struct.fitsvect[ii] << " : " << endl ;

			if (return_value_tot[ii] == 0 ) {
				switch(format_tot[ii]){
				case 1:
					cout << "- RefPos/offset Format" << endl;
					break;
				case 2:
					cout << "- Lon/Lat Format" << endl;
					break;
				case 3:
					cout << "- Hybrid Format" << endl;
					break;
				default:
					cout << "- Unknown Format" << endl;
					break;
				}


				cout << "- " << (nNan_tot[ii]>0 ? StringOf(nNan_tot[ii]) : "No") << " unflagged NaNs" << endl;
				cout << "- Sampling Freq. : " << sampFreqs_tot[ii] << " Hz" << endl;
				cout << "- " << nTimeGaps_tot[ii] << " time gaps" << endl;
				cout << "- " << nFullyFlagged_tot[ii] << " fully flagged bolometers" << endl;
				cout << "- " << nMostlyFlagged_tot[ii] << " mostly flagged bolometers" << endl;

			} else {
				cout << "- Bad file format (missing HDUs)" << endl;
			}

			cout << endl;

			/* -------------------------------------------------------------------------- */
		}

		cout << "\nPlease run saneFix\n";

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// Close previously openened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])){
				cerr << "EE - error closing dirfile...";
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

	//clean up

	delete [] nTimeGaps;
	delete [] format;
	delete [] nFullyFlagged;
	delete [] nMostlyFlagged;
	delete [] nNan;
	delete [] sampFreqs ;


#ifdef USE_MPI

	if(rank==0){

		delete [] nTimeGaps_tot;
		delete [] format_tot;
		delete [] nFullyFlagged_tot ;
		delete [] nMostlyFlagged_tot;
		delete [] nNan_tot;
		delete [] sampFreqs_tot;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_free(&MPI_COMM_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_NODE);

	MPI_Finalize();

#endif

	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;
}
