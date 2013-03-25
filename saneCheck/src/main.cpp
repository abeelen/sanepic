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

	struct param_sanePos Pos_param;
	struct param_saneProc Proc_param;
	struct param_sanePic Pic_param;
	struct param_saneInv Inv_param;
	struct param_sanePS PS_param;
	struct param_saneCheck Check_param;

	string outname; // Ouput log files name
	string parser_output = "";
	string bolo_gain_filename="";

	uint32_t parsed=0x0000; // parser error status

	uint32_t mask_saneCheck = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM | BOLOFILE_NOT_FOUND | FSAMP_PROBLEM | FITS_FILELIST_NOT_FOUND; // 0x411f

	uint32_t compare_to_mask; // parser error status

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
		parsed=parse_saneCheck_ini_file(argv[1], parser_output, dir, samples_struct, Pos_param, Proc_param,
				PS_param, Inv_param, Pic_param, Check_param, size, rank);

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
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_param);
		print_param_saneCheck(Check_param);
	}

#ifdef USE_MPI

	// Extra BCast for saneCheck structure...
	MPI_Bcast(&Check_param.checkNAN,  1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Check_param.checkTime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Check_param.checkGain, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Check_param.checkFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
	int *testExt, *testExt_tot=NULL;

	long *nFullyFlagged, *nFullyFlagged_tot =NULL;
	long *nMostlyFlagged, *nMostlyFlagged_tot=NULL;
	long *nNan,     *nNan_tot=NULL;
	long *nTimeGaps, *nTimeGaps_tot = NULL;

	long nFrames = samples_struct.iframe_max-samples_struct.iframe_min;

	double *sampFreqs , *sampFreqs_tot=NULL;


	long *nSpeed,             *nSpeed_tot=NULL;
	double *meanSpeeds,   *meanSpeeds_tot=NULL;
	double *aboveSpeeds, *aboveSpeeds_tot=NULL;
	double *belowSpeeds, *belowSpeeds_tot=NULL;
	int *speedIndices;

	int ntotscan = samples_struct.ntotscan;

	return_value   = new    int[ntotscan];
	testExt        = new    int[ntotscan];
	sampFreqs      = new double[ntotscan];
	meanSpeeds     = new double[ntotscan];
	aboveSpeeds    = new double[ntotscan];
	belowSpeeds    = new double[ntotscan];
	nSpeed         = new   long[ntotscan];


	nTimeGaps      = new long[ntotscan];
	nFullyFlagged  = new long[ntotscan];
	nMostlyFlagged = new long[ntotscan];
	nNan           = new long[ntotscan];

	fill(return_value,   return_value   + ntotscan, 0);
	fill(nTimeGaps,      nTimeGaps      + ntotscan, 0);
	fill(nFullyFlagged,  nFullyFlagged  + ntotscan, 0);
	fill(nMostlyFlagged, nMostlyFlagged + ntotscan, 0);
	fill(nNan,           nNan           + ntotscan, 0);

	fill(testExt,        testExt        + ntotscan, 0);
	fill(sampFreqs ,     sampFreqs      + ntotscan, 0.0);
	fill(meanSpeeds,     meanSpeeds     + ntotscan, 0.0);
	fill(aboveSpeeds,   aboveSpeeds     + ntotscan, 0.0);
	fill(belowSpeeds,   belowSpeeds     + ntotscan, 0.0);
	fill(nSpeed,         nSpeed         + ntotscan, 0);


	for (long iframe = samples_struct.iframe_min +  floor(bolo_rank*nFrames*1.0/bolo_size); iframe < samples_struct.iframe_min + floor((bolo_rank+1)*nFrames*1.0/bolo_size); iframe++){

		long ns = samples_struct.nsamples[iframe];
		string fits_filename = samples_struct.fitsvect[iframe];

		std::vector<string> bolo_fits;
		long ndet_fits;

		long *bolo_flag;
		double *percent_tab;


		std::vector<long> gapIndices;

		double computedFsamp=0.0;
		double meanSpeed=0.0, aboveSpeed = -1., belowSpeed = -1.;

		long init_flag_num=0;
		long end_flag_num=0;

		// initialize default values
		Check_param.Check_it.checkLAT=1;
		Check_param.Check_it.checkLON=1;
		Check_param.Check_it.checkREFERENCEPOSITION=1;
		Check_param.Check_it.checkOFFSETS=1;

#ifdef DEBUG
		cout << endl << endl << "[" << rank <<  "] Checking : " << fits_filename << endl << endl;
#endif

		testExt[iframe]=testExtensions(fits_filename);

		if ( ( (testExt[iframe] & EXT_NECESSARY )   != EXT_NECESSARY) || \
				( ( (testExt[iframe] & HIPE_FORMAT)    != HIPE_FORMAT)  && \
						( (testExt[iframe] & SANEPIC_FORMAT) != SANEPIC_FORMAT) )  ){
			// They are missing extensions that are necessary so skip this file...
			continue;
		}

		read_bolo_list(fits_filename, bolo_fits, ndet_fits); // read fits file detector list

		std::vector<string> det_vect=samples_struct.bolo_list[iframe];
		long ndet_vect = (long)det_vect.size();

		check_detector_is_in_fits(det_vect, ndet_vect, bolo_fits, fits_filename); // check whether used detector user list is correct

		if(Check_param.checkNAN){

			nNan[iframe] += check_NAN_commonHDU(fits_filename,ns,bolo_fits, ndet_fits,Check_param.Check_it);

			switch(testExt[iframe] & EXT_POS) {
			case SANEPIC_FORMAT: // Check RefPos/offsets tables only
				nNan[iframe] += check_NAN_positionHDU(fits_filename,ns,bolo_fits, ndet_fits,Check_param.Check_it);
				break;
			case BOTH_FORMAT: // Check RefPos/offsets tables
				nNan[iframe] += check_NAN_positionHDU(fits_filename,ns,bolo_fits, ndet_fits,Check_param.Check_it);
				// and also ...
				/* no break */
			case HIPE_FORMAT: // Check lon/lat tables
				nNan[iframe] += check_NAN_altpositionHDU(fits_filename,ns,bolo_fits, ndet_fits,Check_param.Check_it);
				break;
			}

		}

		if(Check_param.checkTime && ( (testExt[iframe] & EXT_TIME) == EXT_TIME) ){
#ifdef DEBUG
			cout << "\n[" << rank <<  "] Checking time gaps in time table\n"; // check for time gaps in time table
#endif
			check_time_gaps(fits_filename,ns, samples_struct.fsamp[iframe], gapIndices, computedFsamp, Check_param.Check_it);
			nTimeGaps[iframe] = (long)gapIndices.size();
		}else{
			computedFsamp     = samples_struct.fsamp[iframe];
			nTimeGaps[iframe] = -1;
		}

		sampFreqs [iframe]=computedFsamp;

		if(Check_param.checkGain){
			//				cout << "\n[" << rank <<  "] Computing and Checking bolometer gain correction in signal table\n"; // check for time gaps in time table
			//				check_bolo_gain(fits_filename,ns, bolo_gain_filename, det, check_struct.Check_it); // TODO : uncomment if needed
			//				getchar();
		}

		if(Check_param.checkFlag && ( (testExt[iframe] & EXT_MASK) == EXT_MASK) ){
			percent_tab = new double[ndet_fits];
			fill(percent_tab,percent_tab+ndet_fits,0.0);

			// Look for fully or more than 80% flagged detectors, also flag singletons

			check_Flag(fits_filename,bolo_fits,ndet_fits, ns, percent_tab, init_flag_num, end_flag_num);

#ifdef DEBUG
			if((init_flag_num>0) && (init_flag_num<ns)){
				cout << "\nInitial Flagged samples will be removed from " << fits_filename << ".\n Considering " << init_flag_num << " samples...\n\n";
			} else
				init_flag_num=0;

			if((end_flag_num>0) && (end_flag_num<ns)){
				cout << "Final Flagged samples will be removed from " << fits_filename << ".\n Considering " << end_flag_num << " samples...\n\n";
			}else
				end_flag_num=0;
#endif

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
			outname = dir.output_dir + FitsBasename(fits_filename) +"_fully_flagged.bolo";
			log_gen(bolo_flag,outname, bolo_fits, ndet_fits, percent_tab);

			// ... next for mostly wrong channels... (more than 80% flagged)
			for (long ii=0; ii< ndet_fits; ii++)
				if (percent_tab[ii] >= 80 && percent_tab[ii] < 99){
					nMostlyFlagged[iframe] += 1;
					bolo_flag[ii] = 1;
				} else {
					bolo_flag[ii] = 0;
				}

			outname = dir.output_dir + FitsBasename(fits_filename) +"_mostly_flagged.bolo";
			log_gen(bolo_flag, outname, bolo_fits, ndet_fits, percent_tab);

			// ... next for good channels ... (less than 80% flagged)
			for (long ii=0; ii< ndet_fits; ii++)
				if (percent_tab[ii] >= 80)
					bolo_flag[ii] = 0;
				else
					bolo_flag[ii] = 1;
			outname = dir.output_dir + FitsBasename(fits_filename) +"_good.bolo";
			log_gen(bolo_flag, outname, bolo_fits, ndet_fits);


			// ... next for all channels ...
			for (long ii=0; ii< ndet_fits; ii++)
				bolo_flag[ii] = 1;
			outname = dir.output_dir + FitsBasename(fits_filename) +"_all.bolo";
			log_gen(bolo_flag,outname, bolo_fits, ndet_fits, percent_tab);


			delete [] bolo_flag;
			delete [] percent_tab;


		}

		speedIndices = new int[ns];
		fill(speedIndices, speedIndices+ns, 0);

		if(Check_param.checkSpeed && ( (testExt[iframe] & EXT_TIME) == EXT_TIME ) && ( ((testExt[iframe] & HIPE_FORMAT) == HIPE_FORMAT) || ((testExt[iframe] & SANEPIC_FORMAT) == SANEPIC_FORMAT) ) ) {

			long initFlagSpeed = 0, endFlagSpeed = 0;
			nSpeed[iframe] = check_Speed(Check_param, fits_filename, testExt[iframe] , ns, bolo_fits[0], meanSpeed, belowSpeed, aboveSpeed, speedIndices, initFlagSpeed, endFlagSpeed);

			meanSpeeds[iframe]  =  meanSpeed;
			aboveSpeeds[iframe] = aboveSpeed;
			belowSpeeds[iframe] = belowSpeed;

			if (initFlagSpeed > init_flag_num)
				init_flag_num = initFlagSpeed;
			if (endFlagSpeed > end_flag_num)
				end_flag_num = endFlagSpeed;

		}

		if( print_to_bin_file(dir.tmp_dir, fits_filename, init_flag_num, end_flag_num, computedFsamp, ns, speedIndices, gapIndices) )
			return 1;

		if ( exportCheckToFits(dir.tmp_dir, fits_filename, init_flag_num, end_flag_num, computedFsamp, ns, speedIndices, gapIndices) )
			return 1;

		delete [] speedIndices;

	}

#ifdef USE_MPI

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){ // only for MPI_reduce

		testExt_tot        = new int[ntotscan];
		return_value_tot   = new int[ntotscan];

		nTimeGaps_tot      = new long[ntotscan];
		nFullyFlagged_tot  = new long[ntotscan];
		nMostlyFlagged_tot = new long[ntotscan];
		nNan_tot           = new long[ntotscan];

		nSpeed_tot         = new long[ntotscan];
		meanSpeeds_tot     = new double[ntotscan];
		aboveSpeeds_tot    = new double[ntotscan];
		belowSpeeds_tot    = new double[ntotscan];

		sampFreqs_tot      = new double[ntotscan];

		fill(testExt_tot,        testExt_tot        + ntotscan, 0);
		fill(return_value_tot,   return_value_tot   + ntotscan, 0);
		fill(nTimeGaps_tot,      nTimeGaps_tot      + ntotscan, 0);
		fill(nFullyFlagged_tot , nFullyFlagged_tot  + ntotscan, 0);
		fill(nMostlyFlagged_tot, nMostlyFlagged_tot + ntotscan, 0);
		fill(nNan_tot,           nNan_tot           + ntotscan, 0);

		fill(nSpeed_tot,		nSpeed_tot		+ ntotscan, 0);
		fill(meanSpeeds_tot,	meanSpeeds_tot	+ ntotscan, 0.0);
		fill(aboveSpeeds_tot,	aboveSpeeds_tot	+ ntotscan, 0.0);
		fill(belowSpeeds_tot,	belowSpeeds_tot	+ ntotscan, 0.0);

		fill(sampFreqs_tot,      sampFreqs_tot      + ntotscan, 0.0);
	}

	// Eventually use the custom MPI_Node_Reduce...

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(testExt,             testExt_tot,   ntotscan, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(return_value,   return_value_tot,   ntotscan, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(nTimeGaps,      nTimeGaps_tot,      ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nFullyFlagged,  nFullyFlagged_tot,  ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nMostlyFlagged, nMostlyFlagged_tot, ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(nNan,           nNan_tot,           ntotscan, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(nSpeed,			 nSpeed_tot,	ntotscan, MPI_LONG,		MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(meanSpeeds,	 meanSpeeds_tot,	ntotscan, MPI_DOUBLE,	MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(aboveSpeeds,	aboveSpeeds_tot,	ntotscan, MPI_DOUBLE,	MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(belowSpeeds,	belowSpeeds_tot,	ntotscan, MPI_DOUBLE,	MPI_SUM, 0, MPI_COMM_WORLD);


	MPI_Reduce(sampFreqs,      sampFreqs_tot,      ntotscan, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#else
	testExt_tot			= testExt;
	return_value_tot	= return_value;

	nTimeGaps_tot		= nTimeGaps;
	nFullyFlagged_tot	= nFullyFlagged;
	nMostlyFlagged_tot	= nMostlyFlagged;
	nNan_tot			= nNan;

	nSpeed_tot			= nSpeed;
	meanSpeeds_tot		= meanSpeeds;
	aboveSpeeds_tot		= aboveSpeeds;
	belowSpeeds_tot		= belowSpeeds;

	sampFreqs_tot		= sampFreqs;
#endif


	if(rank==0){

		cout << endl;
		/* -------------------------- Screen parser_output report ------------------------------ */

		cout << "Report :" << endl;

		for(int iframe=0;iframe<ntotscan;iframe++){

			cout << endl << samples_struct.fitsvect[iframe] << " : " << endl ;

			if ( ( (testExt_tot[iframe] & EXT_NECESSARY ) != EXT_NECESSARY) || ( ( (testExt_tot[iframe] & HIPE_FORMAT) != HIPE_FORMAT)  &&  ( (testExt_tot[iframe] & SANEPIC_FORMAT) !=  SANEPIC_FORMAT) ) ){
				// They are missing extensions that are necessary so report...
				cout << "EE - missing necessary fits extensions (";
				if ((testExt_tot[iframe] & EXT_SIGNAL) == 0)
					cout << " signals";
				if ((testExt_tot[iframe] & EXT_MASK) == 0)
					cout << " masks";
				if ((testExt_tot[iframe] & EXT_CHAN) == 0)
					cout << " channels";
				if ((testExt_tot[iframe] & EXT_POS) == 0)
					cout << " positions (either RefPos/offset or Lon/lat)";
				cout << ")" << endl;
			}

			cout << "\t- ";
			switch(testExt_tot[iframe] & EXT_POS){
			case SANEPIC_FORMAT:
				cout << "RefPos/offset";
				break;
			case HIPE_FORMAT:
				cout << "Lon/Lat";
				break;
			case BOTH_FORMAT:
				cout << "Hybrid";
				break;
			default:
				cout << "Unknown";
				break;
			}
			cout << " Format" << endl;

			if	(Check_param.checkNAN) {
				cout << "\t- " << (nNan_tot[iframe]>0 ? StringOf(nNan_tot[iframe]) : "No") << " unflagged NaNs" << endl;
			}

			if (Check_param.checkTime) {
				cout << "\t- Sampling Freq. : " << sampFreqs_tot[iframe] << " Hz" << endl;
				cout << "\t- " << nTimeGaps_tot[iframe] << " time gaps" << endl;
			}

			if (Check_param.checkFlag){
				cout << "\t- " << nFullyFlagged_tot[iframe] << " fully flagged bolometers" << endl;
				cout << "\t- " << nMostlyFlagged_tot[iframe] << " mostly flagged bolometers" << endl;
			}

			if (Check_param.checkSpeed){
				cout << "\t- meanSpeed :" << meanSpeeds_tot[iframe] << " arcsec/sec" << endl;
				cout << "\t- " << nSpeed_tot[iframe] << " flags below " << belowSpeeds_tot[iframe] << " arcsec/sec or above "<< aboveSpeeds_tot[iframe] << " arcsec/sec ("<< nSpeed_tot[iframe]*100./samples_struct.nsamples[iframe] << " %)" << endl;
			}
		}

		cout << endl;

		/* -------------------------------------------------------------------------- */

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
	delete [] testExt;
	delete [] nFullyFlagged;
	delete [] nMostlyFlagged;
	delete [] nNan;
	delete [] sampFreqs ;
	delete [] nSpeed;
	delete [] meanSpeeds;
	delete [] aboveSpeeds;
	delete [] belowSpeeds;

#ifdef USE_MPI

	if(rank==0){

		delete [] nTimeGaps_tot;
		delete [] testExt_tot;
		delete [] nFullyFlagged_tot ;
		delete [] nMostlyFlagged_tot;
		delete [] nNan_tot;
		delete [] sampFreqs_tot;
		delete [] nSpeed_tot;
		delete [] meanSpeeds_tot;
		delete [] aboveSpeeds_tot;
		delete [] belowSpeeds_tot;
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
