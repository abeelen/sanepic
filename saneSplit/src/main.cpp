
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneSplit.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>


extern "C" {
#include "nrutil.h"
}




using namespace std;

//! Usage function, print to stdout a standard command line that explains how to run correctly saneMerge
/*!
 * \param name A char* containing program name
 */
void usage(char *name)
/* You must give,  as an input, the inifile, the fits filename, and one or more pairs of -m -M options */
{
	cerr << "USAGE: " << name << " inifile.ini [-f <filename>] [-m <min time>] [-M <max time>]" << endl;
	cerr << "USAGE: You can use multiple -m and -M options, each -m followed by a -M\n";

}


/*!
 *  This is organized as :
 *
 *  - parse user command line
 *  - check for existence of directory/files given in the command line
 *  - Print parser output to screen
 *
 *  - Get time limits (min and max) to determine which part of the fits file samples are kept
 *  - Read time table in fits file and check whether time limits are correct
 *
 *  - Get fits file format (HIPE or SANEPIC) and channel list
 *  - Determines samples indexes (min and max) corresponding to time limits (min and max)
 *
 *  - Create the output fits file
 *  - Copy input fits file HDU by HDU, creating a new table in the output fits for each, and without copying samples outside of samples indexes
 *
 */

int main(int argc, char *argv[])
{

	/* This project is able to split a fits file (Sanepic or Hipe format) into multiple fits files */

	int parsed = -1;

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_common dir; /* structure that contains output input temp directories */
	string fname; /* fits filename variable */
	std::vector< double > min_time, max_time; /* vectors to store time limits for the output fits files */
	string output = "";

	int m_count = 0, f_count = 0; /* counter for m, f and M options */
	int mM_count = 0; /* counter for m, f and M options */
	int format_fits=0; /* 0 = Hipe format, 1 = SanePic format */
	double *time=NULL; /* input fits file time vector */
	double time_min=0, time_max=0; /* input file min and max time */
//	struct detectors det; /* A structure that contains everything about the detectors names and number */
	std::vector<string> det;
	long ndet;

	printf("\nBeginning of saneSplit:\n\n");


	if (argc<2) /* not enough argument */
		parsed=1;
	else{
		// Parse ini file
		parsed=parse_saneSplit_ini_file(argv[1], output, dir);

		// print parser warning and/or errors
		cout << endl << output << endl;
	}

	if (parsed>0){ /* error during parsing phase */
		switch (parsed){

		case 1: printf("Please run %s using a *.ini file\n",argv[0]);
		usage(argv[0]);
		break;

		case 2 : printf("Wrong program options or argument. Exiting !\n");
		break;

		default :;
		}

		return(EX_CONFIG);
	}


	int retval;
	// read -f, -m and -M options
	while ( (retval = getopt(argc, argv, "f:m:M:")) != -1) {
		switch (retval) {
		case 'f': /* read the fits file name and the number of samples */
			samples_struct.fitsvect.push_back(optarg);
			readFrames(dir.data_dir, samples_struct.fitsvect, samples_struct.nsamples);
#ifdef DEBUG
			cout << "Scan.            : " << samples_struct.fitsvect[0] << endl;
			cout << "Containing.      : " << samples_struct.nsamples[0] << " samples. " << endl;
#endif
			f_count++;
			break;
		case 'm': /* bottom limit for a new fits file extracted from the -f fitsfile */
			min_time.push_back(atof(optarg));
			m_count++;
			mM_count++;
			break;

		case 'M': /* top limit for a new fits file extracted from the -f fitsfile */
			max_time.push_back(atof(optarg));
			mM_count--;
			break;

		default:;
		}
	}


	/************** wrong options  ***************/

	if((f_count>1)||(f_count==0)){
		cout << "\nError ! You must give 1 file AND only one ! Exiting\n";
		usage(argv[0]);
		return (EX_CONFIG);
	}

	// if no m or no M, return
	if((m_count==0)||(mM_count!=0)){
		cout << "\nError ! You must give at least 1 min and max value AND each -M must follow a -m option ! Exiting\n";
		usage(argv[0]);
		return (EX_CONFIG);
	}


	// if m<0 or M<m, return
	for(int ii=0;ii<m_count;ii++){
		if((min_time[ii]<0)||(max_time[ii]-min_time[ii])<=0){
			cout << "Warning : you must provide a crescent order of Positives cut limits : M>m !\nExiting\n";
			return (EX_CONFIG);
		}
	}

	/*********************************************/

#ifdef DEBUG
	cout << setprecision(20) << "bot time limit.  : " << min_time[0] << "\ntop time limit.  : " << max_time[0] << endl;
#endif

	read_time_from_fits(dir.data_dir + samples_struct.fitsvect[0], time, samples_struct.nsamples[0]);
	time_min=time[0];
	time_max=time[samples_struct.nsamples[0]-1];


	int ii=0;
	/* checking that user has given decent time limits according to original fits file time vector */
	while(ii<m_count){
		int indic=0;
		if(max_time[ii]>time_max){
			max_time[ii]=time_max;
			indic++;
		}

		// if given min value (-M) is shorter than the min time (fits), set m to this min value
		if(min_time[ii]<time_min){
			min_time[ii]=time_min;
			indic++;
		}

		if(indic==2){
			cout << "Warning : You are trying to split a file into the exact same file... \nSwitching to next cut limits...\n";
			max_time.erase (max_time.begin()+ii);
			min_time.erase (min_time.begin()+ii);
			m_count--;
			ii--;
		}
		ii++;
	}

	/* in case the time limits were not correct and if no more limits have been given by user : EXIT */
	if(m_count==0){
		cout << "Warning : no more cut limits to apply ! Exiting ...\n";
		return (EX_CONFIG);
	}



	/* get input fits file format : Sanepic or HIPE */
	format_fits=test_format(dir.data_dir + samples_struct.fitsvect[0]);
	if(format_fits==0){
		cerr << "input fits file format is undefined : " << dir.data_dir + samples_struct.fitsvect[ii] << " . Exiting...\n";
	}

	/* read the bolo list in the fits file */
	read_bolo_list(dir.data_dir + samples_struct.fitsvect[0], det, ndet);


	int status=0; /* fits error status number */
	fitsfile *fptr; /* input fits file pointer */
	fitsfile *outfptr; /* output fits file pointer */

	std::ostringstream oss; // we need to store the string in a stringstream because of numbers min_time, max_time
	string fname2 = FitsBasename(samples_struct.fitsvect[0]) + "_split_";

	fname=dir.data_dir + samples_struct.fitsvect[0];
	for(int ii=0; ii < m_count ; ii++){ // for each correct time limits

		// generate a name for the output fits files
		oss << dir.output_dir << fname2 << setprecision(14) << min_time[ii] << "_" << max_time[ii] << ".fits";
		std::string temp = oss.str();

#ifdef DEBUG
		cout << "\nCreating output file :\n" << temp << endl;
#endif

		long min_sample=0; /* bottom and top samples index according to given time limits */
		long max_sample=0;

		// find samples index using time :
		if(max_time[ii]==time[samples_struct.nsamples[0]-1])
			max_sample=samples_struct.nsamples[0]-1;
		else{
			long jj=0;
			while((jj<samples_struct.nsamples[0])){
				if(max_time[ii]<time[jj]){
					max_sample=jj;
					break;
				}
				jj++;
			}
		}

		if(min_time[ii]==time[0])
			min_sample=0;
		else{
			long jj=samples_struct.nsamples[0]-1;
			while(jj>=0){
				if(min_time[ii]>time[jj]){
					min_sample=jj;
					break;
				}
				jj--;
			}
		}


#ifdef DEBUG
		// print to std
		long ns_final= max_sample - min_sample+1; // total number of samples in output file
		cout << "min sample -> max sample : " << min_sample << " -> " << max_sample << endl;
		cout << "Total number of samples : " << ns_final << endl;
#endif

		temp = "!" + temp;

		// open input fits and create output fits
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		if (fits_create_file(&outfptr, temp.c_str(), &status))
			fits_report_error(stderr, status);

		// Copy primary Header
		fits_copy_header(fptr, outfptr, &status);




		if(format_fits==1){ // HIPE format
			cout << "HIPE format found\n";

			// 1 signal
			copy_signal(fptr, outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample, det, ndet);

			// 2 LON 3 LAT
			copy_LON_LAT(fptr,outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample, det, ndet);

			// 4 mask
			copy_mask(fptr, outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample, det, ndet);

			// 5 time
			copy_time(fptr, outfptr, time, min_sample, max_sample);

			// 6 channels
			copy_channels(fptr, outfptr);

			// 7 ref pos
			copy_ref_pos(fptr, outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample);

			// 8 offsets
			copy_offsets(fptr, outfptr);





		}else{ // sanepic format
			cout << "SANEPIC format found\n";


			// 1 ref pos
			copy_ref_pos(fptr,outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample);

			// 2 offsets
			copy_offsets(fptr, outfptr);

			// 3 channels
			copy_channels(fptr, outfptr);

			// 4 time
			copy_time(fptr, outfptr, time, min_sample, max_sample);

			// 5 signal
			copy_signal(fptr, outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample, det, ndet);

			// 6 mask
			copy_mask(fptr, outfptr, dir.data_dir + samples_struct.fitsvect[0], min_sample, max_sample, det, ndet);
		}

		// close both fits files
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		if (fits_close_file(outfptr, &status))
			fits_report_error(stderr, status);

		oss.str("");
		cout << endl;
	}

	// clean up
	delete [] time;

	cout << "END OF SANESPLIT\n";

	return EXIT_SUCCESS;

}




