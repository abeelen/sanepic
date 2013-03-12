#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_math.h>
#include <sysexits.h>

#include "ImageIO.h"
#include "TemporaryIO.h"
#include "CorrPreprocess.h"
#include "NoCorrPreprocess.h"
#include "ParserFunctions.h"
#include "StructDefinition.h"
#include "MPIConfiguration.h"
#include "StructDefinition.h"
#include "SanePicIO.h"
#include "Crc.h"
#include "InputFileIO.h"
#include "ErrorCode.h"
#include "Utilities.h"

#include "getopt.h"

extern "C" {
#include "wcslib/wcshdr.h"
#include "getdata.h"
}

//#define GD_NO_C99_API 1

#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;


/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *	- Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *  - For each file :
 *      - Generate or clear the dirfile parts that will be filled : fData
 *
 *	- Get input fits META DATA
 *
 *  - Read mapHeader, pixel's indice, masked pixel's indice from disk
 *
 *  - Restore incomplete work with previous saved session if needed
 *
 *	- Compute Checksum for crash recovery (if save_data is ON)
 *
 *	- Compute conjugate gradient until convergence or max iteration is reached
 *		- Compute and save temporary maps (and data if save_data is ON) every iterW iterations
 *		- Print convergence criterias to screen every loop
 *
 *	- Write final map to disk
 *
 *
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


	if(rank == 0)
		cout << endl << "sanePic : conjugate gradient descent" << endl;

	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */

	struct param_saneProc Proc_param; /* A structure that contains user options about preprocessing properties */
	struct param_sanePos Pos_param; /* A structure that contains user options about map projection and properties */
	struct param_common dir; /* structure that contains output input temp directories */

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_sanePS PS_param;
	struct param_sanePic Pic_param;
	struct param_saneInv Inv_param;
	string parser_output = "";

	//	iterw = sanePic writes a temporary fits file (map) to disk each iterw iterations (conjugate gradient)
	int flagon = 0; /*  if at least one sample is rejected, flagon=1 */
	int factdupl = 1; /* map duplication factor */
	long long addnpix = 0; /* number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /* number of pix in box constraint */

	// map making parameters
	long long indpix_size; /* indpix read size */
	long long indpsrc_size; /* indpsrc read size */
	long long npix; /* nn = side of the map, npix = number of filled pixels */

	int nwcs=1;             // We will only deal with one wcs....
	struct wcsprm * wcs;    // wcs structure of the image
	long NAXIS1, NAXIS2;  // size of the image

	char * subheader;       // Additionnal header keywords
	int nsubkeys=0;           //

	double *PNdtot = NULL; /* to deal with mpi parallelization : Projected noised data */
	double *PNd = NULL; // (At N-1 d)
	long long *indpix, *indpsrc; /* pixels indices, mask pixels indices */

	string field; /* actual boloname in the bolo loop */

	// main loop variables
	double *S; /* Pure signal */

	// parallel scheme file
	string fname; /* parallel scheme filename */

	uint16_t mask_sanePic = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | NAPOD_WRONG_VALUE | FSAMP_PROBLEM |
			FHP_PROBLEM | FITS_FILELIST_NOT_FOUND | FCUT_PROBLEM; // 0xc39f

	// parser variables
	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	// Default value
	Pic_param.restore = 0;

	// Using getopt to retrieve command line options
	int c;

	static struct option long_options[] = {
			{"restore",     no_argument,       0, 'r'},
			{0, 0, 0, 0}
	};

	/* getopt_long stores the option index here. */
	int option_index = 0;

	while ( (c = getopt_long (argc, argv, "r", long_options, &option_index) ) != -1) {
		switch (c)
		{
		case 'r':
			Pic_param.restore=1;
			break;
		case '?':
			/* getopt_long already printed an error message. */
			break;
		default:
			abort ();
		}
	}

	if ( (argc - optind) != 1){ // less or more than 1 argument
		if (rank == 0) {
			cout << argc << " " << optind << " : " << argc-optind << endl;
			cerr << "EE - Please run  " << argv[0] << " with a .ini file" << endl;
			cerr << "     or with an optional --restore option" << endl;
		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);
	} else {

		/* parse ini file and fill structures */
		parsed = parser_function(argv[optind], parser_output, dir,
				samples_struct, Pos_param, Proc_param, PS_param, Inv_param,
				Pic_param, size, rank);

		compare_to_mask = parsed & mask_sanePic;

		// print parser warning and/or errors
		if (rank == 0)
			cout << parser_output << endl;

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


	// Start...

	if(rank == 0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_param);

	}

	/* ------------------------------------------------------------------------------------*/


#ifdef USE_MPI
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

	if ( cleanup_dirfile_sanePic(dir.tmp_dir, samples_struct, bolo_rank) ) {
		cerr << "EE - Error in initializing dirfile - Did you run sanePre ?" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

	// Read file sizes once for all
	if (bolo_rank == 0) {
		if ( readFramesFromDirfile(dir.tmp_dir, samples_struct)) {
			cerr << "EE - Error in frame size - Did you run sanePre ?" << endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, EX_CONFIG);
#endif
			return EX_CONFIG;
		}
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
	std::vector<long> nsamples_tot;
	nsamples_tot.assign(samples_struct.ntotscan, 0);

	MPI_Allreduce(&(samples_struct.nsamples[0]), &nsamples_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);

	samples_struct.nsamples = nsamples_tot;
	nsamples_tot.clear();

	//	MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_SUB);
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	//	ostringstream t_stream; // fits files filename string stream
	//	t_stream << rank << "\t ";
	//	for (int ii=0; ii < samples_struct.ntotscan; ii++)
	//		t_stream << " " << samples_struct.nsamples[ii];
	//	t_stream << endl;
	//	cout << t_stream.str();

	if (bolo_rank == 0) {
		if ( readNoiseBinSizeFromDirfile(dir.tmp_dir, samples_struct) ) {
			cerr << "EE - Error in reading Noise sizes - Did you run saneInv ?" << endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, EX_CONFIG);
#endif
			return EX_CONFIG;
		}
	}

#ifdef USE_MPI
	std::vector<long> ndet_tot, nbins_tot;
	ndet_tot.assign(samples_struct.ntotscan,0);
	nbins_tot.assign(samples_struct.ntotscan,0);

	MPI_Barrier(MPI_COMM_NODE);

	MPI_Allreduce(&(samples_struct.ndet[0]),   &ndet_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);
	MPI_Allreduce(&(samples_struct.nbins[0]), &nbins_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);

	samples_struct.ndet  =  ndet_tot;
	samples_struct.nbins = nbins_tot;
	ndet_tot.clear();
	nbins_tot.clear();

	//	MPI_Bcast_vector_long(samples_struct.ndet, 0, MPI_COMM_SUB);
	//	MPI_Bcast_vector_long(samples_struct.nbins, 0, MPI_COMM_SUB);
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	//	read pointing header
	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, &subheader, &nsubkeys, rank)){ // read keyrec file
		cerr << "EE - Error reading saved keyrec -- Did you run sanePre ?" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_IOERR;
	}

	//
	//
	//	if (rank == 0){
	//		print_MapHeader(wcs);
	//		cout << "nkeys" << nsubkeys << endl;
	//		char * hptr;
	//		// .... print it
	//		hptr = subheader;
	//		printf("\n\n Sub Header :\n");
	//		for (int ii = 0; ii < nsubkeys; ii++, hptr += 80) {
	//			printf("%.80s\n", hptr);
	//		}
	//	}

	if (Pos_param.flgdupl)
		factdupl = 2; // default 0 : if flagged data are put in a duplicated map

	// Read indpsrc and checks...

	if (rank == 0){

		if (readIndexCCR(indpsrc_size, npixsrc, indpsrc, dir.tmp_dir)) { // read mask index
			cerr << "EE - Please run sanePos" << endl;

#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		if (indpsrc_size != NAXIS1 * NAXIS2) { // check size compatibility
			if (rank == 0)
				cerr << "EE - indpsrc size is not the right size : Check indpsrc.bin file or run sanePos"
				<< endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		addnpix = samples_struct.ntotscan * npixsrc;

		// read indpix
		if (read_indpix(indpix_size, npix, indpix, dir.tmp_dir, flagon)) { // read map index
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		if (indpix_size != (factdupl * NAXIS1 * NAXIS2 + 2 + addnpix + 1)) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpix size is not the right size : Check Indpix_*.bi file or run sanePos"
				<< endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

	} // rank ==0

#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&npix,         1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&npixsrc,      1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&addnpix,      1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&indpix_size,  1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&indpsrc_size, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank!=0){
		indpix  = new long long[indpix_size];
		indpsrc = new long long[indpsrc_size];
	}

	MPI_Bcast(indpix,  indpix_size,  MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(indpsrc, indpsrc_size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

	if (rank == 0) {
		cout << endl;
		cout << "Map Size         : " << NAXIS1 << " x " << NAXIS2 << " pixels" << endl; // print map size
		cout << "Mem. per process : " << prettyPrintSize((indpix_size+npix*9)*8.) << endl;
	}
	if (bolo_rank == 0) {
		if (getAvailableSystemMemory() < (indpix_size+npix*9)*8.*bolo_size){
			cerr << endl;
			cerr << "WW --------------------------------------------" << endl;
			cerr << "WW - Available physical memory may be too low -" << endl;
			cerr << "WW --------------------------------------------" << endl;

			// #ifdef USE_MPI
			// 			MPI_Abort(MPI_COMM_WORLD, 1);
			// #endif
			// 			return (EX_IOERR);

		}
	}

	/*************************************************************/

#ifdef USE_MPI
	if(rank == 0){	// global (At N-1 D) malloc for mpi
		PNdtot = new double[npix];
		fill(PNdtot, PNdtot+npix, 0.0);
	}
#endif


	if (Pic_param.restore) { // restore incomplete work with previous saved data
		if (rank == 0){

			cout << "II - Checking previous session" << endl;
			struct checksum chk_t, chk_t2;
			compute_checksum(dir, Pos_param, Proc_param, Inv_param, PS_param, Pic_param, samples_struct, npix,
					indpix, indpsrc, indpsrc_size, chk_t);
			read_checksum(dir.tmp_dir, chk_t2, "sanePic"); // read previous checksum
			if (compare_checksum(chk_t, chk_t2)) { // compare them
				cout << "Checksums are different !!! Exiting..." << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CONFIG;
			}
		}

		if(rank == 0){
			read_PNd(PNd, npix, dir.tmp_dir, "PNd.bi");
#ifdef USE_MPI
			for(long ii=0; ii<npix;ii++)
				PNdtot[ii]=PNd[ii];
#else
			PNdtot=PNd;
#endif
		}
		else{
			PNd = new double[npix];
			fill(PNd, PNd+npix, 0.0);
		}

	} else {

		//************************************************************************//
		//************************************************************************//
		//Pre-processing of the data : compute PNdtot !
		//************************************************************************//
		//************************************************************************//

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if (rank == 0)
			cout << endl << "Computing Pre Conditioner..." << endl;

		PNd = new double[npix];
		fill(PNd, PNd+npix, 0.0);


		// loop over the scans
		for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){

			// if there is correlation between detectors
			if (Proc_param.CORRon){
				int pb=0;
				/* write_ftrProcesdata parameters : */
				// var double *S =  NULL
				// indpix = seen pixel indice
				// indpsrc = box crossing constraint removal pixel indices
				// NAXIS1, NAXIS2 = sizes of map (pixels)
				// npix = number of pixels that are seen
				// npixsrc = number of pixels in the mask
				// addnpix = number of added pixels in the map
				// tmp_dir = temporary directory
				// det = bolo names array + number of bolo
				// fhp_pix = filter freq in term of sample
				// ns = number of sample for this scan
				// iframe = scan number : 0=> ntotscan if non-MPI

				pb=write_ftrProcesdata(NULL,Proc_param,samples_struct, Pos_param,dir.tmp_dir, indpix,indpsrc,
						NAXIS1, NAXIS2,npix, npixsrc,addnpix, iframe,
						bolo_rank, bolo_size);

				if(pb>0){
					cout << "Problem in write_ftrProcesdata. Exiting ...\n";
#ifdef USE_MPI
					MPI_Abort(MPI_COMM_WORLD, EX_SOFTWARE);
#endif
					return EX_SOFTWARE;
				}

#ifdef USE_MPI
				if ( samples_struct.parallel_scheme == 0) // bolo case
					MPI_Barrier(MPI_COMM_NODE);
#endif

				/* do_PtNd parameters : */
				// PNd => npix (dimension), initialised to 0.0 : (At N-1 d)
				// det = bolo names array + bolo number
				// fcut_pix = freq threshold noise (in term of samples)
				// fsamp = sampling frequency
				// ns = number of samples in the scan
				// size. cf mpi
				// rank. cf mpi
				// indpix = pixel indice double[NAXIS1*NAXIS2]
				// NAXIS1, NAXIS = taille de la carte (1 coté)
				// npix = total number of filled pixels
				// iframe = scan indice
				// *Mp = Null :
				// *Hits = Null (map hits)

				pb+=do_PtNd(samples_struct, PNd, "fData_", bolo_rank,bolo_size,indpix,
						NAXIS1, NAXIS2,npix,iframe, NULL, NULL);
				// Returns Pnd = (At N-1 d), Mp and hits

				if(pb>0){
					cout << "Problem after do_PtNd. Exiting...\n";
#ifdef USE_MPI
					MPI_Abort(MPI_COMM_WORLD, 1);
#endif
					return -1;
				}


			} else { // No correlation case

				do_PtNd_nocorr(PNd, dir.tmp_dir,Proc_param,Pos_param,
						samples_struct, addnpix, indpix,indpsrc,
						NAXIS1, NAXIS2,npix,npixsrc,
						iframe,NULL,rank,size);
				// fillgaps + butterworth filter + fourier transform and PNd generation

			}

		} // end of iframe loop

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
		PNdtot=PNd;
#endif

	} // end of else (restore = 0)


	if (Pic_param.save_data) {
		if (rank == 0) {
			struct checksum chk_t;
			/* Compute Checsum for crash recovery ! */
			compute_checksum(dir, Pos_param, Proc_param, Inv_param, PS_param, Pic_param, samples_struct, npix,
					indpix, indpsrc, indpsrc_size, chk_t);
			if(write_checksum(dir.tmp_dir, chk_t, "sanePic")){ // write down on disk the checksum values
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CANTCREAT;
			}
		}
	}

	/*  END OF CHECKSUM   */

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	/* END PARAMETER PROCESSING */

	if (rank == 0)
		cout << "Starting Conjugate Gradient Descent... " << endl << endl;

	//////////////////////////////////// Computing of sanePic starts here
	string testfile; // log file to follow evolution of both criteria
	ostringstream temp_stream; // fits files filename string stream

	// inititialisation of the Conjugate gradient with preconditioner
	// see (for a complete description of the following variables, go to section B3 for algorithm) :
	// http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
	double *PtNPmatS= NULL, *PtNPmatStot = NULL, *r, *q, *qtot = NULL, *d, *s; // =NULL to avoid warnings
	double *Mp,	*Mptot = NULL; // =NULL to avoid warnings
	// Mp = M in the paper = preconditioner


	double var0 = 0.0, var_n = 0.0, delta0 = 0.0, delta_n = 0.0, alpha = 0.0; // conjugate gradient convergence criteria
	double delta_o, rtq, beta; // conjugate gradient needed parameters (see paper for a complete description)

	int iter=0; // conjugate gradient loop counter
	long long npixeff; // number of filled pixels
	int idupl = 0;

	// memory allocs
	S = new double[npix];
	fill(S, S + npix, 0.0);

	r = new double[npix];
	q = new double[npix];
	d = new double[npix];
	Mp = new double[npix];
	s = new double[npix];

	if (Pos_param.projgaps || !flagon) {
		npixeff = npix;
	} else {
		npixeff = npix - 1;
	}

	if (Pic_param.restore) {
		if (rank == 0){
			cout << "loading idupl\n";
			load_idupl(dir.tmp_dir, idupl);

			// idupl compatibility
			if((idupl>0) && !Pos_param.flgdupl){
				cerr << "EE - idupl cannot be >0 if flgdupl is False ! Exiting...\n";
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, EX_CONFIG);
#endif
				return EX_CONFIG;
			}
		}
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&idupl,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0) { // malloc only for first processor that reduces the data
		Mptot = new double[npix];
		qtot = new double[npix];
		fill(qtot,qtot+npix,0.0);
		fill(Mptot,Mptot+npix,0.0);
	}
#endif

	// in case flagged pixels are put in a duplicated map
	while(idupl <= Pos_param.flgdupl){

		if (Pic_param.save_data && rank == 0)
			write_PNd(PNdtot,npix,dir.tmp_dir, "PNd.bi");


		if((Pic_param.restore)){ // Restore or ...

			if (rank == 0){
				cout << " EE - Loading session" << endl;

				// fill S, d, r, indpix, npixeff, var_n, delta_n and iter with previously saved on disk values
				load_from_disk(dir.tmp_dir, S, d, r, npixeff, var0, var_n, delta0,
						delta_n, iter, Mp);
			}
#ifdef USE_MPI
			MPI_Bcast(&iter,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(S,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);

			//copy Mp to Mptot and PtNPmatS to PtNPmatStot
			if(rank == 0)
				for(long ii = 0; ii< npix; ii++){
					Mptot[ii] = Mp[ii];
				}
#else
			Mptot = Mp;
#endif

			Pic_param.restore=0; // set to 0 because we don't want to load again on idupl=1 loop !


		} else { // Compute

			PtNPmatS = new double[npix];

			fill(PtNPmatS, PtNPmatS + npix, 0.0);
			fill(Mp, Mp + npix, 0.0);
			fill(r, r + npix, 0.0);
			fill(d, d + npix, 0.0);
			fill(s, s + npix, 0.0);

			for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

				// preconditioner computation : Mp
				if (Proc_param.CORRon) {

					write_tfAS(samples_struct, S, indpix, NAXIS1, NAXIS2, npix,
							Pos_param.flgdupl, iframe, bolo_rank, bolo_size);

#ifdef USE_MPI
					if ( samples_struct.parallel_scheme == 0) // bolo case
						MPI_Barrier(MPI_COMM_NODE);
#endif
					do_PtNd(samples_struct, PtNPmatS, "fPs_",
							bolo_rank, bolo_size, indpix,
							NAXIS1, NAXIS2, npix, iframe, Mp, NULL);

				} else {

					do_PtNPS_nocorr(samples_struct, S, dir, Pos_param.flgdupl,
							indpix, NAXIS1, NAXIS2, npix, iframe,
							PtNPmatS, Mp, NULL, rank, size);
				}

			} // end of iframe loop

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);

			if(rank == 0){
				PtNPmatStot = new double[npix];

				fill(Mptot,Mptot+npix,0.0);
				fill(PtNPmatStot,PtNPmatStot+npix,0.0);
			}
			MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			PtNPmatStot = PtNPmatS; // in case non-MPI : pointer equality
			Mptot = Mp;
#endif

			// inititialisation of the Conjugate gradient with preconditioner
			// see : http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
			if (rank == 0) {

				for (long ii = 0; ii < npixeff; ii++)
					if (Mptot[ii] == 0)
						printf("ERROR: Mp[%ld] has elements = 0\n", ii);
					else
						Mptot[ii] = 1.0 / Mptot[ii]; // M : preconditioner

				for (long ii = 0; ii < npixeff; ii++)
					r[ii] = PNdtot[ii] - PtNPmatStot[ii]; // r = b - Ax

				for (long ii = 0; ii < npixeff; ii++)
					d[ii] = Mptot[ii] * r[ii]; // d = M-1 * r


				delta_n = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					delta_n += r[ii] * d[ii]; // delta_new = rT * d


				var_n = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					var_n += r[ii] * r[ii];


				delta0 = delta_n; // delta_0 <= delta_new
				var0 = var_n;

			}

			delete [] PtNPmatS;

#ifdef USE_MPI
			if(rank == 0)
				delete [] PtNPmatStot;
#endif

		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&var0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(d,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

		//start loop
		//max iter = 2000 (default), but ~100 iterations are required to achieve convergence

		// while i<imax and var_new > epsilon² * var_0 : epsilon² = 1e-10 => epsilon = 1e-5
		while ((iter < Pic_param.itermax) && (((var_n / var0 > Pic_param.tolerance) && (idupl
				|| !Pos_param.flgdupl)) || (!idupl && (var_n / var0 > Pic_param.subtolerance)))) {

			fill(q, q + npixeff, 0.0); // q <= A*d


			for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

				if (Proc_param.CORRon) {

					write_tfAS(samples_struct, d, indpix, NAXIS1, NAXIS2, npix,
							Pos_param.flgdupl, iframe, bolo_rank, bolo_size);

#ifdef USE_MPI
					if ( samples_struct.parallel_scheme == 0) // bolo case
						MPI_Barrier(MPI_COMM_NODE);
#endif
					do_PtNd(samples_struct, q, "fPs_",
							bolo_rank, bolo_size, indpix,
							NAXIS1, NAXIS2, npix, iframe, NULL, NULL);

				} else {

					do_PtNPS_nocorr(samples_struct, d,  dir, Pos_param.flgdupl,
							indpix, NAXIS1, NAXIS2, npix, iframe,
							q, NULL, NULL, rank, size);
				}
			} // end of iframe loop

#ifdef USE_MPI
			if(rank == 0)
				fill(qtot,qtot+npix,0.0);
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			qtot = q;
#endif
			if (rank == 0) {

				rtq = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					rtq += qtot[ii] * d[ii]; // rtq = (dT * q)

				alpha = delta_n / rtq; // alpha <= delta_new / (dT * q)


				for (long ii = 0; ii < npixeff; ii++)
					S[ii] += alpha * d[ii]; // x = x + alpha * d, x = S = signal
			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			// every 10 iterations do ....
			if ((iter % 10) == 0) { // if iter is divisible by 10, recompute PtNPmatStot

				PtNPmatS = new double[npix];
				fill(PtNPmatS, PtNPmatS + npix, 0.0);

				for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

					if (Proc_param.CORRon) {

						write_tfAS(samples_struct, S, indpix, NAXIS1, NAXIS2, npix,
								Pos_param.flgdupl, iframe, bolo_rank, bolo_size);

#ifdef USE_MPI
						if ( samples_struct.parallel_scheme == 0) // bolo case
							MPI_Barrier(MPI_COMM_NODE);
#endif

						do_PtNd(samples_struct, PtNPmatS, "fPs_",
								bolo_rank, bolo_size, indpix,
								NAXIS1, NAXIS2, npix, iframe, NULL, NULL);

					} else {

						do_PtNPS_nocorr(samples_struct, S, dir,	Pos_param.flgdupl,
								indpix, NAXIS1, NAXIS2,	npix, iframe,
								PtNPmatS, NULL, NULL, rank, size);
					}

				} // end of iframe loop

#ifdef USE_MPI
				if(rank == 0) { // malloc only for first processor that reduces the data
					PtNPmatStot = new double[npix];
					fill(PtNPmatStot,PtNPmatStot+npix,0.0);
				}
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
				PtNPmatStot = PtNPmatS;
#endif

				if (rank == 0) {

					for (long ii = 0; ii < npixeff; ii++)
						r[ii] = PNdtot[ii] - PtNPmatStot[ii]; //r = b - Ax
				}

				delete [] PtNPmatS;

#ifdef USE_MPI
				if(rank == 0)
					delete [] PtNPmatStot;
#endif

			} else { // if iter is not divisible by 10 ...
				if (rank == 0) {
					for (long ii = 0; ii < npixeff; ii++)
						r[ii] -= alpha * qtot[ii]; // r = r - alpha * q
				}
			}

			if (rank == 0) {

				for (long ii = 0; ii < npixeff; ii++)
					s[ii] = Mptot[ii] * r[ii]; // s = M-1 * r


				delta_o = delta_n; // delta_old <= delta_new

				delta_n = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					delta_n += r[ii] * s[ii]; // delta_new = rT * s

				var_n = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					var_n += r[ii] * r[ii];

				beta = delta_n / delta_o; // beta = delta_new / delta_old
				for (long ii = 0; ii < npixeff; ii++)
					d[ii] = s[ii] + beta * d[ii]; // d = s + beta * d

				// saving iterated maps
				if (Pic_param.iterw && (iter % Pic_param.iterw) == 0) {

					if ((Pic_param.save_data > 0) && ((iter > 0) || (idupl > 0))){
						write_disk(dir.tmp_dir, d, r, S, npixeff, var0, var_n, delta0, delta_n, iter, idupl, Mptot);
					}


					temp_stream << dir.output_dir << Pic_param.map_prefix << "_" << (idupl>0 ? (std::string)"b" : (std::string)"a") << iter << ".fits";
					fname = temp_stream.str();
					temp_stream.str("");

					// Write the raw maps (without flag and without crossing constrains)
					if (writeRawMapToFits(fname, S, NAXIS1, NAXIS2, indpix,	wcs, subheader, nsubkeys, false))
						cerr << "EE - Error writing map... " << endl;

					// Write the flagged data map
					if (Pos_param.flgdupl)
						if ( writeFlagMapToFits(fname, S, NAXIS1, NAXIS2, indpix, wcs, subheader, nsubkeys, true))
							cerr << "EE - Error Writing Flagged Data map... \n";

					// Write the crossing constrains only maps
					if (addnpix)
						writeCCRMapToFits(fname, S, Mptot, NAXIS1, NAXIS2, indpix, indpsrc, npixsrc, factdupl, samples_struct.ntotscan,
								wcs, subheader, nsubkeys, true);

					//					if(write_fits_inifile(fname, dir, Proc_param, Pos_param,
					//							PS_param, Pic_param, Inv_param)) // write saneProc parameters in naive Map fits file header
					//						cerr << "WARNING ! No history will be included in the file : " << fname << endl;
					//					if( write_fits_inputfile(fname, samples_struct))
					//						cerr << "WARNING ! No input files will be included in the file : " << fname << endl;


				} // end of saving iterated maps

				time_t rawtime;
				time ( &rawtime );

				// Now in ISO Standard
				char mytime[20];
				strftime(mytime,20, "%Y-%m-%dT%X", localtime(&rawtime));

				temp_stream <<  mytime << " -- " << "iter_" << (idupl>0 ? (std::string)"b" : (std::string)"a");
				temp_stream << "= "     << setw(4) << iter;
				temp_stream << " crit= "      << setiosflags(ios::scientific) << setprecision (2) << var_n / var0;
				temp_stream << " crit2= "     << setiosflags(ios::scientific) << setprecision (2) << delta_n / delta0;

				// Output to screen ...
				cout << temp_stream.str() << "\r" << flush;
				//cout << endl; // temp : to test for dirfiles !

				// ... and to a logfile
				ofstream logfile;
				string filename = dir.output_dir + "ConvPic.txt";
				logfile.open(filename.c_str(),  ios::out | ios::app);
				logfile << temp_stream.str() << endl;
				logfile.close();
				temp_stream.str("");

			} // end of if (rank == 0)


#ifdef USE_MPI
			// BCast updated criteria and map d
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++; // i = i +1

		} // end of while loop

		// If the gaps are projected on a second map, redo the preprocessing/writing of fdata,
		// which will now replace the flagged data by data from the signal map

		if ((Pos_param.projgaps || (Pos_param.flgdupl)) && !idupl) {

			fill(PNd, PNd + npix, 0.0); // correct : has to be reset to 0 !

			for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

				if (Proc_param.CORRon) {

					write_ftrProcesdata(S, Proc_param, samples_struct,
							Pos_param, dir.tmp_dir, indpix, indpsrc,
							NAXIS1, NAXIS2, npix, npixsrc, addnpix,
							iframe, bolo_rank, bolo_size);

#ifdef USE_MPI
					if ( samples_struct.parallel_scheme == 0) // bolo case
						MPI_Barrier(MPI_COMM_NODE);
#endif

					do_PtNd(samples_struct, PNd, "fData_",
							bolo_rank, bolo_size, indpix,
							NAXIS1, NAXIS2, npix, iframe, NULL, NULL);

				} else {

					do_PtNd_nocorr(PNd, dir.tmp_dir, Proc_param, Pos_param,
							samples_struct, addnpix, indpix, indpsrc,
							NAXIS1, NAXIS2, npix, npixsrc,
							iframe, S, rank, size);
				}

			} // end of iframe loop

#ifdef USE_MPI
			if(rank == 0)
				fill(PNdtot,PNdtot+npix,0.0);
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			PNdtot = PNd; // useless here ??
#endif
		}

		idupl++;
		iter=0;

	}// end of idupl loop


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif



	//******************************  Write final maps in file ********************************

	if (rank == 0)
		cout << endl << "done." << endl;


	// Finally compute hits map & find chart....

	long *hits,      *hitstot=NULL;
	long *findchart, *findcharttot=NULL;

	long long *samptopix;

	hits  = new long[npix];
	fill(hits,hits+npix,0);

	findchart  = new long[npix];
	fill(findchart,findchart+npix,0);

	std::vector<long>  byteID;
	byteID.assign(samples_struct.ntotscan,1);
	for (unsigned long ii=0; ii < byteID.size(); ii++)
		byteID[ii] = (byteID[ii] << ii);

	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){
		long ns = samples_struct.nsamples[iframe];
		samptopix=new long long [ns];

		std::vector<string> det_vect = samples_struct.bolo_list[iframe];
		long ndet = (long)det_vect.size();

		for (long idet1 = floor(bolo_rank*ndet*1.0/bolo_size); idet1<floor((bolo_rank+1)*ndet*1.0/bolo_size); idet1++){

			string field1 = det_vect[idet1];

			if(readSampToPix(samples_struct.dirfile_pointers[iframe], samples_struct.basevect[iframe], field1, samptopix, ns))
				return 1;
			//compute hit counts
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
				findchart[indpix[samptopix[ii]]] |= byteID[iframe];
			}
		}
		delete [] samptopix;
	}
#ifdef USE_MPI

	if(rank==0){
		hitstot = new long[npix];
		fill(hitstot,hitstot+npix,0);
		findcharttot = new long[npix];
		fill(findcharttot,findcharttot+npix,0);
	}

	MPI_Reduce(hits,           hitstot, npix, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(findchart, findcharttot, npix, MPI_LONG,   MPI_BOR, 0, MPI_COMM_WORLD);

#else
	hitstot = hits;
	findcharttot = findchart;
#endif




	if (rank == 0){
		string fname = dir.output_dir + Pic_param.map_prefix + "_sanePic.fits";

		cout << endl << "Output file      : " << fname << endl;

		if(writeMapsToFits(fname, S, Mptot, addnpix, NAXIS1, NAXIS2,
				indpix, indpsrc, npixsrc, factdupl, samples_struct.ntotscan,
				wcs, subheader, nsubkeys, false)) {
			cerr << "EE - Error in writing the maps into fits file." << endl;
			cerr << "     Exiting ..." << endl;
		}


		if ( writeHitMapToFits(fname, hitstot, findcharttot, addnpix, NAXIS1, NAXIS2,
				indpix, indpsrc, npixsrc, factdupl, samples_struct.ntotscan,
				wcs, subheader, nsubkeys, true) ){
			cerr << "WW - Error in writing the hits map into fits file." << endl;
		}

		if( exportExtraToFits(fname, dir,  samples_struct, Proc_param, Pos_param, PS_param, Pic_param, Inv_param)){
			cerr << "EE - Error in writing extra information into fits file."<< endl;
			cerr << "     Exiting... " << endl;
		}


	}// end of rank == 0

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// Close previously opened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])){
				cerr << "EE - error closing dirfile..." << endl;
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}


	// clean up of conjugate gradient variables
	delete[] r;
	delete[] q;
	delete[] d;
	delete[] Mp;
	delete[] s;
	delete[] PNd;

	// clean up map variable
	delete[] S;
	delete[] hits;
	delete[] indpsrc;
	delete[] indpix;


	wcsvfree(&nwcs, &wcs);


#ifdef USE_MPI
	if (rank == 0 && size > 1){
		delete [] qtot;
		delete [] Mptot;
		delete [] PNdtot;
		delete [] hitstot;
	}

	MPI_Comm_free(&MPI_COMM_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_NODE);
	MPI_Finalize();

#endif


	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;

}

