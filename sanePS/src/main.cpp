#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <sysexits.h>

#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "estimPS.h"
#include "struct_definition.h"
#include "inputFileIO.h"

extern "C" {
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}


#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;

int main(int argc, char *argv[])
{
	int size;
	int rank;
#ifdef USE_MPI
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
#endif

	if (rank == 0)
		cout << endl << "sanePS :  Noise Power Spectra Estimation" << endl;

	struct param_sanePre proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */
	std::vector<detectors> detector_tab; /*! A structure that contains everything about the detectors names and their number */

	// map making parameters
	int flagon; /*!  if one sample is rejected, flagon=1 */
	long long ind_size; // indpix size
	long long *indpix; // map index
	long NAXIS1, NAXIS2; // map size
	long long npix; // npix = number of filled pixels

	string field; // actual boloname in the bolo loop
	string prefixe; // prefix used for temporary name file creation

	long iframe_min = 0, iframe_max = 0;

	// main loop variables
	double *S = NULL; // signal
	struct param_sanePS structPS;

	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
	std::vector<double> fcut_vector;
	struct param_sanePic struct_sanePic;
	string output = "";

	//	ncomp = number of noise component to estimate
	//	fcut = cut-off freq : dont focus on freq larger than fcut for the estimation !

	int parsed = 0; // parser error code
	if (argc < 2) { // too few arguments
		printf("Please run %s using a *.ini file\n", argv[0]);
		parsed = 1;
	} else {

		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], output, dir, samples_struct, pos_param, proc_param, fcut_vector,
				structPS, struct_sanePic, rank, size);

		// print parser warning and/or errors
		cout << endl << output << endl;
	}

	if (parsed > 0) { // error during parser phase
		if (rank == 0)
			switch (parsed) {

			case 1:
				printf("Please run %s using a *.ini file\n", argv[0]);
				break;

			case 2:
				printf("Wrong program options or argument. Exiting !\n");
				break;

			case 3:
				printf("Exiting...\n");
				break;

			default:
				break;
			}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

	// parser print screen function
	parser_printOut(dir, samples_struct, pos_param,  proc_param,
			structPS, struct_sanePic, rank);


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table = new int[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];

	fill_sanePS_struct(structPS, samples_struct);

	//First time run S=0, after sanepic, S = Pure signal
	if (structPS.signame != "NOSIGFILE") {

		//		long NAXIS1_read=0, NAXIS2_read=0;
		long long addnpix = 0;
		int factdupl = 1;
		//		double *PNdtot;
		int nwcs = 1;
		long long npixsrc = 0;
		long long *indpsrc;
		//		long long npix2=0;
		struct wcsprm * wcs;
		read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2); // read keyrec file
		if (rank == 0)
			cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl << endl; // print map size


		if(pos_param.maskfile!=""){ // in case a mask have been used
			long long test_size;

			if (read_indpsrc(test_size, npixsrc, indpsrc, dir.tmp_dir)) { // read mask index
#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return (EX_IOERR);
			}

			if (test_size != NAXIS1 * NAXIS2) { // check size compatibility
				if (rank == 0)
					cout
					<< "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos"
					<< endl;
#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return (EX_IOERR);
			}
		}

		if (pos_param.flgdupl)
			factdupl = 2; // default 0 : if flagged data are put in a duplicated map

		if (read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon)) { // read map indexes
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return (EX_IOERR);
		}

		if (ind_size != (factdupl * NAXIS1 * NAXIS2 + 2 + addnpix)) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpix size is not the right size : Check Indpix_*.bi file or run sanePos"
				<< endl;
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return (EX_IOERR);
		}

		S = new double[npix]; // pure signal
		fill(S,S+npix,0.0); // TODO : maybe not needed !
		// if map argument build S from map
		if (rank == 0){
			cout << "Reading model map : " << structPS.signame << endl;

			// read pure signal
			if (read_fits_signal(structPS.signame, S, indpix, NAXIS1, NAXIS2, wcs)) {
#ifdef USE_MPI
				MPI_Finalize();
#endif
				return (EX_IOERR);
			}
		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(S,npix,MPI_DOUBLE,0,MPI_COMM_WORLD); // broadcast it to the other procs
#endif

		wcsvfree(&nwcs, &wcs);

#ifdef DEBUG
		FILE * fp;
		fp = fopen("reconstructed_1dsignal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);
		fclose(fp);
#endif


	}


#ifdef USE_MPI

	ofstream file;
	if(samples_struct.scans_index.size()==0) {

		int test=0;
		string fname = dir.output_dir + parallel_scheme_filename;
		if(rank==0)
			cout << "Getting configuration and frame order from file : " << fname << endl;
		test = define_parallelization_scheme(rank,fname,dir.input_dir,samples_struct,size, iframe_min, iframe_max);

		if(test==-1) {
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return EX_SOFTWARE;
		}

	}else{
		int test=0;
		test = verify_parallelization_scheme(rank,dir.output_dir,samples_struct, size, iframe_min, iframe_max);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&test,1,MPI_INT,0,MPI_COMM_WORLD);

		if(test>0){
			MPI_Finalize();
			return EX_SOFTWARE;

		}

	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min) {
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";

	}

	MPI_Barrier(MPI_COMM_WORLD);

#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index, samples_struct.index_table);
#endif

	fill_noisevect(samples_struct);
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);

	if(read_bolo_for_all_scans(detector_tab, dir, samples_struct, rank, size)){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return(EX_IOERR);
	}
	printf("Number of bolometers : \n");
	for(long iframe=0;iframe<samples_struct.ntotscan;iframe++)
		printf("Scan number %ld : %s %ld\n", iframe,(char*)(FitsBasename(samples_struct.fits_table[iframe]).c_str()), detector_tab[iframe].ndet);


	for (long iframe = iframe_min; iframe < iframe_max; iframe++) { // proceed scan by scan

		struct detectors det = detector_tab[iframe];

		if(EstimPowerSpectra(det, proc_param, dir, pos_param, structPS,
				samples_struct, NAXIS1, NAXIS2, npix, iframe, indpix, S, rank)){
			cout << "Error in EstimPowerSpectra procedure. Exiting ...\n";
			break;
		}
		// ns = number of samples in the "iframe" scan
		// npix = total number of filled pixels
		// iframe = scan number
		// indpix = pixels index
		// MixMatfile = this file contains the number of components in the common-mode component of the noise
		// and the value of alpha, the amplitude factor which depends on detectors but not on time (see formulae (3) in "Sanepic:[...], Patanchon et al.")
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	//clean up
	delete[] samples_struct.nsamples;
	delete[] samples_struct.fits_table;
	delete[] samples_struct.index_table;
	delete[] samples_struct.noise_table;

	fftw_cleanup();

	if (structPS.signame != "NOSIGFILE") {
		delete[] S;
		delete[] indpix;
	}

	if (rank == 0)
		cout << endl << "done." << endl;

	return EXIT_SUCCESS;
}
