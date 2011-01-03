#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <sysexits.h>
#include <cmath>

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


#ifdef PARA_FRAME
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;

int main(int argc, char *argv[])
{
	int size;
	int rank;
#ifdef PARA_FRAME
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
	struct param_saneInv saneInv_struct;
	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
//	std::vector<double> fcut_vector;
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
		parsed=parser_function(argv[1], output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, struct_sanePic);

		if(rank==0)
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

#ifdef PARA_FRAME
//		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

	fill_sanePS_struct(structPS, samples_struct, dir);

	//First time run S=0, after sanepic, S = Pure signal
	if (structPS.signame != "") {

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
#ifdef PARA_FRAME
//				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return (EX_IOERR);
			}

			if (test_size != NAXIS1 * NAXIS2) { // check size compatibility
				if (rank == 0)
					cout
					<< "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos"
					<< endl;
#ifdef PARA_FRAME
//				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return (EX_IOERR);
			}
		}

		if (pos_param.flgdupl)
			factdupl = 2; // default 0 : if flagged data are put in a duplicated map

		if (read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon)) { // read map indexes
#ifdef PARA_FRAME
//			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return (EX_IOERR);
		}

		if (ind_size != (factdupl * NAXIS1 * NAXIS2 + 2 + addnpix)) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpix size is not the right size : Check Indpix_*.bi file or run sanePos"
				<< endl;
#ifdef PARA_FRAME
//			MPI_Barrier(MPI_COMM_WORLD);
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
#ifdef PARA_FRAME
//				MPI_Sendrecv_replace(data_adr,count,datatype, // TODO : dire aux autres rank qu'il faut sortir sinon ils attendent a MPI_BARRIER ...
//				destproc,sendtag,srcproc,recvtag,
//				comm,status_adr)
				cout << "RANK 0 goes out !\n";
				MPI_Finalize();
#endif
				exit (EX_IOERR);
			}

			// 	Check for unobserved Pixels
			int badPix = 0;
			for (long ii=0; ii<npix; ii++) {
				if ( isnan(S[ii]) || isinf(S[ii]) ){
					badPix++;
					S[ii] = 0;
				}
			}
			if (badPix > 0)
				cout << "WW - Some observed pixels fall outside the given map... set to 0" << endl;



#ifdef DEBUG
		FILE * fp;
		fp = fopen("reconstructed_1dsignal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);
		fclose(fp);
#endif

		}

#ifdef PARA_FRAME
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(S,npix,MPI_DOUBLE,0,MPI_COMM_WORLD); // broadcast it to the other procs
#endif

		wcsvfree(&nwcs, &wcs);


	}


#ifdef PARA_FRAME

	if(configure_PARA_FRAME_samples_struct(dir.output_dir, samples_struct, rank, size, iframe_min, iframe_max)){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // ifram_min=iframe_max => This processor will not do anything
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
	}
#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;
#endif


	if(rank==0)
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic);

	for (long iframe = iframe_min; iframe < iframe_max; iframe++) { // proceed scan by scan

		std::vector<string> det;

		string output_read = "";
		if(read_channel_list(output_read, samples_struct.bolovect[iframe], det)){
			cout << output_read << endl;
			return 1;
		}

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

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	//clean up
//	delete[] samples_struct.nsamples;

	fftw_cleanup();

	if (structPS.signame != "NOSIGFILE") {
		delete[] S;
		delete[] indpix;
	}

	if (rank == 0)
		cout << endl << "done." << endl;

	return EXIT_SUCCESS;
}
