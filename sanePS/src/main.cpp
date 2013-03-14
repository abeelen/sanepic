#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <sysexits.h>
#include <cmath>
#include <ctime>

#include "ImageIO.h"
#include "TemporaryIO.h"
#include "MPIConfiguration.h"
#include "ParserFunctions.h"
#include "StructDefinition.h"
#include "InputFileIO.h"
#include "Crc.h"
#include "ErrorCode.h"
#include "TodProcess.h"
#include "MapMaking.h"

#include "sanePSIO.h"
#include "EstimPSSteps.h"

#include "getopt.h"

extern "C" {
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>


#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#include "Utilities.h"
#endif

using namespace std;

/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - for each file :
 *      - Generate or clear the dirfile parts that will be filled : fData
 *
 *  - Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *	For each scan :
 *
 *
 */
int main(int argc, char *argv[]){


	int      rank,      size; /* MPI processor rank and MPI total number of used processors */
	int bolo_rank, bolo_size; /* As for parallel scheme */
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
		cout << endl << "sanePS : noise-noise power spectra estimation" << endl;

	struct param_common dir; /*! structure that contains output input temp directories */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_saneProc Proc_param; /*! A structure that contains user options about preprocessing properties */
	struct param_sanePos Pos_param; /*! A structure that contains user options about map projection and properties */

	struct param_sanePS PS_param;
	struct param_saneInv Inv_param;
	struct param_sanePic Pic_param;
	string parser_output = "";


	// map making parameters
	int flagon; /*!  if one sample is rejected, flagon=1 */
	long long indpsrc_size;
	long long indpix_size; /*! indpix read size */
	long long *indpix; // map index
	long long npix; // npix = number of filled pixels

	int nwcs=1;             // We will only deal with one wcs....
	struct wcsprm * wcs;    // wcs structure of the image
	long NAXIS1, NAXIS2;  // size of the image

	char * subheader;       // Additionnal header keywords
	int nsubkeys;           //

	string prefixe; // prefix used for temporary name file creation

	// main loop variables
	double *S = NULL; // signal

	string filename; // for filename output
	std::ostringstream temp_stream; // used to remove sprintf horror


	uint16_t mask_sanePS = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | FSAMP_PROBLEM | NCOMP_WRONG_VALUE | ELL_FILE_NOT_FOUND | MIX_FILE_NOT_FOUND |
			FITS_FILELIST_NOT_FOUND | FCUT_PROBLEM; // 0xdd1f

	//	ncomp = number of noise component to estimate
	//	fcut = cut-off freq : dont focus on freq larger than fcut for the estimation !

	// parser variables
	int indice_argv = 1;
	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	PS_param.restore = 0; //default


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
			PS_param.restore=1;
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
		parsed=parser_function(argv[indice_argv], parser_output, dir, samples_struct, Pos_param, Proc_param,
				PS_param, Inv_param, Pic_param, size, rank);


		compare_to_mask = parsed & mask_sanePS;

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

	if(rank==0){
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


	fill_sanePS_struct(PS_param, samples_struct, dir); // all ranks can do it !

	long long addnpix = 0;
	int factdupl = 1;
	long long npixsrc = 0;
	long long *indpsrc;


	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, &subheader, &nsubkeys, rank)){ // read keyrec file
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return (EX_IOERR);
	}


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

	if (rank == 0){
		cout << endl << "Map Size         : " << NAXIS1 << " x " << NAXIS2 << " pixels" << endl; // print map size
	}


	//First time run S=0, after sanepic, S = Pure signal
	if (PS_param.signame != "") {

		S = new double[npix]; // pure signal
		fill(S,S+npix,0.0);

		// if map argument build S from map
		if (rank == 0){
#ifdef DEBUG
			cout << "Reading map.     : " << PS_param.signame << endl;
#endif

			// read pure signal
			if (read_fits_signal(PS_param.signame, S, indpix, NAXIS1, NAXIS2, wcs)) {
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#else
				return (EX_IOERR);
#endif
			}else{

				// 	Check for unobserved Pixels
				int badPix = 0;
				for (long ii=0; ii<npix; ii++) {
					if ( isnan(S[ii]) || isinf(S[ii]) ){
						badPix++;
						S[ii] = 0.0;
					}
				}
				if (badPix > 0)
					cout << "WW - Some observed pixels fall outside the given map... set to 0" << endl;
			}

		} //rank ==0

#ifdef USE_MPI

		MPI_Barrier(MPI_COMM_WORLD);

		if(rank!=0){
			indpix = new long long[indpix_size];
			S = new double[npix];
			//			fill(S,S+npix,0.0);
		}

		MPI_Bcast(indpix,indpix_size,MPI_LONG_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(S,npix,MPI_DOUBLE,0,MPI_COMM_WORLD); // broadcast it to the other procs

#endif
		wcsvfree(&nwcs, &wcs);

	} // PS_param.signame != ""

	if (PS_param.restore) { // restore incomplete work with previous saved data
		if (rank == 0){
			cout << "Checking previous session\n";
			struct checksum chk_t, chk_t2;
			compute_checksum(dir, Pos_param, Proc_param, Inv_param, PS_param, Pic_param, samples_struct, npix,
					indpix, indpsrc, NAXIS1 * NAXIS2, chk_t);
			read_checksum(dir.tmp_dir, chk_t2, "sanePS"); // read previous checksum
			if (compare_checksum(chk_t, chk_t2)) { // compare them
				cout << "Checksums are different !!! Exiting..." << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CONFIG;
			}
		}
	}

	if (PS_param.save_data) {
		if (rank == 0) {
			struct checksum chk_t;
			/* Compute Checsum for crash recovery ! */
			compute_checksum(dir, Pos_param, Proc_param, Inv_param, PS_param, Pic_param, samples_struct, npix,
					indpix, indpsrc, NAXIS1 * NAXIS2, chk_t);
			if(write_checksum(dir.tmp_dir, chk_t, "sanePS")){ // write down on disk the checksum values
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CANTCREAT;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(rank==0)
		cout << endl << "Noise Power Spectra Estimation started : " << endl;

	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) { // proceed scan by scan

		//TODO : in progress

		// ns = number of samples in the "iframe" scan
		// npix = total number of filled pixels
		// iframe = scan number
		// indpix = pixels index
		// MixMatfile = this file contains the number of components in the common-mode component of the noise
		// and the value of alpha, the amplitude factor which depends on detectors but not on time (see formulae (3) in "Sanepic:[...], Patanchon et al.")


		long      ns = samples_struct.nsamples[iframe];
		double fsamp = samples_struct.fsamp[iframe];
		int    ncomp = PS_param.ncomp;

		string basevect = samples_struct.basevect[iframe];
		DIRFILE * dirfile_pointer = samples_struct.dirfile_pointers[iframe];

		//TODO: dummy fcut for the moment...
		double fcut  = samples_struct.fsamp[iframe];

		std::vector<string> det = samples_struct.bolo_list[iframe];
		long ndet = (long)det.size();

		string field; // detector name in the loops


		long nbins; // number of bins for the averaging

		double *SPref;
		double **commonMode; // ?
		double **P, **N;
		double **Rellth = NULL, **Rellexp = NULL;	// Rellth = theorical covariance matrix, Rellexp = experimental cov mat
		//	Nk = Noise power spectrum , Nell = binned noise PS
		double *ell, *km;						// bins values
		double **mixmat; 						// mixing matrix

		int factdupl = 1;
		if(Pos_param.flgdupl==1) factdupl = 2;


		// Precomputation ...

		//	Butterworth filter Precomputation ...
		double *bfilter;
		bfilter = new double[ns / 2 + 1];
		//TODO : f_lp_pix is hard fixed to 1.0 -> avoid errors and wasted time ?
		butterworth_filter(ns, 1.0, 8, bfilter);


		// Apodization Precomputation ...
		double *apodwind;
		apodwind = apodwindow(ns,Proc_param.napod);
		double factapod= 0.0; // apodization factor
		for (long ii=0;ii<ns;ii++)
			factapod += gsl_pow_2(apodwind[ii]);
		factapod /= ns; // apodization factor

		// If we input a map S
		double *Ps           = NULL;		// Projected map, if any....
		long long *samptopix = NULL;	// sample to pixel projection matrix
		if (S != NULL){
			samptopix = new long long[ns]; // sample to pixel proj matrix
			Ps = new double[ns];
		}


		//TODO: Only bolo_rank 0 to read the file....
		// Read out the ells ...
		std::vector<double> dummy;
		std::string output;
		if ( read_file(output, samples_struct.ell_names[iframe], dummy) )
			return FILE_PROBLEM;
		vDouble2carray(dummy, &ell, &nbins);
		dummy.clear();
		nbins   = nbins-1;

		long nbins2=0;
		while ((ell[nbins2] < fcut) && (nbins2 < nbins)){
			nbins2++;
		}

		// .. transform them into km once for all
		km = new double[nbins];
		for (int ii=0;ii<nbins;ii++)
			km[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0)*ns/fsamp;


		// data buffer allocated only once...
		//	fdata*    = fourier transform
		//	data/flag = raw data
		fftw_plan fftplan;
		fftw_complex *fdata1, *fdata2;
		double *data;
		int *flag;
		fdata1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ns/2+1));
		fdata2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ns/2+1));
		data   = (double *) fftw_malloc(sizeof(double)*ns);
		flag   = new int[ns];


		Rellth  = dmatrix(0,(ndet)*(ndet)-1,0,nbins-1);		// covariance matrix (Theory) = pattern : invertible but not Rellexp
		N       = dmatrix(0,ndet-1,0,nbins-1);				// uncorrelated part of the noise
		P       = dmatrix(0,ncomp-1,0,nbins-1);	            // component power spectra

		commonMode = dmatrix(0,ncomp-1,0,ns-1);

		SPref   = new double[nbins];						// first detector PS to avoid numerical problems


		double * Nell, *Nk;

		Nell = new double[nbins]; // binned noise PS
		Nk   = new double[ns/2+1]; // noise PS

		// sign = sigma of the noise


		// Working Matrices and vectors, set once for all loops....
		// Cov = AtN-1A
		// iCov = inverted AtN-1A
		gsl_matrix *iCov, *Cov;
		gsl_vector *ivec, *uvec;
		Cov  = gsl_matrix_alloc(ncomp, ncomp);
		iCov = gsl_matrix_alloc(ncomp, ncomp);
		uvec = gsl_vector_alloc(ncomp);
		ivec = gsl_vector_alloc(ncomp);

		// One has to initialize the two matrices for each iteration...
		init2D_double(Rellth,0,0, (ndet)*(ndet),nbins ,0.0);  // big
		init2D_double(N,0,0,ndet,nbins,0.0);
		init2D_double(P,0,0,ncomp,nbins,0.0);
		init2D_double(commonMode,0,0,ncomp,ns,0.0);

		// Restore variables...
		int goto_step=0;

		if(PS_param.restore)
			if(restore_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), goto_step, commonMode, N, P, Rellexp, Rellth, SPref, ndet, ncomp, nbins, ns)){
				cout << "ERROR restoring data... Exiting ...\n";
				return 1;
			}


		switch(goto_step){

		case 0: {
			//----------------------------------- READ MIXMAT PART -------------------------------//

			mixmat = dmatrix(0, ndet - 1, 0, ncomp - 1);

			if (bolo_rank == 0) {
				rndInitMixmat(ndet, ncomp, mixmat);

				if ( assignMixMat(samples_struct.mix_names[iframe], det, ncomp, mixmat) )
					return FILE_PROBLEM;

				//				if(readMixmatTxt(samples_struct.mix_names[iframe], ndet, ncomp, mixmat))
				//					return FILE_PROBLEM;
			}

#ifdef USE_MPI
			MPI_Bcast_dmatrix(mixmat, ndet, ncomp, 0, MPI_COMM_NODE);
#endif

			goto_step++;

		}
		/* no break */
		case 1: {
			//----------------------------------- COMMON MODE -------------------------------//

			{
				//*************************** Read data and compute components

				double sign0=1;

				double *sign;
				sign    = new double[ndet];
				fill(sign,sign+ndet,0.0);


				double **commonModeTmp;    // common mode per mode
				commonModeTmp = dmatrix(0,ncomp,0,ns-1);
				init2D_double(commonModeTmp,0,0,ncomp,ns,0.0);

				fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata1, FFTW_ESTIMATE);


				// TODO : avoid this stupid first iteration to normalize by the first detector
				if (bolo_rank == 0){

					long idet = 0;

					field = det[idet];

					if(readDataFromDirfile(dirfile_pointer, basevect, field, data, ns))
						return 1;
					if(readFlagFromDirfile(dirfile_pointer, basevect, field, flag, ns))
						return 1;

					if (S != NULL){
						//Read pointing data
						if(readSampToPix(dirfile_pointer, basevect, field, samptopix, ns))
							return EX_NOINPUT;

						deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,Pos_param.flgdupl,factdupl);

						for(long ii=0;ii<ns;ii++)
							data[ii] = data[ii] - Ps[ii];
					}

					MapMakePreProcessData(data,flag,ns,Proc_param ,bfilter,NULL);


					// should apodisation be part of MapMakePreProcess ? : no, not in Corr_preprocess !
					for (long ii=0;ii<ns;ii++)
						data[ii] = data[ii]*apodwind[ii];

					sign0 = crudeSigma(data, ns);

				}

#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_NODE);
				MPI_Bcast(&sign0, 1, MPI_DOUBLE, 0, MPI_COMM_NODE);
#endif


				// loop over detectors
				for (long idet = floor(bolo_rank*ndet*1.0/bolo_size); idet < floor((bolo_rank+1)*ndet*1.0/bolo_size); idet++){

					field = det[idet];

					if(readDataFromDirfile(dirfile_pointer, basevect, field, data, ns))
						return 1;
					if(readFlagFromDirfile(dirfile_pointer, basevect, field, flag, ns))
						return 1;

					if (S != NULL){
						//Read pointing data
						if(readSampToPix(dirfile_pointer, basevect, field, samptopix, ns))
							return EX_NOINPUT;

						deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,Pos_param.flgdupl,factdupl);

						for(long ii=0;ii<ns;ii++)
							data[ii] = data[ii] - Ps[ii];
					}

					MapMakePreProcessData(data,flag,ns,Proc_param ,bfilter,NULL);


					// should apodisation be part of MapMakePreProcess ? : no, not in Corr_preprocess !
					for (long ii=0;ii<ns;ii++)
						data[ii] = data[ii]*apodwind[ii];

					//       fdata are used in cross power spectrum estimation...
					//       BUT it is done differently than power spectrum estimation WHY ?

					// compute fft and save data to disk for later
					fftw_execute_dft_r2c(fftplan, data, fdata1);

					if(writeFdata(dirfile_pointer, ns, fdata1,  "fData_", idet, basevect, det))
						return EX_DATAERR;

					sign[idet] = crudeSigma(data, ns)/sign0;

					// common mode computation
					for (long iComp=0;iComp<ncomp;iComp++)
						for (long ii=0;ii<ns;ii++)
							commonModeTmp[iComp][ii] += mixmat[idet][iComp]/gsl_pow_2(sign[idet])*data[ii];

				}

#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_NODE);
				MPI_Reduce_dmatrix(commonModeTmp, ncomp, ns , MPI_SUM, 0, MPI_COMM_NODE);
				if (bolo_rank == 0)
					MPI_Reduce(MPI_IN_PLACE, sign, ndet, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_NODE);
				else
					MPI_Reduce(sign, NULL, ndet, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_NODE);
#endif

				fftw_destroy_plan(fftplan);

				if (bolo_rank == 0) {

					for (long iComp=0; iComp<ncomp; iComp++){
						for (long ii= 0 ;ii<ns;ii++){
							if ( isnan(commonModeTmp[iComp][ii]) || isinf(commonModeTmp[iComp][ii]) ){
								cerr << basevect << " : something went wrong in the computation of the common mode " << ncomp << endl;
								exit(EXIT_FAILURE);
							}
						}
					}


					//***************************************************************************

					gsl_matrix_set_zero(Cov);

					/////////// AtN-1A
					for (long jComp=0;jComp<ncomp;jComp++){
						for (long iComp=0;iComp<ncomp;iComp++){
							for (long iDet=0;iDet<ndet;iDet++){
								gsl_matrix_set(Cov,jComp,iComp, \
										gsl_matrix_get(Cov,jComp,iComp)  + mixmat[iDet][jComp] * mixmat[iDet][iComp]/gsl_pow_2(sign[iDet]));
							}
						}
					}


					// invert AtN-1A
					gsl_linalg_cholesky_decomp (Cov);
					for (long iComp=0;iComp<ncomp;iComp++){
						gsl_vector_set_basis(uvec,iComp);
						gsl_linalg_cholesky_solve(Cov, uvec, ivec);
						gsl_matrix_set_row(iCov, iComp, ivec);
					}

					//	printf("noise var det 0 =  %10.15g\n",sign0*sign0);

					for (long ii=0;ii<ns;ii++)
						for (long iComp=0;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<ncomp;jComp++)
								commonMode[iComp][ii] += gsl_matrix_get(iCov,iComp,jComp) * commonModeTmp[jComp][ii]; //common mode * (AtN-1A)-1


					// clean up

				}


				// clean up
				delete [] sign;
				free_dmatrix(commonModeTmp,0,ncomp,0,ns-1);

#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_NODE);
				MPI_Bcast_dmatrix(commonMode, ncomp, ns, 0, MPI_COMM_NODE);
#endif

			}


			goto_step++;

			if(PS_param.save_data && bolo_rank == 0){
				save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), goto_step, commonMode, N, P, Rellexp, Rellth, SPref, ndet, ncomp, nbins, ns);
			}
		}
		/* no break */

		case 2: {
			//----------------------------------- ESTIMATE NOISE PS -------------------------------//

			{
				for (long idet = floor(bolo_rank*ndet*1.0/bolo_size); idet < floor((bolo_rank+1)*ndet*1.0/bolo_size); idet++){

					field = det[idet];

					if(readDataFromDirfile(dirfile_pointer, basevect, field, data, ns))
						return 1;
					if(readFlagFromDirfile(dirfile_pointer, basevect, field, flag, ns))
						return 1;

					//TODO : This computation is already done when computing the common mode ->
					//       reuse the fdata if possible?
					//******************************* subtract signal

					if (S != NULL){
						//Read pointing data
						if(readSampToPix(dirfile_pointer,  basevect, field, samptopix, ns))
							return EX_NOINPUT;

						deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,Pos_param.flgdupl,factdupl);

						for(long ii=0;ii<ns;ii++)
							data[ii] = data[ii] - Ps[ii];

					}

					MapMakePreProcessData(data,  flag, ns, Proc_param, bfilter, NULL);


					for (long ii=0;ii<ns;ii++) {
						// Subtract components
						for (long iComp=0; iComp<ncomp; iComp++) {
							data[ii] -= mixmat[idet][iComp]*commonMode[iComp][ii];
						}
						// apply apodization
						data[ii] *= apodwind[ii];
					}

					// Noise power spectra
					/// measure power spectrum of the uncorrelated part of the noise
					noisepectrum_estim(data,ns,km,(int)nbins,fsamp,NULL,Nell,Nk);

					//TODO : normalization by factapod is also done in noisespectrum_estim ?? DONE TWICE ??

					for (long iBin=0;iBin<nbins;iBin++){
						Rellth[idet*ndet+idet][iBin] += Nell[iBin]/factapod; // uncorrelated part added in covariance matrix ??
						N[idet][iBin] = Nell[iBin]/factapod; // uncorrelated part
					}

				}

				delete [] data;
				delete [] flag;
				delete [] bfilter;

				// clean up
				if (S != NULL){
					delete [] Ps ;
					delete [] samptopix;
				}

#ifdef USE_MPI

				MPI_Barrier(MPI_COMM_NODE);
				MPI_Allreduce_dmatrix(N, ndet, nbins , MPI_SUM, MPI_COMM_NODE);

				// Only bolo_rank needs Rellth at this stage....
				MPI_Reduce_dmatrix(Rellth, (ndet)*(ndet),nbins ,MPI_SUM,0,MPI_COMM_NODE);
#endif

				////*********************** Component power spectra

				if (bolo_rank == 0) {

					double  *commonMode_perMode;
					commonMode_perMode = (double *) fftw_malloc(sizeof(double)*ns);
					fill(commonMode_perMode,commonMode_perMode+ns,0.0);

					for (long iComp=0;iComp<ncomp;iComp++){

						for (long ii=0;ii<ns;ii++)
							commonMode_perMode[ii]=commonMode[iComp][ii] * apodwind[ii];

						noisepectrum_estim(commonMode_perMode,ns,km,(int)nbins,fsamp,NULL,Nell,Nk);

						for (long ii=0;ii<nbins;ii++)
							P[iComp][ii] = Nell[ii]/factapod;
					}

					delete [] commonMode_perMode;


					// Write down the Power spectra of the common modes...
					double *data1d; // buffer used to write down 1d array
					data1d = new double[ncomp*nbins];

					for (long iComp=0; iComp< ncomp; iComp++)
						for (long iBin=0; iBin<nbins; iBin++)
							data1d[iComp*nbins+iBin] = P[iComp][iBin];

					// get filename
					temp_stream << dir.output_dir + "Pinit_" << FitsBasename(samples_struct.fitsvect[iframe]) << ".fits";
					filename= temp_stream.str();
					temp_stream.str("");
					write_psd_tofits(filename.c_str(),ncomp,   nbins,'d', data1d, (char *) "CorrelatedNoise", false);
					write_psd_tofits(filename.c_str(),    1, nbins+1,'d',    ell, (char *) "Frequency", true);
					delete [] data1d;


					for (long idet1=0;idet1<ndet;idet1++)
						for (long idet2=0;idet2<ndet;idet2++)
							for (long iComp=0;iComp<ncomp;iComp++)
								for (long iBin=0;iBin<nbins;iBin++)
									Rellth[idet1*ndet+idet2][iBin] += mixmat[idet1][iComp] * mixmat[idet2][iComp] * P[iComp][iBin]; // add correlated part to covariance matrix


				}

#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_NODE);
				MPI_Bcast_dmatrix(P, ncomp, nbins, 0, MPI_COMM_NODE);
				MPI_Bcast_dmatrix(Rellth, ndet*ndet,nbins, 0, MPI_COMM_NODE);
#endif

				// has to be on all rank...
				delete [] apodwind;

			}

			goto_step++;

			if(PS_param.save_data && bolo_rank == 0)
				save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), goto_step, commonMode, N, P, Rellexp, Rellth, SPref, ndet, ncomp, nbins, ns);

			free_dmatrix(Rellth,0,(ndet)*(ndet)-1,0,nbins-1);

		}
		/* no break */

		case 3: {
			//----------------------------------- ESTIMATE COVMAT of the DATA R_exp -------------------------------//

			{
				Rellexp = dmatrix(0,(ndet)*(ndet)-1,0,nbins-1);		// covariance estimated from signal (measured)
				init2D_double(Rellexp,0,0, (ndet)*(ndet),nbins ,0.0); // big

				/////////////////////////////////////// loop again over detectors
				for (long idet1=0;idet1<ndet;idet1++){

					// read data from disk
					if(readFdata(dirfile_pointer, basevect, det[idet1], "fData_", fdata1, ns ))
						return EX_NOINPUT;

					//	TODO: In principle half the loop (-> idet1) should be sufficient...
					for (long idet2 = floor(bolo_rank*(ndet)*1.0/bolo_size); idet2 < floor((bolo_rank+1)*(ndet)*1.0/bolo_size); idet2++){

						// read data from disk
						if(readFdata(dirfile_pointer, basevect, det[idet2],"fData_", fdata2, ns))
							return EX_NOINPUT;

						noisecrosspectrum_estim(fdata1,fdata2,ns,km,(int)nbins,fsamp,NULL,Nell,Nk);

						for (long iBin=0;iBin<nbins;iBin++){
							Rellexp[idet1*ndet+idet2][iBin] += Nell[iBin]/factapod; // noise cross PS ?
							//							Rellexp[idet2*ndet+idet1][iBin] = Rellexp[idet1*ndet+idet2][iBin];
						}

					}

				}

#ifdef USE_MPI
				MPI_Allreduce_dmatrix(Rellexp, ndet*ndet, nbins, MPI_SUM, MPI_COMM_NODE);
#endif


				delete [] fdata1;
				delete [] fdata2;


				////// normalize to the first detector power spectrum in order to avoid numerical problems
				for (long iBin=0;iBin<nbins;iBin++){
					SPref[iBin] = Rellexp[0][iBin]; // first detector PS

					if (bolo_rank == 0) {
						if ( isnan(SPref[iBin]) || isinf(SPref[iBin]) || SPref[iBin] == 0){
							cerr << "EE - Problem in reference power spectrum : " << det[0] << " in " << basevect << endl;
							return 1;
						}
					}
				}

				for (long idet1=0;idet1<ndet;idet1++)
					for (long idet2=0;idet2<ndet;idet2++)
						for (long iBin=0;iBin<nbins;iBin++)
							Rellexp[idet1*ndet+idet2][iBin] /= SPref[iBin]; // normalize to first detector
				for (long iComp=0;iComp<ncomp;iComp++)
					for (long iBin=0;iBin<nbins;iBin++)
						P[iComp][iBin] /= SPref[iBin]; // normalize common mode part
				for (long idet=0;idet<ndet;idet++)
					for (long iBin=0;iBin<nbins;iBin++)
						N[idet][iBin] /= SPref[iBin]; //normalize uncorrelated part


				//	// write Rellexp to disk and also first guess of parameters
				//	temp_stream << dir.output_dir + "Rellexp_" << basename << ".txt";
				//
				//	// get filename
				//	filename= temp_stream.str();
				//	temp_stream.str("");
				//
				//	fp = fopen(filename.c_str(),"w");
				//	for (long jj=0;jj<nbins;jj++)
				//		for (long ii=0;ii<ndet;ii++)
				//			for (long kk=0;kk<ndet;kk++)
				//				fprintf(fp,"%10.15g \t",Rellexp[ii*ndet+kk][jj]); // cross power spectrum
				//	fprintf(fp,"\n");
				//	fclose(fp);
				//
				//	temp_stream << dir.output_dir + "Ninit_" << basename << ".txt";
				//
				//	// get filename
				//	filename= temp_stream.str();
				//	temp_stream.str("");
				//
				//	fp = fopen(filename.c_str(),"w");
				//	for (long ii=0;ii<ndet;ii++){
				//		for (long jj=0;jj<nbins;jj++)
				//			fprintf(fp,"%10.15g \t",N[ii][jj]); // uncorralated part of the noise
				//		fprintf(fp,"\n");
				//	}
				//	fclose(fp);


				if (bolo_rank == 0) {
					double *data1d; // buffer used to write down 1d array
					data1d = new double[ndet*nbins];

					for (long idet=0; idet< ndet; idet++)
						for (long iBin=0; iBin<nbins; iBin++)
							data1d[idet*nbins+iBin] = N[idet][iBin];

					// get filename
					temp_stream << dir.output_dir + "Ninit_" << FitsBasename(samples_struct.fitsvect[iframe]) << ".fits";
					filename= temp_stream.str();
					temp_stream.str("");
					write_psd_tofits(filename.c_str(),ndet,nbins,'d', data1d, (char *) "Noise", false); //resized uncorrelated part
					write_psd_tofits(filename.c_str(),    1, nbins+1,'d',    ell, (char *) "Frequency", true);
					delete [] data1d;
				}
				//	temp_stream << dir.output_dir + "Pinit_" << basename << ".txt";
				//
				//	// get filename
				//	filename= temp_stream.str();
				//	temp_stream.str("");
				//
				//	fp = fopen(filename.c_str(),"w");
				//	for (long ii=0;ii<ncomp;ii++)
				//		for (long jj=0;jj<nbins;jj++)
				//			fprintf(fp,"%10.15g \t",P[ii][jj]); // common mode part of the noise
				//	fprintf(fp,"\n");
				//	fclose(fp);
				//
				//
				//	temp_stream << dir.output_dir + "Ainit_" << basename << ".txt";
				//
				//	// get filename
				//	filename= temp_stream.str();
				//	temp_stream.str("");
				//
				//	fp = fopen(filename.c_str(),"w");
				//	for (long ii=0;ii<ndet;ii++)
				//		for (long jj=0;jj<ncomp;jj++){
				//			fprintf(fp,"%10.15g \t",mixmat[ii][jj]); // mixing matrix
				//			fprintf(fp,"\n");
				//		}
				//	fclose(fp);

			}

			goto_step++;

			if(PS_param.save_data && bolo_rank == 0 )
				save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), goto_step, commonMode, N, P, Rellexp, Rellth, SPref, ndet, ncomp, nbins, ns);
		}
		/* no break */

		case 4: {

			//*********************** fit component and noise power spectra, and mixing matrix *************//
			//********* Using Expectation/Maximization algorithm

			{

				double tottest=0.0;

				double f;

				double *iN, *Pr, *w;
				double **Rxs, **Rxsq, **RnRxsb, **Rxx, **Rxxq, **Rss, **Rssq, **RnRssb;
				double **Pr2, **AiNA, **Mattmp, **ImDR, **ACq, **Cq, **Wq;
				double **new_P, **new_mixmat, **new_N; // Variable for MPI


				iN = new double[ndet];
				Pr = new double[ncomp];
				w  = new double[nbins2];

				Rxs    = dmatrix(0,ndet-1 ,0,ncomp-1);
				Rxsq   = dmatrix(0,ndet-1 ,0,ncomp-1);
				RnRxsb = dmatrix(0,ndet-1 ,0,ncomp-1);
				Rxx    = dmatrix(0,ndet-1 ,0,ndet-1) ;
				Rxxq   = dmatrix(0,ndet-1 ,0,ndet-1) ;
				Rss    = dmatrix(0,ncomp-1,0,ncomp-1);
				Rssq   = dmatrix(0,ncomp-1,0,ncomp-1);
				RnRssb = dmatrix(0,ncomp-1,0,ncomp*ndet-1);
				Pr2    = dmatrix(0,ncomp-1,0,ncomp-1);
				AiNA   = dmatrix(0,ncomp-1,0,ncomp-1);
				Mattmp = dmatrix(0,ndet-1 ,0,ndet-1) ;
				ImDR   = dmatrix(0,ndet-1, 0,ndet-1) ;
				ACq    = dmatrix(0,ndet-1, 0,ncomp-1);
				Cq     = dmatrix(0,ncomp-1,0,ncomp-1);
				Wq     = dmatrix(0,ndet-1, 0,ncomp-1);

				//// Compute weights
				for (long ii=0;ii<nbins2;ii++)
					w[ii] = (ell[ii+1] - ell[ii])*ns/fsamp;

#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_NODE);
				f = fdsf_MPI(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2, bolo_rank, bolo_size, MPI_COMM_NODE) ;
#else
				f = fdsf(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2) ;
#endif

#ifdef DEBUG
				printf("Pre em:   obj: %10.15g\n", f) ;

				cout << "nbins  " << nbins  << endl;
				cout << "nbins2 " << nbins2 << endl;
				cout << "ndet   " << ndet << endl;
				cout << "ncomp  " << ncomp << endl << endl;
#endif

				for (long iter=1;iter <= PS_param.niter;iter++){

					fill(iN,iN+ndet,0.0);
					fill(Pr,Pr+ncomp,0.0);
					init2D_double(Rxs,   0,0,ndet,  ncomp,0.0);
					init2D_double(Rxsq,  0,0,ndet,  ncomp,0.0);
					init2D_double(RnRxsb,0,0,ndet,  ncomp,0.0);
					init2D_double(Rxx,   0,0,ndet,  ndet,0.0); // big
					init2D_double(Rxxq,  0,0,ndet,  ndet,0.0); // big
					init2D_double(Rss,   0,0,ncomp, ncomp,0.0);
					init2D_double(Rssq,  0,0,ncomp, ncomp,0.0);
					init2D_double(RnRssb,0,0,ncomp, ncomp*ndet,0.0);
					init2D_double(Mattmp,0,0,ndet,  ndet,0.0); // big
					init2D_double(ImDR,  0,0,ndet,  ndet,0.0); // big
					init2D_double(ACq,   0,0,ndet,  ncomp,0.0);
					init2D_double(Cq,    0,0,ncomp, ncomp,0.0);
					init2D_double(Wq,    0,0,ndet,  ncomp,0.0);


					new_P       = dmatrix(0,ncomp-1,0,nbins-1);			// MPI component power spectra
					init2D_double(new_P,0,0,ncomp,nbins,0.0);

					for (long iBin = floor(bolo_rank*nbins2*1.0/bolo_size); iBin < floor((bolo_rank+1)*nbins2*1.0/bolo_size); iBin++){

						for (long idet=0;idet<ndet;idet++)
							iN[idet] = 1.0/N[idet][iBin];

						for (long idet1=0;idet1<ndet;idet1++)
							for (long idet2=0;idet2<ndet;idet2++)
								Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][iBin];

						// Robust wrt Pq=0
						for (long iComp=0;iComp<ncomp;iComp++)
							Pr[iComp] = sqrt(P[iComp][iBin]);

						for (long iComp=0;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<ncomp;jComp++)
								Pr2[iComp][jComp] = Pr[iComp]*Pr[jComp];


						for (long iComp=0;iComp<ncomp;iComp++){
							for (long jComp=0;jComp<ncomp;jComp++){
								AiNA[iComp][jComp] = 0.0;
								for (long idet=0;idet<ndet;idet++)
									AiNA[iComp][jComp] += mixmat[idet][jComp] * mixmat[idet][iComp] * iN[idet];
							}
						}

						for (long iComp=0;iComp<ncomp;iComp++){
							for (long jComp=0;jComp<ncomp;jComp++)
								gsl_matrix_set(Cov,iComp,jComp, Pr2[iComp][jComp] * AiNA[iComp][jComp]);

							gsl_matrix_set(Cov,iComp,iComp, gsl_matrix_get(Cov,iComp,iComp) + 1.0);
						}

						// invert matrix (ncomp x ncomp)
						gsl_linalg_cholesky_decomp (Cov);
						for (long iComp=0;iComp<ncomp;iComp++){
							gsl_vector_set_basis(uvec,iComp);
							gsl_linalg_cholesky_solve(Cov, uvec, ivec);
							gsl_matrix_set_row(iCov, iComp, ivec);
						}

						for (long iComp=0;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<ncomp;jComp++)
								Cq[iComp][jComp] = Pr2[iComp][jComp] * gsl_matrix_get(iCov,iComp,jComp);

						for (long idet=0;idet<ndet;idet++)
							for (long iComp=0;iComp<ncomp;iComp++){
								Wq[idet][iComp] = 0.0;
								for (long jComp=0;jComp<ncomp;jComp++)
									Wq[idet][iComp] += mixmat[idet][jComp] * iN[idet] *Cq[jComp][iComp] ;
							}
						for (long idet1=0;idet1<ndet;idet1++)
							for (long iComp=0;iComp<ncomp;iComp++){
								Rxsq[idet1][iComp] = 0.0;
								for (long idet2=0;idet2<ndet;idet2++)
									Rxsq[idet1][iComp] += Rxxq[idet1][idet2] * Wq[idet2][iComp];
							}

						for (long iComp=0;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<ncomp;jComp++){
								Rssq[iComp][jComp] = Cq[iComp][jComp];
								for (long idet=0;idet<ndet;idet++)
									Rssq[iComp][jComp] += Rxsq[idet][iComp] * Wq[idet][jComp];
							}

						for (long iComp=1;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<iComp;jComp++){
								Rssq[jComp][iComp] = 0.5*(Rssq[jComp][iComp]+Rssq[iComp][jComp]);
								Rssq[iComp][jComp] = Rssq[jComp][iComp];
							}

						// update power spectra
						for (long iComp=0;iComp<ncomp;iComp++)
							new_P[iComp][iBin] = abs(Rssq[iComp][iComp]);

						for (long idet=0;idet<ndet;idet++)
							for (long iComp=0;iComp<ncomp;iComp++)
								RnRxsb[idet][iComp] += w[iBin] * iN[idet]*Rxsq[idet][iComp];

						for (long iComp=0;iComp<ncomp;iComp++)
							for (long idet=0;idet<ndet;idet++)
								for (long jComp=0;jComp<ncomp;jComp++)
									RnRssb[iComp][jComp+idet*ncomp] += w[iBin] * iN[idet] * Rssq[iComp][jComp] ;

					}

#ifdef USE_MPI
					MPI_Barrier(MPI_COMM_NODE);
					//TODO : More a Allgather...
					MPI_Allreduce_dmatrix(new_P,   ncomp, nbins,      MPI_SUM, MPI_COMM_NODE);
					MPI_Allreduce_dmatrix(RnRxsb,  ndet, ncomp,      MPI_SUM, MPI_COMM_NODE);
					MPI_Allreduce_dmatrix(RnRssb, ncomp, ncomp*ndet, MPI_SUM, MPI_COMM_NODE);

#endif
					free_dmatrix(P,0,ncomp-1,0, nbins-1);
					P = new_P;


					// update mixing matrix
					new_mixmat  = dmatrix(0, ndet - 1, 0, ncomp - 1);	// MPI mixing matrix
					init2D_double(new_mixmat,0,0,ndet,ncomp,0.0);
					for (long idet = floor(bolo_rank*ndet*1.0/bolo_size); idet < floor((bolo_rank+1)*ndet*1.0/bolo_size); idet++){

						for (long iComp=0;iComp<ncomp;iComp++){
							gsl_vector_set(uvec,iComp, RnRxsb[idet][iComp] );

							for (long jComp=0;jComp<ncomp;jComp++){
								gsl_matrix_set(Cov,iComp,jComp, RnRssb[iComp][jComp+idet*ncomp]);
							}
						}

						// solving the linear system
						// invert matrix (ncomp x ncomp)
						gsl_linalg_cholesky_decomp (Cov);
						gsl_linalg_cholesky_solve(Cov, uvec, ivec);

						for (long iComp=0;iComp<ncomp;iComp++)
							new_mixmat[idet][iComp] = gsl_vector_get(ivec,iComp);

					}

#ifdef USE_MPI
					MPI_Barrier(MPI_COMM_NODE);
					MPI_Allreduce_dmatrix(new_mixmat, ndet, ncomp, MPI_SUM, MPI_COMM_NODE);
#endif

					free_dmatrix(mixmat, 0, ndet-1, 0, ncomp-1);
					mixmat=new_mixmat;

					// EM Step with respect to N, with the new values of A and P
					new_N       = dmatrix(0,ndet-1,0,nbins-1);			// MPI uncorrelated part of the noise
					init2D_double(new_N, 0, 0, ndet, nbins, 0.0);
					for (long iBin = floor(bolo_rank*nbins2*1.0/bolo_size); iBin < floor((bolo_rank+1)*nbins2*1.0/bolo_size); iBin++) {

						for (long idet = 0;idet<ndet;idet++)
							iN[idet] = 1.0/N[idet][iBin];

						for (long idet1=0;idet1<ndet;idet1++)
							for (long idet2=0;idet2<ndet;idet2++)
								Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][iBin];

						// Robust wrt Pq=0
						for (long iComp=0;iComp<ncomp;iComp++)
							Pr[iComp] = sqrt(P[iComp][iBin]);

						for (long iComp=0;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<ncomp;jComp++)
								Pr2[iComp][jComp] = Pr[iComp]*Pr[jComp];

						for (long iComp=0;iComp<ncomp;iComp++) {
							for (long jComp=0;jComp<ncomp;jComp++){
								AiNA[iComp][jComp] = 0.0;
								for (long idet=0;idet<ndet;idet++)
									AiNA[iComp][jComp] += mixmat[idet][jComp] * mixmat[idet][iComp] * iN[idet];
							}
						}

						for (long iComp=0;iComp<ncomp;iComp++){
							for (long jComp=0;jComp<ncomp;jComp++){
								gsl_matrix_set(Cov, iComp,jComp, Pr2[iComp][jComp] * AiNA[iComp][jComp]);
							}
							gsl_matrix_set(Cov, iComp, iComp, gsl_matrix_get(Cov,iComp,iComp) + 1.0 );
						}

						// invert matrix
						gsl_linalg_cholesky_decomp (Cov);

						for (long iComp=0;iComp<ncomp;iComp++){
							gsl_vector_set_basis(uvec,iComp);
							gsl_linalg_cholesky_solve(Cov, uvec, ivec);
							gsl_matrix_set_row(iCov, iComp, ivec);
							//							for (long jj=0;jj<ncomp;jj++)
							//								iCov[iComp][jj] = gsl_vector_get(ivec, jj);
						}

						for (long iComp=0;iComp<ncomp;iComp++)
							for (long jComp=0;jComp<ncomp;jComp++)
								Cq[iComp][jComp] = Pr2[iComp][jComp] * gsl_matrix_get(iCov,iComp,jComp);

						for (long idet=0;idet<ndet;idet++)
							for (long iComp=0;iComp<ncomp;iComp++){
								ACq[idet][iComp] = 0.0;
								for (long jComp=0;jComp<ncomp;jComp++)
									ACq[idet][iComp] += mixmat[idet][jComp]*Cq[jComp][iComp];
							}

						for (long idet2=0;idet2<ndet;idet2++)
							for (long idet1=0;idet1<ndet;idet1++){
								ImDR[idet1][idet2] = 0.0;
								if (idet1 == idet2)
									ImDR[idet1][idet2] = 1.0;
								for (long iComp=0;iComp<ncomp;iComp++)
									ImDR[idet1][idet2] -= ACq[idet1][iComp]*iN[idet2]*mixmat[idet2][iComp] ;
							}

						for (long idet1=0;idet1<ndet;idet1++)
							for (long idet2=0;idet2<ndet;idet2++){
								Mattmp[idet1][idet2] = 0.0;
								for (long jj=0;jj<ndet;jj++)
									Mattmp[idet1][idet2] += Rellexp[idet1+jj*ndet][iBin] * ImDR[idet2][jj];
							}

						for (long idet1=0;idet1<ndet;idet1++){
							for (long idet2=0;idet2<ndet;idet2++)
								new_N[idet1][iBin] += ImDR[idet1][idet2] * Mattmp[idet2][idet1];
						}

						for (long idet=0;idet<ndet;idet++){
							for (long iComp=0;iComp<ncomp;iComp++){
								new_N[idet][iBin] += ACq[idet][iComp]*mixmat[idet][iComp];
							}
							new_N[idet][iBin] = abs(new_N[idet][iBin]);
						}

					}


#ifdef USE_MPI
					MPI_Barrier(MPI_COMM_NODE);
					MPI_Allreduce_dmatrix(new_N, ndet, nbins, MPI_SUM, MPI_COMM_NODE);
#endif

					free_dmatrix(N, 0, ndet-1, 0, nbins-1);
					N=new_N;


					tottest = 0.0;
					for (long iBin=0;iBin<nbins2;iBin++){
						for (long idet=0;idet<ndet;idet++)
							tottest += N[idet][iBin];
						tottest/=ndet;

						for (long idet=0;idet<ndet;idet++)
							if (N[idet][iBin] < tottest*1e-8)
								N[idet][iBin] = tottest*1e-8;
					}

					///// here is the problem


#ifdef USE_MPI
					f = fdsf_MPI(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2, bolo_rank, bolo_size, MPI_COMM_NODE) ;
#else
					f = fdsf(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2) ;
#endif


					if (bolo_rank == 0){

						if (isnan(f) || isinf(f)) {
							cerr << "Nan after fdsf........." << endl;
							return EX_SOFTWARE;
						}

						//#ifdef DEBUG
						//						if ((iter % 10) ==  0) {
						time_t rawtime;
						time ( &rawtime );
						char mytime[20];
						strftime(mytime,20, "%Y-%m-%dT%X", localtime(&rawtime));
						temp_stream <<  mytime << " -- " << "iframe" << iframe << " iter " << setw(4) << iter;
						temp_stream << " f= "      << setiosflags(ios::scientific) << setprecision (12) << f;
						// Output to screen ...
						cout << temp_stream.str() << "\r" << flush;
						ofstream logfile;
						string filename = dir.output_dir + "ConvPS.txt";
						logfile.open(filename.c_str(),  ios::out | ios::app);
						logfile << temp_stream.str() << endl;
						logfile.close();
						temp_stream.str("");
						//						}
						//#endif
					}

				}



				//cleaning up


				delete [] w;

				gsl_matrix_free(Cov);
				gsl_matrix_free(iCov);
				gsl_vector_free(ivec);
				gsl_vector_free(uvec);


				delete [] iN;
				delete [] Pr;
				free_dmatrix(Rxs,0,ndet-1 ,0,ncomp-1);
				free_dmatrix(Rxsq,0,ndet-1 ,0,ncomp-1);
				free_dmatrix(RnRxsb,0,ndet-1 ,0,ncomp-1);
				free_dmatrix(Rxx,0,ndet-1 ,0,ndet-1);
				free_dmatrix(Rxxq,0,ndet-1 ,0,ndet-1);
				free_dmatrix(Rss,0,ncomp-1,0,ncomp-1);
				free_dmatrix(Rssq,0,ncomp-1,0,ncomp-1);
				free_dmatrix(RnRssb,0,ncomp-1,0,ncomp*ndet-1);
				free_dmatrix(Pr2,0,ncomp-1,0,ncomp-1);
				free_dmatrix(AiNA,0,ncomp-1,0,ncomp-1);
				free_dmatrix(Mattmp,0,ndet-1,0,ndet-1);
				free_dmatrix(ImDR,0,ndet-1,0,ndet-1);
				free_dmatrix(ACq,0,ndet-1,0,ncomp-1);
				free_dmatrix(Cq,0,ncomp-1,0,ncomp-1);
				free_dmatrix(Wq,0,ndet-1,0,ncomp-1);

			}



			goto_step++;

			if(PS_param.save_data && bolo_rank == 0)
				save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), goto_step, commonMode, N, P, Rellexp, Rellth, SPref, ndet, ncomp, nbins, ns);
		}
		/* no break */
		case 5: {
			//----------------------------------- WRITE TO DISK -------------------------------//

			if (bolo_rank == 0){

				// Fixing the indeterminacies.  Is it useful here?
				rescaleAP(mixmat, P, ndet, ncomp, nbins2) ;

				//****************************** Compute covariance matrix from the fitted model

				Rellth  = dmatrix(0,(ndet)*(ndet)-1,0,nbins-1);		// covariance matrix (Theory) = pattern : invertible but not Rellexp
				init2D_double(Rellth,0,0,ndet*ndet,nbins,0.0);

				for (long iBin=0;iBin<nbins2;iBin++){
					for (long idet=0;idet<ndet;idet++){
						Rellth[idet*ndet+idet][iBin] += N[idet][iBin]*SPref[iBin];
						for (long idet2=0;idet2<ndet;idet2++)
							for (long iComp=0;iComp<ncomp;iComp++)
								Rellth[idet*ndet+idet2][iBin] += mixmat[idet][iComp] * mixmat[idet2][iComp] * P[iComp][iBin]*SPref[iBin];
					}
				}


				if (nbins2 < nbins)
					for (long idet=0;idet<ndet;idet++)
						for (long iBin=nbins2;iBin<nbins;iBin++)
							Rellth[idet*ndet+idet][iBin] = Rellexp[idet*ndet+idet][iBin]*SPref[iBin];


				if(write_to_disk(dir.output_dir, samples_struct.fitsvect[iframe], PS_param, det, nbins, ell, mixmat, Rellth,
						Rellexp, N, SPref,P))
					return 1;

				free_dmatrix(Rellth,0,(ndet)*(ndet)-1,0,nbins-1);

			}


			goto_step++;

		}
		break;
		}
		//----------------------------------- END OF ESTIMPS -------------------------------//

		// clean up
		delete [] SPref;
		delete [] ell;
		delete [] km;

		free_dmatrix(Rellexp,0,(ndet)*(ndet)-1,0,nbins-1);
		free_dmatrix(mixmat,0,ndet-1,0,ncomp-1);

		free_dmatrix(commonMode,0,ncomp,0,ns-1);
		free_dmatrix(P,0,ncomp-1,0,nbins-1);
		free_dmatrix(N,0,ndet-1,0,nbins-1);


	}


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (rank == 0)
		cout << endl << "done." << endl;

	fftw_cleanup();

	if (PS_param.signame != "") {
		delete[] S;
		delete[] indpix;
	} else {
		if(rank==0)
			delete[] indpix;
	}

	if(rank==0)
		delete[] indpsrc;

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


#ifdef USE_MPI

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_free(&MPI_COMM_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_NODE);
	MPI_Finalize();
#endif


	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;

}
