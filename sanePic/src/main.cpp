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
#include "getopt.h"

#include "imageIO.h"
#include "temporary_IO.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "parser_functions.h"
#include "struct_definition.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "write_maps_to_disk.h"
#include "crc.h"
#include "inputFileIO.h"

extern "C" {
#include "wcslib/wcshdr.h"
#include "getdata.h"
}

//#define GD_NO_C99_API 1

//#if defined(PARA_BOLO) || defined(PARA_FRAME)
//#define USE_MPI
//#endif

#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;

//**********************************************************************************//
//**********************************************************************************//
//***************** Beginning of conjugate gradient program ************************//
//**********************************************************************************//
//**********************************************************************************//


int main(int argc, char *argv[]) {


	int size;
	int rank;

#ifdef USE_MPI
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank==0)
		printf("\nsanepic_conjugate_gradient:\n\n");

#else
	size = 1;
	rank = 0;

#endif

	if(rank==0)
		printf("\nBeginning of sanePic:\n\n");

	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//

	struct param_sanePre proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */

	int nwcs = 1; // number of wcs : 1
	//	iterw = sanePic writes a temporary fits file (map) to disk each iterw iterations (conjugate gradient)
	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if at least one sample is rejected, flagon=1 */
	int factdupl = 1; /*! map duplication factor */
	long long addnpix = 0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /*! number of pix in box constraint */

	// map making parameters
	//	long long npix2; /*! used to check PNd reading was correct */
	long long indpix_size; /*! indpix read size */
	long long indpsrc_size; /*! indpsrc read size */
	long NAXIS1, NAXIS2; // map dimensions
	long long npix; /*! nn = side of the map, npix = number of filled pixels */

	double *PNdtot = NULL; /*! to deal with mpi parallelization : Projected noised data */
	double *PNd = NULL; // (At N-1 d)
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */


	double f_lppix_Nk, f_lppix; // noise cut-off frequency (in terms of samples number), filter cut-off freq (samples)
	long ns; // number of samples for the considered scan

	string field; /*! actual boloname in the bolo loop */
	string prefixe; /*! prefix used for temporary name file creation */
	struct param_sanePic struct_sanePic;
	//	std::vector<double> fcut; /*! noise cutting frequency vector */

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	string parser_output = "";

	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	string fname; /*! parallel scheme filename */

	uint16_t mask_sanePic = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | NAPOD_WRONG_VALUE | FSAMP_WRONG_VALUE |
			F_LP_WRONG_VALUE | FITS_FILELIST_NOT_FOUND | FCUT_FILE_PROBLEM; // 0xc39f

	std::vector<string> key;
	std::vector<int> datatype;
	std::vector<string> val;
	std::vector<string> com;

	std::vector<std::vector<std::string> > bolo_list; // this vector contains all bolonames for all the scans

	int para_bolo_indice = 0;
	int para_bolo_size = 1;

#ifdef PARA_BOLO

	para_bolo_indice = rank;
	para_bolo_size = size;

#endif

	// parser variables
	int indice_argv = 1;
	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	//	int retval;
	//	while ( (retval = getopt(argc, argv, "r")) != -1) {
	//		switch (retval) {
	//		case 'r': /* read the fits file name and the number of samples */
	//			struct_sanePic.restore = 1;
	//			cout << "restore!\n";
	//			break;
	//		default :
	//			cout << "default!\n";
	//			break;
	//		}
	//	}

	//TODO: getopt or not ??
	if ((argc < 2) || (argc > 3)) // no enough arguments
		compare_to_mask = 0x0001;
	else {
		struct_sanePic.restore = 0; //default

		if (argc == 3) {

			//			int c;
			//
			//			static struct option long_options[] =
			//			{
			//					/* These options don't set a flag.
			//			                  We distinguish them by their indices. */
			//					{"restore",     no_argument,       0, 'r'},
			//					{0, 0, 0, 0}
			//			};
			//			/* getopt_long stores the option index here. */
			//			int option_index = 0;
			//
			//			c = getopt_long (argc, argv, "r",
			//					long_options, &option_index);
			//
			//			if(c!=-1)
			//				switch (c)
			//				{
			//				case 'r':
			//					struct_sanePic.restore = 1;
			//					break;
			//				case '?':
			//					/* getopt_long already printed an error message. */
			//					break;
			//
			//				default:
			//					break;
			//				}
			//
			//
			//			cout << "struct_sanePic.restore : " << struct_sanePic.restore << endl;
			//			cout << "optind : " << optind << " vs argc : " << argc << endl;
			//
			//			while (optind < argc)
			//				printf ("%s\n", argv[optind++]);
			//
			//			getchar();

			struct_sanePic.restore = 1;
			if (strcmp(argv[1], (char*) "--restore") != 0) {
				if (strcmp(argv[2], (char*) "--restore") != 0)
					indice_argv = -1;
				else
					indice_argv = 1;
			} else {
				indice_argv = 2;
			}
		}

		if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors
			if (indice_argv > 0){
				/* parse ini file and fill structures */
				parsed = parser_function(argv[indice_argv], parser_output, dir,
						samples_struct, pos_param, proc_param, structPS, saneInv_struct,
						struct_sanePic, size, rank);

				compare_to_mask = parsed & mask_sanePic;

				// print parser warning and/or errors
				cout << endl << parser_output << endl;

			}else
				compare_to_mask = 0x0001;


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
		}


	}


	string name_rank = dir.output_dir + "debug_sanePic.txt"; // log file name

#ifdef DEBUG

	std::ostringstream oss;
	oss << dir.output_dir + "debug_sanePic_" << rank << ".txt";
	name_rank = oss.str();

	ofstream file_rank;
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	file_rank.open(name_rank.c_str(), ios::out | ios::trunc); //& ios::trunc
	if(!file_rank.is_open()) {
		cerr << "File [" << file_rank << "] Invalid." << endl;
		return EX_CANTCREAT;
	}
	file_rank << "Opening file for writing debug at " << asctime (timeinfo) << endl;
	file_rank.close();

#endif

	/********************* Define parallelization scheme   *******/

#ifdef USE_MPI

	MPI_Datatype message_type;
	struct ini_var_strings ini_v;
	int ntotscan;

	if(rank==0){
		fill_var_sizes_struct(dir, pos_param, proc_param,
				saneInv_struct, structPS, samples_struct, ini_v);

		ntotscan = ini_v.ntotscan;
	}

	MPI_Bcast(&ntotscan, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank!=0){
		ini_v.fitsvect=new int[ntotscan];
		ini_v.noisevect=new int[ntotscan];
		ini_v.bolovect=new int[ntotscan];
	}

	ini_v.ntotscan=ntotscan;

	Build_derived_type_ini_var (&ini_v,	&message_type);

	MPI_Bcast(&ini_v, 1, message_type, 0, MPI_COMM_WORLD);

	commit_struct_from_root(dir, pos_param, proc_param, saneInv_struct, struct_sanePic, structPS, samples_struct, ini_v, rank);

	MPI_Barrier(MPI_COMM_WORLD);
#endif


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


	/* ------------------------------------- READ bolo list ----------------------------*/

	if(channel_list_to_vect_list(samples_struct, bolo_list, rank)){
		cout << "error in channel_list_to_vect_list" << endl;
		return EX_CONFIG;
	}

	/* ------------------------------------------------------------------------------------*/

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic, saneInv_struct);

		cleanup_dirfile_fdata(dir.tmp_dir, samples_struct, bolo_list);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
#endif

	if(get_noise_bin_sizes(dir.tmp_dir, samples_struct, rank)){
#ifdef USE_MPI
		MPI_Finalize();
#endif
		return (EX_IOERR);
	}

	// Open the dirfile to read temporary files
	string filedir = dir.tmp_dir + "dirfile";
	samples_struct.dirfile_pointer = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED | GD_BIG_ENDIAN);

	if (gd_error(samples_struct.dirfile_pointer) != 0) {
		cout << "error opening dirfile : " << filedir << endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return 1;
	}

	// get input fits META DATA
	if(rank==0)
		if(get_fits_META(samples_struct.fitsvect[0], key, datatype, val, com))
			cout << "\nProblem while getting fits META... Continue but the map header will not be full...\n\n";

	//	read pointing informations
	struct wcsprm * wcs;
	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, rank)){ // read keyrec file
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return (EX_IOERR);
	}

	if (pos_param.flgdupl)
		factdupl = 2; // default 0 : if flagged data are put in a duplicated map

	if (rank == 0){
		cout << "Map Size         : " << NAXIS1 << " x " << NAXIS2 << " pixels\n" << endl; // print map size

		if (read_indpsrc(indpsrc_size, npixsrc, indpsrc, dir.tmp_dir)) { // read mask index
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}
		if (indpsrc_size != NAXIS1 * NAXIS2) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos"
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

		if (indpix_size != (factdupl * NAXIS1 * NAXIS2 + 2 + addnpix)) { // check size compatibility
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

	MPI_Bcast(&npix,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&npixsrc,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);

	if(rank!=0){
		// each frame contains npixsrc pixels with index indsprc[] for which
		// crossing constraint are removed
		// thus
		// addnpix = number of pix to add in pixon
		//         = number of scans * number of pix in box crossing constraint removal
		addnpix = samples_struct.ntotscan * npixsrc;

		indpix_size = factdupl * NAXIS1 * NAXIS2 + 2 + addnpix;
		indpsrc_size = NAXIS1 * NAXIS2;

		indpix = new long long[indpix_size];
		indpsrc = new long long[indpsrc_size];
	}

	MPI_Bcast(indpix,indpix_size,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(indpsrc,indpsrc_size,MPI_LONG_LONG,0,MPI_COMM_WORLD);
#endif


	/*************************************************************/

	if (iframe_min < 0 || iframe_min > iframe_max || iframe_max	> samples_struct.ntotscan) {
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting"
				<< endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return (EX_SOFTWARE);
	}


#ifdef USE_MPI
	if(rank==0){	// global (At N-1 D) malloc for mpi
		PNdtot = new double[npix];
		fill(PNdtot, PNdtot+npix, 0.0);
	}
#endif


	if (struct_sanePic.restore>0) { // restore incomplete work with previous saved data
		if (rank == 0){
			cout << "Checking previous session\n";
			struct checksum chk_t, chk_t2;
			//			compute_checksum(argv[indice_argv], dir.tmp_dir, npix, indpix,
			//					indpsrc, indpsrc_size, chk_t); // compute input data checksum to ensure they haven't changed since the previous run
			compute_checksum(dir, pos_param, proc_param, saneInv_struct, structPS, struct_sanePic, samples_struct, npix,
					indpix, indpsrc, indpsrc_size, chk_t); // TODO : test it
			read_checksum(dir.tmp_dir, chk_t2, "sanePic"); // read previous checksum
			if (compare_checksum(chk_t, chk_t2)) { // compare them
				cout << "Checksums are different !!! Exiting..." << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CONFIG;
			}
		}

		if(rank==0)
			read_PNd(PNd, npix, dir.tmp_dir, "PNdCorr.bi");
		else{
			PNd = new double[npix];
			fill(PNd, PNd+npix, 0.0);
		}

	}else {

		//************************************************************************//
		//************************************************************************//
		//Pre-processing of the data : compute PNdtot !
		//************************************************************************//
		//************************************************************************//

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		PNd = new double[npix];
		fill(PNd, PNd+npix, 0.0);

		prefixe = "fdata_"; // Fourier transform of the data file prefixe

		// loop over the scans
		for (long iframe=iframe_min;iframe<iframe_max;iframe++){

			ns = samples_struct.nsamples[iframe]; // number of samples for this scan
			f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp; // knee freq of the filter in terms of samples in order to compute fft
			f_lppix_Nk = samples_struct.fcut[iframe]*double(ns)/proc_param.fsamp; // noise PS threshold freq, in terms of samples

			std::vector<string> det_vect = bolo_list[iframe];
			long ndet = (long)det_vect.size();

			if(ndet!=samples_struct.ndet[iframe]){ // check here to avoid problems between saneInv and sanePic
				if(rank==0){
					cout << "Error. The number of detector in noisePower Spectra file must be egal to input bolofile number\n";
					cout << "Did you forgot to run saneInv ??? Exiting..." << endl;
				}
#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CONFIG;
			}

			// if there is correlation between detectors
			if (proc_param.CORRon){
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
				// f_lppix = filter freq in term of sample
				// ns = number of sample for this scan
				// iframe = scan number : 0=> ntotscan if non-MPI

				pb=write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det_vect,ndet,indpix,indpsrc,NAXIS1, NAXIS2,npix,
						npixsrc,addnpix,f_lppix,ns,	iframe,para_bolo_indice, para_bolo_size, name_rank);

				if(pb>0){
					cout << "Problem in write_ftrProcesdata. Exiting ...\n";
#ifdef USE_MPI
					MPI_Abort(MPI_COMM_WORLD, 1);
#endif
					return EX_SOFTWARE;
				}

#ifdef DEBUG
				time ( &rawtime ); // print processing time for each proc in the log file
				timeinfo = localtime ( &rawtime );
				file_rank.open(name_rank.c_str(), ios::out | ios::app);
				file_rank << "rank " << rank << " a fini et attend a " << asctime (timeinfo) << " \n";
				file_rank.close();
#endif
#ifdef PARA_BOLO
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				/* do_PtNd parameters : */
				// PNd => npix (dimension), initialised to 0.0 : (At N-1 d)
				// prefixe = "fdata"; => write/read fdata files prefix
				// det = bolo names array + bolo number
				// f_lppix_Nk = freq threshold noise (in term of samples)
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

				pb+=do_PtNd(samples_struct, PNd,dir.tmp_dir,prefixe,det_vect,ndet,f_lppix_Nk,
						proc_param.fsamp,ns,para_bolo_indice,para_bolo_size,indpix,NAXIS1, NAXIS2,npix,iframe,samples_struct.fitsvect[iframe], NULL, NULL, name_rank);
				// Returns Pnd = (At N-1 d), Mp and hits

				if(pb>0){
					cout << "Problem after do_PtNd. Exiting...\n";
#ifdef USE_MPI
					MPI_Abort(MPI_COMM_WORLD, 1);
#endif
					return -1;
				}


			} else { // No correlation case


				do_PtNd_nocorr(PNd, dir.tmp_dir,proc_param,pos_param,samples_struct,
						det_vect,ndet,f_lppix,f_lppix_Nk,addnpix,
						ns,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,iframe,NULL,rank,size);
				// fillgaps + butterworth filter + fourier transform and PNd generation

			}

		} // end of iframe loop

	} // end of else (restore = 0)

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if(struct_sanePic.restore>0){
		if(rank==0)
			for(long ii=0; ii<npix;ii++)
				PNdtot[ii]=PNd[ii];
	}else
		MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
	PNdtot=PNd;
#endif


	if (struct_sanePic.save_data > 0) {
		if (rank == 0) {
			struct checksum chk_t;
			/* Compute Checsum for crash recovery ! */
			//			compute_checksum(argv[indice_argv], dir.tmp_dir, npix,
			//					indpix, indpsrc, indpsrc_size, chk_t);
			compute_checksum(dir, pos_param, proc_param, saneInv_struct, structPS, struct_sanePic, samples_struct, npix,
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
		printf("\nStarting Conjugate Gradient Descent... \n\n");

	//////////////////////////////////// Computing of sanePic starts here
	string testfile; // log file to follow evolution of both criteria
	ostringstream temp_stream; // fits files filename string stream

	// inititialisation of the Conjugate gradient with preconditioner
	// see (for a complete description of the following variables) : http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
	double *PtNPmatS= NULL, *PtNPmatStot = NULL, *r, *q, *qtot = NULL, *d, *s; // =NULL to avoid warnings
	double *Mp,	*Mptot = NULL; // =NULL to avoid warnings
	// Mp = M in the paper = preconditioner


	double var0 = 0.0, var_n = 0.0, delta0 = 0.0, delta_n = 0.0, alpha = 0.0; // conjugate gradient convergence criteria
	double delta_o, rtq, beta; // conjugate gradient needed parameters (see paper for a complete description)

	long mi;
	double *map1d = NULL; // buffer used to store maps before exporting to fits

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

	if (pos_param.projgaps || !flagon) {
		npixeff = npix;
	} else {
		npixeff = npix - 1;
	}


	if ((struct_sanePic.restore > 0)) {
		if (rank == 0){
			cout << "loading idupl\n";
			load_idupl(dir.tmp_dir, dir.output_dir, idupl);

			// idupl compatibility
			if((idupl>0) && !pos_param.flgdupl){
				cout << "Error. idupl cannot be >0 if flgdupl is False ! Exiting...\n";
#ifdef USE_MPI
				MPI_Finalize();
#endif
				return EX_CONFIG;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if ((struct_sanePic.restore > 0))
		MPI_Bcast(&idupl,1,MPI_INT,0,MPI_COMM_WORLD);
	if(rank==0) { // malloc only for first processor that reduces the data
		//		PtNPmatStot = new double[npix];
		Mptot = new double[npix];

		qtot = new double[npix];
		fill(qtot,qtot+npix,0.0);
		//		fill(PtNPmatStot,PtNPmatStot+npix,0.0);
		fill(Mptot,Mptot+npix,0.0);
	}
#endif


	// in case flagged pixels are put in a duplicated map
	while(idupl <= pos_param.flgdupl){

		if (struct_sanePic.save_data > 0)
			if(rank==0)
				write_PNd(PNdtot,npix,dir.tmp_dir, "PNdCorr.bi");


		if((struct_sanePic.restore == 0)){

#ifdef USE_MPI
			if(rank==0) { // malloc only for first processor that reduces the data
				PtNPmatStot = new double[npix];
				fill(PtNPmatStot,PtNPmatStot+npix,0.0);
			}
#endif

			PtNPmatS = new double[npix];

			fill(PtNPmatS, PtNPmatS + npix, 0.0);
			fill(Mp, Mp + npix, 0.0);
			fill(r, r + npix, 0.0);
			fill(d, d + npix, 0.0);
			fill(s, s + npix, 0.0);


			for (long iframe = iframe_min; iframe < iframe_max; iframe++) {

				ns = samples_struct.nsamples[iframe];
				f_lppix_Nk = samples_struct.fcut[iframe] * double(ns) / proc_param.fsamp;

				std::vector<string> det_vect = bolo_list[iframe];
				long ndet = (long)det_vect.size();

				// preconditioner computation : Mp
				if (proc_param.CORRon) {


					write_tfAS(samples_struct, S, det_vect, ndet, indpix, NAXIS1, NAXIS2, npix,
							pos_param.flgdupl, dir.tmp_dir, ns,
							samples_struct.fitsvect[iframe], para_bolo_indice, para_bolo_size);

#ifdef DEBUG
					time ( &rawtime );
					timeinfo = localtime ( &rawtime );
					file_rank.open(name_rank.c_str(), ios::out | ios::app);
					file_rank << "rank " << rank << " a fini write_tfAS et attend a " << asctime (timeinfo) << " \n";
					file_rank.close();
#endif
#ifdef PARA_BOLO
					MPI_Barrier(MPI_COMM_WORLD);
#endif

					do_PtNd(samples_struct, PtNPmatS, dir.tmp_dir,
							"fPs_", det_vect, ndet, f_lppix_Nk, proc_param.fsamp, ns, para_bolo_indice,
							para_bolo_size, indpix, NAXIS1, NAXIS2, npix, iframe,
							samples_struct.fitsvect[iframe], Mp, NULL, name_rank);

				} else {

					do_PtNPS_nocorr(samples_struct, S, samples_struct.noisevect, dir, det_vect, ndet,
							f_lppix_Nk, proc_param.fsamp, pos_param.flgdupl, ns,
							indpix, NAXIS1, NAXIS2, npix, iframe,
							samples_struct.fitsvect[iframe], PtNPmatS, Mp, NULL,
							rank, size);
				}

			} // end of iframe loop

#ifdef USE_MPI
			if(rank==0){
				fill(Mptot,Mptot+npix,0.0);
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
			if(rank==0)
				delete [] PtNPmatStot;
#endif

		}else{ //restore

			if (rank == 0){
				cout << "loading session\n";

				// fill S, d, r, indpix, npixeff, var_n, delta_n and iter with previously saved on disk values
				load_from_disk(dir.tmp_dir, dir.output_dir, S, d, r, npixeff, var0, var_n, delta0,
						delta_n, iter, Mp);
			}
#ifdef USE_MPI
			MPI_Bcast(&iter,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(S,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);

			//copy Mp to Mptot and PtNPmatS to PtNPmatStot
			if(rank==0)
				for(long ii = 0; ii< npix; ii++){
					Mptot[ii] = Mp[ii];
				}
#else
			Mptot = Mp;
#endif

			cout << iter << " " << npixeff << " " << var0 << " " << var_n << " " << delta0 << " "
					<< delta_n << endl;

			struct_sanePic.restore=0; // set to 0 because we don't want to load again on idupl=1 loop !


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
		while ((iter < struct_sanePic.itermax) && (((var_n / var0 > 1e-10) && (idupl
				|| !pos_param.flgdupl)) || (!idupl && (var_n / var0 > 1e-6)))) {

			fill(q, q + npixeff, 0.0); // q <= A*d

			for (long iframe = iframe_min; iframe < iframe_max; iframe++) {

				ns = samples_struct.nsamples[iframe];
				f_lppix_Nk = samples_struct.fcut[iframe] * double(ns) / proc_param.fsamp;

				std::vector<string> det_vect = bolo_list[iframe];
				long ndet = (long)det_vect.size();

				if (proc_param.CORRon) {

					write_tfAS(samples_struct, d, det_vect, ndet, indpix, NAXIS1, NAXIS2, npix,
							pos_param.flgdupl, dir.tmp_dir, ns,
							samples_struct.fitsvect[iframe], para_bolo_indice, para_bolo_size);

#ifdef DEBUG
					time ( &rawtime );
					timeinfo = localtime ( &rawtime );
					file_rank.open(name_rank.c_str(), ios::out | ios::app);
					file_rank << "rank " << rank << " a fini write_tfAS et attend a " << asctime (timeinfo) << " \n";
					file_rank.close();
#endif
#ifdef PARA_BOLO
					MPI_Barrier(MPI_COMM_WORLD);
#endif

					do_PtNd(samples_struct, q, dir.tmp_dir, "fPs_",
							det_vect, ndet, f_lppix_Nk, proc_param.fsamp, ns, para_bolo_indice, para_bolo_size,
							indpix, NAXIS1, NAXIS2, npix, iframe,
							samples_struct.fitsvect[iframe], NULL, NULL,
							name_rank);

				} else {

					do_PtNPS_nocorr(samples_struct, d, samples_struct.noisevect, dir, det_vect, ndet,
							f_lppix_Nk, proc_param.fsamp, pos_param.flgdupl,
							ns, indpix, NAXIS1, NAXIS2, npix, iframe,
							samples_struct.fitsvect[iframe], q, NULL, NULL,
							rank, size);
				}
			} // end of iframe loop


#ifdef USE_MPI
			if(rank==0)
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
#ifdef USE_MPI
				if(rank==0) { // malloc only for first processor that reduces the data
					PtNPmatStot = new double[npix];
					fill(PtNPmatStot,PtNPmatStot+npix,0.0);
				}
#endif

				PtNPmatS = new double[npix];
				fill(PtNPmatS, PtNPmatS + npix, 0.0);

				for (long iframe = iframe_min; iframe < iframe_max; iframe++) {
					ns = samples_struct.nsamples[iframe];
					f_lppix_Nk = samples_struct.fcut[iframe] * double(ns) / proc_param.fsamp;

					std::vector<string> det_vect = bolo_list[iframe];
					long ndet = (long)det_vect.size();

					if (proc_param.CORRon) {

						write_tfAS(samples_struct, S, det_vect, ndet, indpix, NAXIS1, NAXIS2, npix,
								pos_param.flgdupl, dir.tmp_dir, ns,
								samples_struct.fitsvect[iframe], para_bolo_indice, para_bolo_size);

#ifdef DEBUG
						time ( &rawtime );
						timeinfo = localtime ( &rawtime );
						file_rank.open(name_rank.c_str(), ios::out | ios::app);
						file_rank << "rank " << rank << " a fini write_tfAS et attend a " << asctime (timeinfo) << " \n";
						file_rank.close();
#endif
#ifdef PARA_BOLO
						MPI_Barrier(MPI_COMM_WORLD);
#endif

						do_PtNd(samples_struct, PtNPmatS, dir.tmp_dir, "fPs_", det_vect, ndet, f_lppix_Nk,
								proc_param.fsamp, ns, para_bolo_indice, para_bolo_size, indpix,
								NAXIS1, NAXIS2, npix, iframe,
								samples_struct.fitsvect[iframe], NULL, NULL,
								name_rank);

					} else {

						do_PtNPS_nocorr(samples_struct, S, samples_struct.noisevect, dir,
								det_vect, ndet, f_lppix_Nk, proc_param.fsamp,
								pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2,
								npix, iframe,
								samples_struct.fitsvect[iframe], PtNPmatS,
								NULL, NULL, rank, size);
					}

				} // end of iframe loop


#ifdef USE_MPI
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
				if(rank==0)
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


				delta_o = delta_n; // delta_0 <= delta_new

				delta_n = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					delta_n += r[ii] * s[ii]; // delta_new = rT * s

				var_n = 0.0;
				for (long ii = 0; ii < npixeff; ii++)
					var_n += r[ii] * r[ii];

				beta = delta_n / delta_o; // beta = delta_new / delta_0
				for (long ii = 0; ii < npixeff; ii++)
					d[ii] = s[ii] + beta * d[ii]; // d = s + beta * d

				// saving iterated maps
				if (struct_sanePic.iterw && (iter % struct_sanePic.iterw) == 0) {

					if ((struct_sanePic.save_data > 0) && ((iter > 0) || (idupl > 0))){
						write_disk(dir.tmp_dir, d, r, S, npixeff, var0, var_n, delta0, delta_n, iter, idupl, Mptot);
					}

					// Every iterw iteration compute the map and save it
					map1d = new double[NAXIS1 * NAXIS2];

					for (long ii = 0; ii < NAXIS1; ii++) {
						for (long jj = 0; jj < NAXIS2; jj++) {
							mi = jj * NAXIS1 + ii;
							if (indpix[mi] >= 0) {
								map1d[mi] = S[indpix[mi]];
							} else {
								map1d[mi] = NAN;
							}
						}
					}

					temp_stream << dir.output_dir + struct_sanePic.map_prefix + "_" << iter	<< (idupl>0 ? (std::string)"b" : (std::string)"a") + ".fits";
					fname = temp_stream.str();
					temp_stream.str("");
					if (write_fits_wcs("!" + fname, wcs, NAXIS1, NAXIS2, 'd',
							(void *) map1d, (char *) "Image", 0,key,datatype,val,com)) {
						cerr << "Error Writing map : EXITING ... \n";
#ifdef USE_MPI
						MPI_Abort(MPI_COMM_WORLD, 1);
#endif
						return EX_CANTCREAT;
					}

					if (pos_param.flgdupl) {
						for (long ii = 0; ii < NAXIS1; ii++) {
							for (long jj = 0; jj < NAXIS2; jj++) {
								mi = jj * NAXIS1 + ii;
								if (indpix[mi] >= 0) { // pixel observed in first map
									if (indpix[mi + NAXIS1 * NAXIS2] >= 0)
										map1d[mi] = S[indpix[mi + NAXIS1
										                     * NAXIS2]]; //-finalmap[ii][jj];
									else
										map1d[mi] = INFINITY;
								} else {
									map1d[mi] = NAN;
								}
							}
						}

						if (write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd',
								(void *) map1d, (char *) "Flagged Data", 1,key,datatype,val,com)) {
							cerr << "Error Writing map : EXITING ... \n";
#ifdef USE_MPI
							MPI_Abort(MPI_COMM_WORLD, 1);
#endif
							return EX_CANTCREAT;
						}
					}
					if (addnpix) {
						// initialize the container
						for (long jj = 0; jj < NAXIS2; jj++) {
							for (long ii = 0; ii < NAXIS1; ii++) {
								mi = jj * NAXIS1 + ii;
								map1d[mi] = 0.0;
							}
						}
						// loop thru frame to coadd all pixels
						for (long ii = 0; ii < NAXIS1; ii++) {
							for (long jj = 0; jj < NAXIS2; jj++) {
								mi = jj * NAXIS1 + ii;
								double b = 0.;
								for (long iframe = 0; iframe
								< samples_struct.ntotscan; iframe++) {
									//TODO : Better : Make a weithted mean with Mptot
									long long ll = factdupl * NAXIS1 * NAXIS2
											+ iframe * npixsrc + indpsrc[mi];
									if ((indpsrc[mi] != -1) && (indpix[ll]
									                                   != -1)) {
										map1d[mi] += S[indpix[ll]] / gsl_pow_2(
												Mptot[indpix[ll]]);
										b += 1. / gsl_pow_2(Mptot[indpix[ll]]);
									}
								}
								if (b >= 1)
									map1d[mi] /= b;
							}
						}
						// replace the non observed pixels by NAN
						for (long ii = 0; ii < NAXIS1; ii++) {
							for (long jj = 0; jj < NAXIS2; jj++) {
								mi = jj * NAXIS1 + ii;
								if (map1d[mi] == 0.0)
									map1d[mi] = NAN;
							}
						}
						if (write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd',
								(void *) map1d, "CCR Image", 1,key,datatype,val,com)) {
							cerr << "Error Writing map : EXITING ... \n";
#ifdef USE_MPI
							MPI_Abort(MPI_COMM_WORLD, 1);
#endif
							return EX_CANTCREAT;
						}

						for (long ii = 0; ii < NAXIS1; ii++) {
							for (long jj = 0; jj < NAXIS2; jj++) {
								mi = jj * NAXIS1 + ii;
								for (long iframe = 0; iframe
								< samples_struct.ntotscan; iframe++) {

									long long ll = factdupl * NAXIS1 * NAXIS2
											+ iframe * npixsrc + indpsrc[mi];
									if ((indpsrc[mi] != -1) && (indpix[ll]
									                                   != -1))
										map1d[mi] += 1. / gsl_pow_2(
												Mptot[indpix[ll]]);
								}
								if (map1d[mi] != 0)
									map1d[mi] = 1. / sqrt(map1d[mi]);
							}
						}
						if (write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd',
								(void *) map1d, (char *) "CCR Error", 1,key,datatype,val,com)) {
							cerr << "Error Writing map : EXITING ... \n";
#ifdef USE_MPI
							MPI_Abort(MPI_COMM_WORLD, 1);
#endif
							return EX_CANTCREAT;
						}
					}
					if(write_fits_hitory2(fname, NAXIS1, NAXIS2, dir, proc_param, pos_param,
							samples_struct.fcut, samples_struct, structPS, struct_sanePic, saneInv_struct)) // write sanePre parameters in naive Map fits file header
						cerr << "WARNING ! No history will be included in the file : " << fname << endl;

					delete[] map1d;
				} // end of saving iterated maps

				char mytime[25];

				time_t rawtime;
				time ( &rawtime );
				strncpy(mytime, ctime(&rawtime),24); // remove the newline character
				mytime[24] = '\0';
				temp_stream.str("");
				temp_stream << "[" << mytime << "] ";
				temp_stream << "iter = "      << setw(4) << iter;
				temp_stream << ", crit = "    << setiosflags(ios::scientific) << setprecision (2) << var_n / var0;
				temp_stream << ", crit2 = "   << setiosflags(ios::scientific) << setprecision (2) << delta_n / delta0;

				// Output to screen ...
				cout << temp_stream.str() << "\r" << flush;
				cout << endl; // temp : to test for dirfiles !

				// ... and to a logfile
				ofstream logfile;
				string filename = dir.output_dir + "ConvFile.txt";
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

		if ((pos_param.projgaps || (pos_param.flgdupl)) && !idupl) {

			fill(PNd, PNd + npix, 0.0); // correct : has to be reset to 0 !

			for (long iframe = iframe_min; iframe < iframe_max; iframe++) {

				ns = samples_struct.nsamples[iframe];
				f_lppix = proc_param.f_lp * double(ns) / proc_param.fsamp;
				f_lppix_Nk = samples_struct.fcut[iframe] * double(ns) / proc_param.fsamp;

				std::vector<string> det_vect = bolo_list[iframe];
				long ndet = (long)det_vect.size();

				if (proc_param.CORRon) {

					write_ftrProcesdata(S, proc_param, samples_struct,
							pos_param, dir.tmp_dir, det_vect, ndet, indpix, indpsrc,
							NAXIS1, NAXIS2, npix, npixsrc, addnpix, f_lppix,
							ns, iframe, para_bolo_indice, para_bolo_size, name_rank);

#ifdef DEBUG
					time ( &rawtime );
					timeinfo = localtime ( &rawtime );
					file_rank.open(name_rank.c_str(), ios::out | ios::app);
					file_rank << "rank " << rank << " a fini et attend a " << asctime (timeinfo) << " \n";
					file_rank.close();
#endif
#ifdef PARA_BOLO
					MPI_Barrier(MPI_COMM_WORLD);
#endif

					do_PtNd(samples_struct, PNd, dir.tmp_dir,
							"fdata_", det_vect, ndet, f_lppix_Nk, proc_param.fsamp, ns,
							para_bolo_indice, para_bolo_size, indpix, NAXIS1, NAXIS2, npix, iframe,
							samples_struct.fitsvect[iframe], NULL, NULL,
							name_rank);

				} else {

					do_PtNd_nocorr(PNd, dir.tmp_dir, proc_param, pos_param,
							samples_struct, det_vect, ndet, f_lppix, f_lppix_Nk, addnpix,
							ns, indpix, indpsrc, NAXIS1, NAXIS2, npix, npixsrc,
							iframe, S, rank, size);
				}

			} // end of iframe loop

#ifdef USE_MPI
			if(rank==0)
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



	//******************************  write final map in file ********************************

	if (rank == 0) {
		cout << "\ndone.\n";
		//		printf(" after CC INVERSION %lld\n", npix * (npix + 1) / 2);

#ifdef DEBUG
		FILE * fp;
		fp = fopen("output_signal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);
		fclose(fp);
#endif

		if(write_maps_to_disk(S, NAXIS1, NAXIS2, npix, dir, indpix, indpsrc,
				Mptot, addnpix, npixsrc, factdupl, samples_struct.ntotscan,
				proc_param, pos_param, samples_struct, samples_struct.fcut, wcs,
				pos_param.maskfile, structPS, struct_sanePic, saneInv_struct, key, datatype, val, com, bolo_list))
			cout << "Error in write_maps_to_disk. Exiting ...\n"; // don't return here ! let the code do the dealloc and return

	}// end of rank==0


	// clean up of conjugate gradient variables
	delete[] r;
	delete[] q;
	delete[] d;
	delete[] Mp;
	delete[] s;
	delete[] PNd;

#ifdef USE_MPI
	delete [] qtot;
	delete [] Mptot;
	delete [] PNdtot;

	delete [] ini_v.fitsvect;
	delete [] ini_v.noisevect;
	delete [] ini_v.bolovect;
#endif

	//******************************************************************//
	//******************************************************************//
	//*********************** End of sanePic ***************************//
	//******************************************************************//
	//******************************************************************//


#ifdef DEBUG
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	file_rank.open(name_rank.c_str(), ios::out | ios::app);
	if(!file_rank.is_open()) {
		cerr << "File [" << file_rank << "] Invalid." << endl;
		return EX_CANTCREAT;
	}

	file_rank << "[ " << rank << " ] Finish Time : " << asctime (timeinfo) << endl;
	file_rank.close();
#endif

	// clean up
	delete[] S;
	delete[] indpsrc;
	delete[] indpix;

	wcsvfree(&nwcs, &wcs);

	fftw_cleanup();

	if (gd_close(samples_struct.dirfile_pointer))
		cout << "error closing dirfile : " << filedir << endl;

#ifdef USE_MPI
	MPI_Finalize();
#endif

	if(rank==0)
		cout << "\nEND OF SANEPIC" << endl;

	return EXIT_SUCCESS;
}

