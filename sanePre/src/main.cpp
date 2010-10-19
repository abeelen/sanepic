#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>


#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}

#ifdef PARA_BOLO
#define USE_MPI
#endif


#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;


//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//


/*! \mainpage Sanepic for SPIRE
 *
 * \section intro_sec Matthieu HUSSON & Alexandre Beelen
 *
 * Sanepic for SPIRE Manual
 */



int main(int argc, char *argv[])
/*! Sanepic preprocess main function */
{

	int size;/*!< number of processors */
	int rank;

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(rank==0)
		printf("\nsanepic_preprocess\n");

#else
	size = 1;
	rank = 0;
	printf("\n sanepic_preprocess\n");
	cout << "Mpi is not used for this step" << endl;
#endif


	struct param_process proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct;  /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_positions pos_param; /*! A structure that contains user options about map projection and properties */
	struct common dir; /*! structure that contains output input temp directories */
	//	struct detectors det; /*! A structure that contains everything about the detectors names and number */
	std::vector<detectors> detector_tab;

	// default parameters
	int nwcs=1; /// number of wcs that will be used
	long iframe_min=0, iframe_max=0; /*!  min and max number of frame (used with mpi) */
	int factdupl=1; /*! map duplication factor */
	int flagon = 0; /*! if a data is flagged */




	long NAXIS1, NAXIS2; // Map size (pixels)
	long long npix; /*! npix = number of filled pixels*/
	long long npixsrc = 0; /*! number of pixels contained in box constraint removal */
	long long addnpix=0; /* number of pixels to add to the final map */
	long long ind_size; /* indpix readed size */

	//internal data params
	long ns; /*! number of samples for this scan, first frame number of this scan*/
	double f_lppix, f_lppix_Nk; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples*/


	//fftw_complex *fdata_buffer; /*! buffer used to store all the fdata arrays instead of writing on disk */


	double *PNd, *PNdtot=NULL; /*!  projected noised data, and global Pnd for mpi utilization */
	double *PNdNaiv, *PNdtotNaiv=NULL; /*!  projected noised data, and global Pnd for mpi utilization */
	double *Mp, *Mptot=NULL;
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices*/
	long *hits, *hitstot=NULL; /*! naivmap parameters : hits count */
	long *hitsNaiv, *hitstotNaiv=NULL; /*! naivmap parameters : hits count */

	string field; /*! actual boloname in the bolo loop*/
	string fname; // parallel scheme file name
	string prefixe; /*! prefix used for temporary name file creation*/


	/* Parser inputs */
	std::vector<double> fcut; /* noise cutting frequency */

	// Processing time estimation
	time_t t2, t3;// t4, t5, dt;

	int parsed=0; // parser error status


	if (argc<2)/* not enough argument */
		parsed=1;
	else {
		// Parse ini file

		// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
//		double fcut_sanePS=0.0;
		//		string MixMatfile, signame;
		//		long ncomp=1;
		//		int iterw=10;
		//		int save_data, restore;
		struct PS structPS;
		struct sanePic struct_sanePic;
		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], dir, detector_tab, samples_struct, pos_param, proc_param, fcut,
				structPS, struct_sanePic, rank, size);
	}



	if (rank==0)
		switch (parsed){/* error during parsing phase */

		case 1: printf("Please run %s using a *.ini file\n",argv[0]);
		break;

		case 2 : printf("Wrong program options or argument. Exiting !\n");
		break;

		case 3 : printf("Exiting...\n");
		break;

		default :;
		}

	// in case there is a parsing error or the dirfile format file was not created correctly
	if ((parsed>0)||(!compute_dirfile_format_fdata(dir.tmp_dir, samples_struct, detector_tab, rank))){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}


#ifdef DEBUG
	std::ostringstream oss;
	string name_rank;
	oss << dir.output_dir + "debug_sanePre_" << rank << ".txt"; // generate a log file
	name_rank = oss.str();
#else
	string name_rank = dir.output_dir + "debug_sanePre.txt";

#endif

	ofstream file_rank;
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime ); // allocate a time structure to print in the log file

	file_rank.open(name_rank.c_str(), ios::out | ios::trunc); // creating log file for debug mode
	if(!file_rank.is_open()){
		cerr << "File [" << file_rank << "] Invalid." << endl;
		return -1;
	}
	file_rank << "Opening file for writing debug at " << asctime (timeinfo)  << endl;
	file_rank.close();

	// processing begins here
	t2=time(NULL);

	// allocate memory for ...
	samples_struct.fits_table  = new string[samples_struct.ntotscan]; // ...input fits file list
	samples_struct.index_table = new int[samples_struct.ntotscan]; // ... user processor index
	samples_struct.noise_table = new string[samples_struct.ntotscan]; // ... input covariance matrices filenames


	if (pos_param.flgdupl) factdupl = 2;// default 0 : if flagged data are put in a duplicated map

	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs,&NAXIS1, &NAXIS2); // read wcs file header generated by sanePos

	if(rank==0)
		cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;


	long long test_size; // used to verify indpsrc size
	if(read_indpsrc( test_size, npixsrc, indpsrc,  dir.tmp_dir)){ // read indpsrc (mask index) file generated by sanePos
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return(EXIT_FAILURE);
	}
	if(test_size != NAXIS1*NAXIS2){
		if(rank==0)
			cout << "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos" << endl;
		exit(0);
	}
	// each frame contains npixsrc pixels with index indsprc[] for which
	// crossing constraint are removed
	// thus
	// addnpix = number of pix to add in pixon
	//         = number of scans * number of pix in box crossing constraint removal
	addnpix = samples_struct.ntotscan*npixsrc;


	//read projection vector from a file
	if(read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon)){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return(EXIT_FAILURE);
	}

	// Check indpix readed size = expected size
	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		if(rank==0){
			cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
			cout << ind_size << " != "  << (factdupl*NAXIS1*NAXIS2+2 + addnpix) << " " << factdupl << " " << addnpix << endl;
		}
		exit(0);
	}

	//TODO: check the noise matrix file at the beginning of each frame/start of the program...


#if defined(USE_MPI) && ! defined(PARA_BOLO)

	ofstream file;

	if(samples_struct.scans_index.size()==0){ // if user has not given a processor index in file_list

		int test=0;
		fname = dir.output_dir + parallel_scheme_filename; // get parallel scheme file name
		cout << fname << endl;

		// and spread scans between processors
		test = define_parallelization_scheme(rank,fname,dir.input_dir,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){ // define_parallelization did not worked : exit program
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
	}else{
		int test=0; // User has given a processor index in file_list
		// Verify its validity
		test = verify_parallelization_scheme(rank,dir.output_dir,samples_struct, size, iframe_min, iframe_max);


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&test,1,MPI_INT,0,MPI_COMM_WORLD); // all proc have to know if this operation worked fine or not

		if(test>0){ // if not, all exit
			MPI_Finalize();
			exit(0);

		}

	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // ifram_min=iframe_max => This processor will not do anything
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//	for(long ii=0;ii<size;ii++){ // print processors indexes
	//		if(rank==ii)
	//			cout << "[ " << rank << " ]. iframemin : " << iframe_min << " iframemax : " << iframe_max << endl;
	//		else
	//			MPI_Barrier(MPI_COMM_WORLD);
	//	}

#else
#if defined(USE_MPI) && defined(PARA_BOLO)

	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

#else
	fname = dir.output_dir + parallel_scheme_filename; // get parallel file name
	int test=0;

	// Check its validity
	test=check_ParallelizationScheme(fname,dir.input_dir,samples_struct,size);
	if (test==-1){
		if(rank==0)
			cerr << "erreur dans check_parallelizationScheme non-MPI " << endl;
		exit(0);
	}

	for(int ii = 0; ii< samples_struct.ntotscan;ii++){ // add input directory to input fits file names
		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	}

#endif

	iframe_min = 0; // single processor, will compute all the scans
	iframe_max = samples_struct.ntotscan;


#endif

	// (At N-1 D) memory allocation
	PNd = new double[npix];
	PNdNaiv= new double[npix];
	//
	hits=new long[npix];
	hitsNaiv=new long[npix];
	Mp = new double[npix];

	// initialize to 0.0
	fill(PNd,PNd+npix,0.0);
	fill(hits,hits+npix,0);
	fill(Mp,Mp+npix,0.0);
	fill(PNdNaiv,PNdNaiv+npix,0.0);
	fill(hitsNaiv,hitsNaiv+npix,0);


#ifdef USE_MPI

	if(rank==0){	// global (At N-1 D) malloc for mpi
		PNdtot = new double[npix];
		hitstot=new long[npix];
		Mptot = new double[npix];
		PNdtotNaiv = new double[npix];
		hitstotNaiv=new long[npix];
	}

#endif


	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//
	if(rank==0)
		printf("\nPre-processing of the data\n");

	prefixe = "fdata_"; // Fourier transform of the data file prefixe

	// loop over the scans
	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

#ifdef LARGE_MEMORY
		fftw_complex  *fdata_buffer; // fdata are saved in buffers instead of written to disk
		//if(rank==0)
		fftw_complex *fdata_buffer_tot=NULL;
#endif

		ns = samples_struct.nsamples[iframe]; // number of samples for this scan
		f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp; // knee freq of the filter in terms of samples in order to compute fft
		f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp; // noise PS threshold freq, in terms of samples

		struct detectors det = detector_tab[iframe];
		//		if(iframe_min!=iframe_max)
		//			printf("[%2.2i] iframe : %ld/%ld\n",rank,iframe+1,iframe_max);


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

#ifdef LARGE_MEMORY
			// A fdata buffer will be used to avoid binary writing
			fdata_buffer = new fftw_complex[det.ndet*(ns/2+1)];


			if (rank==0){ // allocate buffer + fill with 0.0
				fdata_buffer_tot = new fftw_complex[det.ndet*(ns/2+1)];
				for (long ii=0;ii<det.ndet*(ns/2+1);ii++){
					fdata_buffer_tot[ii][0] = 0.0;
					fdata_buffer_tot[ii][1] = 0.0;
				}


			}

			pb=write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,rank,size,name_rank,fdata_buffer);
			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)
#else

#if defined(USE_MPI) && ! defined(PARA_BOLO)
			pb=write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,0,1, name_rank);
			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)
#else
			pb=write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,rank,size, name_rank);
			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)
#endif
#endif

			if(pb>0){
				cout << "Problem in write_ftrProcesdata. Exiting ...\n";
#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return -1;
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
#ifdef LARGE_MEMORY
			// Rank 0 collects the fourier transform buffer...
			MPI_Reduce(fdata_buffer,fdata_buffer_tot,(ns/2+1)*2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Bcast(fdata_buffer,(ns/2+1)*2,MPI_DOUBLE,0,MPI_COMM_WORLD); // ...and broadcast it to the other procs
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
#ifdef LARGE_MEMORY
			pb+=do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank,fdata_buffer);
			// Returns Pnd = (At N-1 d), Mp and hits
#else

#if defined(USE_MPI) && ! defined(PARA_BOLO)
			pb+=do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,0,1,indpix,NAXIS1, NAXIS2,npix,iframe,samples_struct.fits_table[iframe],Mp,hits, name_rank);

			pb+=do_PtNd_Naiv(PNdNaiv, dir.tmp_dir, samples_struct.fits_table, det,proc_param.poly_order, proc_param.napod, f_lppix, ns, 0, 1, indpix, iframe, hitsNaiv);
			// Returns Pnd = (At N-1 d), Mp and hits
#else
			pb+=do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,samples_struct.fits_table[iframe],Mp,hits, name_rank);
			// Returns Pnd = (At N-1 d), Mp and hits

			pb+=do_PtNd_Naiv(PNdNaiv, dir.tmp_dir, samples_struct.fits_table, det, proc_param.poly_order, proc_param.napod, f_lppix, ns, rank, size, indpix, iframe, hitsNaiv);

#endif
#endif


			if(pb>0){
				cout << "Problem after do_PtNd. Exiting...\n";
#ifdef USE_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
#endif
				return -1;
			}

#ifdef LARGE_MEMORY
			delete [] fdata_buffer;
			if(rank==0)
				delete [] fdata_buffer_tot;
#endif
		} else { // No correlation case


			do_PtNd_nocorr(PNd, dir.tmp_dir,proc_param,pos_param,samples_struct,
					det,f_lppix,f_lppix_Nk,addnpix,
					ns,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,iframe,NULL,rank,size);
			// fillgaps + butterworth filter + fourier transform and PNd generation

		}


	} // end of iframe loop



#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(PNdNaiv,PNdtotNaiv,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(hitsNaiv,hitstotNaiv,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);

#else
	hitstot=hits;
	PNdtot=PNd;
	Mptot=Mp;
	PNdtotNaiv=PNdNaiv;
	hitstotNaiv=hitsNaiv;

#endif

	if(rank==0)
		printf("\nEnd of Pre-Processing\n");

	if (rank == 0){
		// write (At N-1 d) in a file
		if(write_PNd(PNdtot,npix,dir.tmp_dir)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}


		cout << "\nNaive step :" << endl;
		string fnaivname;
		double *map1d;
		long long mi;
		map1d = new double[NAXIS1*NAXIS2];

		//TODO:   addnpix pixel
		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					if(PNdtot[indpix[mi]]==NAN) // replace NAN by 0
						map1d[mi] = 0;
					else
						map1d[mi] = PNdtot[indpix[mi]]/Mptot[indpix[mi]];
				} else {
					map1d[mi] = NAN; // replace map[hits[mi] = 0] = NAN
				}
			}
		}

		fnaivname = '!' + dir.output_dir + "naivMap.fits";
		cout << "Output file : " << fnaivname << endl;
		if(write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Image",0)){ // create naive Map fits file with first de-noised map image
			cerr << "Error Writing map : EXITING ... \n";
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return 0;
		}
		// TODO: Save the map of flag data if needed

		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					//					if(hitsNaiv[indpix[mi]]>0)
					map1d[mi] = PNdtotNaiv[indpix[mi]]/(double)hitstotNaiv[indpix[mi]];
				} else {
					map1d[mi] = NAN;
				}
			}
		}

		fnaivname = dir.output_dir + "naivMap.fits";
		if(write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Ultra Naiv",1)){ // open naive Map fits file and fill ultra naive map image
			cerr << "Error Writing Ultra Naiv map ... \n";
		}

		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = hitstot[indpix[mi]];
				} else {
					map1d[mi] = 0;
				}
			}
		}

		if (addnpix){
			for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++){
				for (long jj=0; jj<NAXIS2; jj++) {
					for (long ii=0; ii<NAXIS1; ii++) {
						mi = jj*NAXIS1 + ii;
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
							map1d[mi] += hitstot[indpix[ll]];
					}
				}
			}
		}

		if(	write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Coverage",1)){ // open naive Map fits file and fill hit (or coverage) image
			cerr << "Error Writing coverage map  ... \n";
		}

		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = Mptot[indpix[mi]];
				} else {
					map1d[mi] = NAN;
				}
			}
		}


		if(write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Error",1)){ // open naive Map fits file and fill Noise error image
			cerr << "Error Writing Error map  ... \n";
		}

		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii;
				map1d[mi] = 0.0;
			}
		}

		//TODO : Treat the 2 pixels properly
		if (addnpix){
			for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
				for (long ii=0; ii<NAXIS1; ii++) {
					for (long jj=0; jj<NAXIS2; jj++) {
						mi = jj*NAXIS1 + ii;
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
							map1d[mi] += 1.0/Mptot[indpix[ll]];
					}
				}
			}

			if(write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Invnoisevaruncpix",1)){ //
				cerr << "Error Writing Uncorralated noise variance map : EXITING ... \n";
			}
		}


		if(write_fits_hitory(fnaivname, NAXIS1, NAXIS2, dir.dirfile, proc_param, pos_param , fcut, detector_tab[0], samples_struct)) // write sanePre parameters in naive Map fits file header
			cerr << "WARNING ! No history will be included in the file : " << fnaivname << endl;
		if (pos_param.maskfile != "")
			if(write_fits_mask(fnaivname, pos_param.maskfile)) // copy mask in naive map file
				cerr << "Warning ! The mask will not be included in naive map fits file ...\n";

		delete [] map1d;
		if(rank==0)
			printf("End of saneNaiv\n\n");

	}
	/* ---------------------------------------------------------------------------------------------*/

	//Get processing time
	t3=time(NULL);



	//debug : computation time
#ifdef USE_MPI
	if(iframe_min!=iframe_max)
		printf("[%2.2i] Time : %d sec\n",rank, (int)(t3-t2));


	if(rank==0){
		// clean up
		delete [] PNdtot;
		delete [] Mptot;
		delete [] hitstot;
		delete [] PNdtotNaiv;
		delete [] hitstotNaiv;
	}
#else
	cout << "Total Time : " << t3-t2 << " sec\n";
#endif

#ifdef DEBUG
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	file_rank.open(name_rank.c_str(), ios::out | ios::app);
	if(!file_rank.is_open()){
		cerr << "File [" << file_rank << "] Invalid." << endl;
		return -1;
	}

	file_rank << "[ " << rank << " ] Finish Time : " << asctime (timeinfo) << endl; // print total processing time in log file
	file_rank.close();
#endif



	// clean up
	delete [] PNd;
	delete [] PNdNaiv;
	delete [] Mp;
	delete [] hits;
	delete [] hitsNaiv;



	delete [] samples_struct.nsamples;
	delete [] indpix;
	delete [] indpsrc;

	delete [] samples_struct.noise_table;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;

	wcsvfree(&nwcs, &wcs); // clean WCS structure

	fftw_cleanup();

#ifdef USE_MPI
	MPI_Finalize();
#endif

	if(rank==0)
		printf("\nEnd of sanePre \n");

	return 0;
}

//******************************************************************//
//******************************************************************//
//**********************  End of init loop *************************//
//******************************************************************//
//******************************************************************//



