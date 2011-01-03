#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>
#include <sysexits.h>


#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "inputFileIO.h"


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}


#if defined(PARA_BOLO) || defined(PARA_FRAME)
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
#else
	size = 1;
	rank = 0;
#endif

	if(rank==0)
		cout << endl << "sanePre :  Pre Processing of the data" << endl;

	struct param_sanePre proc_param; /*! contains user options about preprocessing properties */
	struct samples samples_struct;  /*  everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! contains user options about map projection and properties */
	struct param_common dir; /*! contains output input temp directories */

	// default parameters
	//	int nwcs=1; /// number of wcs that will be used
	long iframe_min=0, iframe_max=0; /*!  min and max number of frame (used with mpi) */
	//	int factdupl=1; /*! map duplication factor */
	//	int flagon = 0; /*! if a data is flagged */


	//	long NAXIS1, NAXIS2; // Map size (pixels)
	//	long long npix; /*! npix = number of filled pixels*/
	//	long long npixsrc = 0; /*! number of pixels contained in box constraint removal */
	//	long long addnpix=0; /* number of pixels to add to the final map */
	//	long long ind_size; /* indpix readed size */

	//internal data params
	//	long ns; /*! number of samples for this scan, first frame number of this scan*/
	//	double f_lppix, f_lppix_Nk; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples*/
	//
	//	double *PNdNaiv, *PNdtotNaiv=NULL; /*!  projected noised data, and global Pnd for mpi utilization */
	////	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices*/
	//	long *hitsNaiv, *hitstotNaiv=NULL; /*! naivmap parameters : hits count */

	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;
	string parser_output = "";

	// Processing time estimation
	time_t t2, t3;

	int parsed=0; // parser error status


	if (argc<2)/* not enough argument */
		parsed=1;
	else {
		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, struct_sanePic);

	}



	if (rank==0){
		// print parser warning and/or errors
		cout << endl << parser_output << endl;
		switch (parsed){/* error during parsing phase */

		case 1: printf("Please run %s using a *.ini file\n",argv[0]);
		break;

		case 2 : printf("Wrong program options or argument. Exiting !\n");
		break;

		case 3 : printf("Exiting...\n");
		break;

		default :;
		}
	}

	// in case there is a parsing error or the dirfile format file was not created correctly
	if (parsed>0){
#ifdef USE_MPI
		//		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
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
	file_rank << "Opening file for writing debug at " << asctime (timeinfo);
	file_rank.close();

	// processing begins here
	t2=time(NULL);

#ifdef PARA_FRAME

	if(configure_PARA_FRAME_samples_struct(dir.output_dir, samples_struct, rank, size, iframe_min, iframe_max)){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(EX_IOERR);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // ifram_min=iframe_max => This processor will not do anything
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
	}
#else

	iframe_min = 0; // single processor, will compute all the scans
	iframe_max = samples_struct.ntotscan;

#endif

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic);

		compute_dirfile_format_file(dir.tmp_dir, samples_struct);
	}

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(write_data_flag_to_dirfile(dir, samples_struct, iframe_min, iframe_max)){
		cout << "error write data !! Exiting ...\n";
#ifdef USE_MPI
		//		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);
	}

	//	if (pos_param.flgdupl) factdupl = 2;// default 0 : if flagged data are put in a duplicated map

	//	struct wcsprm * wcs;
	//	read_keyrec(dir.tmp_dir,wcs,&NAXIS1, &NAXIS2); // read wcs file header generated by sanePos
	//
	//	if(rank==0)
	//		cout << "\nMap size :" << NAXIS1 << "x" << NAXIS2 << endl;
	//
	//
	//	long long test_size; // used to verify indpsrc size
	//	if(read_indpsrc( test_size, npixsrc, indpsrc,  dir.tmp_dir)){ // read indpsrc (mask index) file generated by sanePos
	//#ifdef USE_MPI
	//		//		MPI_Barrier(MPI_COMM_WORLD);
	//		MPI_Finalize();
	//#endif
	//		return(EX_IOERR);
	//	}
	//	if(test_size != NAXIS1*NAXIS2){
	//		if(rank==0)
	//			cout << "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos" << endl;
	//#ifdef USE_MPI
	//		//		MPI_Barrier(MPI_COMM_WORLD);
	//		MPI_Finalize();
	//#endif
	//		return(EX_IOERR);
	//	}
	//	// each frame contains npixsrc pixels with index indsprc[] for which
	//	// crossing constraint are removed
	//	// thus
	//	// addnpix = number of pix to add in pixon
	//	//         = number of scans * number of pix in box crossing constraint removal
	//	addnpix = samples_struct.ntotscan*npixsrc;
	//
	//
	//	//read projection vector from a file
	//	if(read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon)){
	//#ifdef USE_MPI
	//		//		MPI_Barrier(MPI_COMM_WORLD);
	//		MPI_Finalize();
	//#endif
	//		return(EX_IOERR);
	//	}
	//
	//	// Check indpix readed size = expected size
	//	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
	//		if(rank==0){
	//			cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
	//			cout << ind_size << " != "  << (factdupl*NAXIS1*NAXIS2+2 + addnpix) << " " << factdupl << " " << addnpix << endl;
	//		}
	//#ifdef USE_MPI
	//		//		MPI_Barrier(MPI_COMM_WORLD);
	//		MPI_Finalize();
	//#endif
	//		return(EX_IOERR);
	//	}

	//
	//	// (At N-1 D) memory allocation
	//	PNdNaiv= new double[npix];
	//	hitsNaiv=new long[npix];
	//
	//	fill(PNdNaiv,PNdNaiv+npix,0.0);
	//	fill(hitsNaiv,hitsNaiv+npix,0);
	//
	//
	//#ifdef USE_MPI
	//
	//	if(rank==0){	// global (At N-1 D) malloc for mpi
	//		PNdtotNaiv = new double[npix];
	//		hitstotNaiv=new long[npix];
	//	}
	//
	//#endif


	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//
	//	if(rank==0)
	//		printf("\nComputing Naïve map\n");
	//
	//#ifdef USE_MPI
	//	MPI_Barrier(MPI_COMM_WORLD);
	//#endif
	//
	//
	//	// loop over the scans
	//	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
	//
	//		ns = samples_struct.nsamples[iframe]; // number of samples for this scan
	//		f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp; // knee freq of the filter in terms of samples in order to compute fft
	//		f_lppix_Nk = samples_struct.fcut[iframe]*double(ns)/proc_param.fsamp; // noise PS threshold freq, in terms of samples
	//
	//		string output_read = "";
	//		std::vector<string> det_vect;
	//		if(read_channel_list(output_read, samples_struct.bolovect[iframe], det_vect)){
	//			cout << output_read << endl;
	//#ifdef USE_MPI
	//			//			MPI_Barrier(MPI_COMM_WORLD);
	//			MPI_Finalize();
	//#endif
	//			return EX_CONFIG;
	//		}
	//
	//		long ndet = (long)det_vect.size();
	//
	//		int pb=0;
	//
	//#ifdef PARA_FRAME
	//		pb+=do_PtNd_Naiv(PNdNaiv, dir.tmp_dir, samples_struct.fitsvect, det_vect,ndet,proc_param.poly_order, proc_param.napod, f_lppix, ns, 0, 1, indpix, iframe, hitsNaiv);
	//		// Returns Pnd = (At N-1 d), Mp and hits
	//#else
	//
	//		pb+=do_PtNd_Naiv(PNdNaiv, dir.tmp_dir, samples_struct.fitsvect, det_vect,ndet, proc_param.poly_order, proc_param.napod, f_lppix, ns, rank, size, indpix, iframe, hitsNaiv);
	//
	//#endif
	//
	//		if(pb>0){
	//			cout << "Problem after do_PtNd_Naiv. Exiting...\n";
	//#ifdef USE_MPI
	//			//				MPI_Barrier(MPI_COMM_WORLD);
	//			MPI_Finalize();
	//#endif
	//			return -1;
	//		}
	//
	//	} // end of iframe loop
	//
	//
	//#ifdef USE_MPI
	//	MPI_Reduce(PNdNaiv,PNdtotNaiv,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	//	MPI_Reduce(hitsNaiv,hitstotNaiv,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	//
	//#else
	//	PNdtotNaiv=PNdNaiv;
	//	hitstotNaiv=hitsNaiv;
	//
	//#endif
	//
	//	if(rank==0)
	//		printf("\nEnd of Pre-Processing\n");
	//
	//	if (rank == 0){
	//
	//		string fnaivname = dir.output_dir + "naivMap.fits";
	//
	//		cout << "Output file : " << fnaivname << endl;
	//
	//		double *map1d;
	//		long long mi;
	//		map1d = new double[NAXIS1*NAXIS2];
	//
	//
	//
	//		for (long jj=0; jj<NAXIS2; jj++) {
	//			for (long ii=0; ii<NAXIS1; ii++) {
	//				mi = jj*NAXIS1 + ii;
	//				if (indpix[mi] >= 0){
	//					//					if(hitsNaiv[indpix[mi]]>0)
	//					map1d[mi] = PNdtotNaiv[indpix[mi]]/(double)hitstotNaiv[indpix[mi]];
	//				} else {
	//					map1d[mi] = NAN;
	//				}
	//			}
	//		}
	//
	//		if(write_fits_wcs("!" + fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Image",0)){ // open naive Map fits file and fill ultra naive map image
	//			cerr << "Error Writing Ultra Naiv map ... \n";
	//		}
	//
	//		for (long jj=0; jj<NAXIS2; jj++) {
	//			for (long ii=0; ii<NAXIS1; ii++) {
	//				mi = jj*NAXIS1 + ii;
	//				if (indpix[mi] >= 0){
	//					map1d[mi] = hitstotNaiv[indpix[mi]];
	//				} else {
	//					map1d[mi] = 0;
	//				}
	//			}
	//		}
	//
	//		if (addnpix){
	//			for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++){
	//				for (long jj=0; jj<NAXIS2; jj++) {
	//					for (long ii=0; ii<NAXIS1; ii++) {
	//						mi = jj*NAXIS1 + ii;
	//						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
	//						if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
	//							map1d[mi] += hitstotNaiv[indpix[ll]];
	//					}
	//				}
	//			}
	//		}
	//
	//		if(	write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Coverage",1)){ // open naive Map fits file and fill hit (or coverage) image
	//			cerr << "Error Writing coverage map  ... \n";
	//		}
	//
	//		if(write_fits_hitory2(fnaivname, NAXIS1, NAXIS2, dir.dirfile, proc_param, pos_param , samples_struct.fcut, samples_struct, structPS.ncomp)) // write sanePre parameters in naive Map fits file header
	//			cerr << "WARNING ! No history will be included in the file : " << fnaivname << endl;
	//		if (pos_param.maskfile != "")
	//			if(write_fits_mask(fnaivname, dir.input_dir + pos_param.maskfile)) // copy mask in naive map file
	//				cerr << "Warning ! The mask will not be included in naive map fits file ...\n";
	//
	//		delete [] map1d;
	//
	//	}
	/* ---------------------------------------------------------------------------------------------*/

	//Get processing time
	t3=time(NULL);

	//debug : computation time
#ifdef USE_MPI
	if(rank==0){
		printf("[%2.2i] Time : %d sec\n",rank, (int)(t3-t2));

		// clean up
		//		delete [] PNdtotNaiv;
		//		delete [] hitstotNaiv;
	}
#else
	cout << "Total Time : " << t3-t2 << " sec\n";
#endif

	// clean up
	//	delete [] PNdNaiv;
	//	delete [] hitsNaiv;
	//	delete [] indpix;
	//	delete [] indpsrc;

	//	wcsvfree(&nwcs, &wcs); // clean WCS structure

#ifdef USE_MPI
	MPI_Finalize();
#endif

	if (rank == 0)
		cout << endl << "done." << endl;

	return EXIT_SUCCESS;

}
