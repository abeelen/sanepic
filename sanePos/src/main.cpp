#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>


#include "mpi_architecture_builder.h"
#include "struct_definition.h"

#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"

extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}


#include "sanePos_map_making.h"
#include "sanePos_preprocess.h"
#include "parser_functions.h"


#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;




//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//


/*! \mainpage Sanepic
 *
 * \section intro_sec Matthieu HUSSON & Alexandre Beelen
 *
 * Sanepic USER Manual
 */



int main(int argc, char *argv[])
{



	int size; /*! size = number of processor used for this step*/
	int rank; /*! rank = processor MPI rank*/

#ifdef USE_MPI
	// int tag = 10;
	//MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size); // get mpi number of processors
	MPI_Comm_rank(MPI_COMM_WORLD,&rank); // each processor is given a number
	if(rank==0)
		printf("sanePos:\n");
#else
	size = 1; // non-MPI usage : 1 processor with number 0
	rank = 0;
	printf("sanePos:\n");
	cout << "Mpi is not used for this step" << endl;
#endif


	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_process proc_param;
	struct param_positions pos_param;
	struct common dir; /*! structure that contains output input temp directories */
	std::vector<detectors> detector_tab;
//	struct detectors det; /*! A structure that contains everything about the detectors names and number */


	long iframe_min=0, iframe_max=0; /*! frame number min and max each processor has to deal with */
	int flagon = 0; /*! if rejectsample [ii]==3, flagon=1*/
	bool pixout = 0; /*! indicates that at least one pixel has been flagged and is out */


	//set coordinate system
	double *coordscorner; /* srccoord = source coordinates, coordscorner = map corners coordinates*/
	coordscorner = new double[4]; // map min/max RA/DEC coords (-N,-t,-T absents)

	long long npix; /*! npix = number of filled pixels */
	long long npixsrc; /*! number of pixels included in Crossing Constraint Removal */
	long long addnpix=0; /*!add a number 'n' of pixels to the map */

	struct wcsprm * wcs;    // wcs structure of the image
	long NAXIS1, NAXIS2;  // size of the image


	// System should be IEEE 754 complient (TODO : add in the doc)
	double ra_min=NAN, ra_max=NAN, dec_min=NAN, dec_max=NAN; /*! ra/dec min/max coordinates of the map*/
	double gra_min, gra_max, gdec_min, gdec_max; /*! global ra/dec min and max (to get the min and max of all ra/dec min/max computed by different processors) */

	string fname; /*! parallel scheme file name */

	short *mask;
	long long *indpix, *indpsrc; /*! pixels indices, CCR mask pixels indices */

	long long *pixon; /*! this array is used to store the rules for pixels : they are seen or not */

	long long *pixon_tot=NULL;

	string field; /*! actual boloname in the bolo loop */
	string bolofield; /*! bolofield = boloname + bextension */
	string flagfield; /*! flagfield = field+fextension;*/

#ifdef DEBUG_PRINT
	time_t t2, t3;//, t3, t4, t5, dt;
#endif



	// -----------------------------------------------------------------------------//
	int parsed=0;
	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		parsed=-1;
	} else {
		std::vector<double> fcut;
		double fcut_sanePS=0.0;
		string MixMatfile, signame;
		long ncomp=1;
		int iterw=10;
		int save_data, restore;
		parsed=parser_function(argv[1], dir, detector_tab, samples_struct, pos_param, proc_param, fcut,
				fcut_sanePS, MixMatfile, signame, ncomp, iterw, save_data, restore, rank, size);
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


	if ((parsed>0)||(!compute_dirfile_format_file(dir.tmp_dir, samples_struct, detector_tab,rank))){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}
	// -----------------------------------------------------------------------------//
#ifdef DEBUG_PRINT
	t2=time(NULL);
#endif

	samples_struct.fits_table  = new string[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];
	samples_struct.index_table = new int[samples_struct.ntotscan];



#ifdef USE_MPI

	ofstream file;


	// User has not given a processor order
	if(samples_struct.scans_index.size()==0){

		int test=0;
		fname = dir.output_dir + parallel_scheme_filename;
		cout << fname << endl;
		test = define_parallelization_scheme(rank,fname,dir.input_dir,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
		// user has given a processor order
	}else{
		int test=0;
		test = verify_parallelization_scheme(rank,dir.output_dir,samples_struct, size, iframe_min, iframe_max);


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&test,1,MPI_INT,0,MPI_COMM_WORLD);

		if(test>0){
			MPI_Finalize();
			exit(0);

		}

	}

	//	if(rank==0){
	//		//				file.close();
	//		cout << "on aura : \n";
	//		cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
	//		cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
	//		cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;
	//		//cout << samples_struct.filename << endl;
	//	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // test
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";

	}

	MPI_Barrier(MPI_COMM_WORLD);

	//	for(long ii=0;ii<size;ii++){
	//		if(rank==ii)
	//			cout << "[ " << rank << " ]. iframemin : " << iframe_min << " iframemax : " << iframe_max << endl;
	//		else
	//			MPI_Barrier(MPI_COMM_WORLD);
	//	}

#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);


	ofstream file;
	string outfile = dir.output_dir + parallel_scheme_filename;
	//	cout << "outfile : " << outfile << endl;
	file.open(outfile.c_str(), ios::out);
	if(!file.is_open()){
		cerr << "File [" << outfile << "] Invalid." << endl;
		exit(0);
	}

	string temp;
	size_t found;

	for(long jj = 0; jj<samples_struct.ntotscan; jj++){

		temp = samples_struct.fits_table[jj];
		found=temp.find_last_of('/');
		file << temp.substr(found+1) << " " << samples_struct.noisevect[jj] << " 0" << endl;

	}

	file.close();

#endif



	/************************ Look for distriBoxution failure *******************************/
	if (iframe_min < 0 || iframe_min > iframe_max || iframe_max > samples_struct.ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}

	if (pos_param.maskfile == ""){

		//		if(iframe_min!=iframe_max)
		if(rank==0)
			printf("\n\nDetermining the size of the map\n");
		//	printf("[%2.2i] Determining the size of the map\n",rank);

		// TODO: Different ways of computing the map parameters :
		// 1 - find minmax of the pointings on the sky -> define map parameters from that
		// 2 - defined minmax of the map -> define map parameters from that
		// (3 - define center of the map and radius -> define map parameters from that)


		if(iframe_min!=iframe_max){
			switch (pos_param.fileFormat) {
			case 0:
				if(computeMapMinima(detector_tab,samples_struct,
						iframe_min,iframe_max,
						ra_min,ra_max,dec_min,dec_max)){
#ifdef USE_MPI
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Finalize();
#endif
					return(EXIT_FAILURE);
				}
				break;
			case 1:
				if(computeMapMinima_HIPE(detector_tab,samples_struct,
						iframe_min,iframe_max,
						ra_min,ra_max,dec_min,dec_max)){
#ifdef USE_MPI
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Finalize();
#endif
					return(EXIT_FAILURE);
				}
				break;
			}
		}

#ifdef USE_MPI

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&ra_min,&gra_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Reduce(&ra_max,&gra_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Reduce(&dec_min,&gdec_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Reduce(&dec_max,&gdec_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&gra_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&gra_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&gdec_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&gdec_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#else
		gra_min=ra_min;
		gra_max=ra_max;
		gdec_min=dec_min;
		gdec_max=dec_max;
#endif
		//set coordinates
		coordscorner[0] = gra_min; // store ra/dec min/max of the final map
		coordscorner[1] = gra_max;
		coordscorner[2] = gdec_min;
		coordscorner[3] = gdec_max;

		if (rank == 0) {
			printf("ra  = [ %7.3f, %7.3f ] \n", gra_min, gra_max);
			printf("dec = [ %7.3f, %7.3f ] \n", gdec_min, gdec_max);
		}

		computeMapHeader(pos_param.pixdeg, (char *) "EQ", (char *) "TAN", coordscorner, wcs, NAXIS1, NAXIS2);


		npixsrc = 0;
		// Initialize the masks
		mask    = new short[NAXIS1*NAXIS2];
		indpsrc = new long long[NAXIS1*NAXIS2];

		for (long long ii=0; ii<NAXIS1*NAXIS2; ii++){
			mask[ii]    =  0;
			indpsrc[ii] = -1;
		}

	} else {
		// Map header is determined from the mask file

		if(rank==0)
			cout << "Reading Mask map : " << pos_param.maskfile << endl;

		string extname="mask";

		if (read_mask_wcs(pos_param.maskfile, extname, /*(char) 's',*/ wcs, NAXIS1, NAXIS2, mask )){
			cerr << "Error Reading Mask file" << endl;
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}


		npixsrc = 0;
		indpsrc = new long long[NAXIS1*NAXIS2];
		long long ll;

		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				ll = NAXIS1*jj+ii;
				if (mask[ll] != 0)
					indpsrc[ll] = npixsrc++;
				else
					indpsrc[ll] = -1;
			}
		}
	}


	if (rank == 0) {
		printf("  Map Size : %ld x %ld pixels\n", NAXIS1, NAXIS2);
		if(save_MapHeader(dir.tmp_dir,wcs, NAXIS1, NAXIS2)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}
	}



	// each frame contains npixsrc pixels with index indsprc[] for which
	// crossing constraint are removed
	// thus
	// addnpix = number of pix to add in pixon
	//         = number of scans * number of pix in box crossing constraint removal
	addnpix = samples_struct.ntotscan*npixsrc;

	// map duplication factor
	int factdupl;
	(pos_param.flgdupl) ? factdupl = 2: factdupl = 1; //  default 1 : if flagged data are put in a duplicated map

	// pixon indicates pixels that are seen
	// factdupl if flagged data are to be projected onto a separate map
	// 1 more pixel for flagged data
	// 1 more pixel for all data outside the map
	long long sky_size = factdupl*NAXIS1*NAXIS2 + 1 + 1 + addnpix;

	pixon = new long long[sky_size];
	fill(pixon,pixon+(sky_size),0); // TODO test si ca marche avec des long long !

	//**********************************************************************************
	// Compute pixels indices
	//**********************************************************************************

	if(rank==0)
		printf("\n\nCompute Pixels Indices\n");

	switch (pos_param.fileFormat) {
	case 0:
		if(computePixelIndex(dir.tmp_dir, detector_tab,samples_struct,
				proc_param, pos_param, iframe_min, iframe_max,
				wcs, NAXIS1, NAXIS2,
				mask,factdupl,
				addnpix, pixon, rank,
				indpsrc, npixsrc, flagon, pixout)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}
		break;
	case 1:
		if(computePixelIndex_HIPE(dir.tmp_dir, detector_tab,samples_struct,
				proc_param, pos_param, iframe_min, iframe_max,
				wcs, NAXIS1, NAXIS2,
				mask,factdupl,
				addnpix, pixon, rank,
				indpsrc, npixsrc, flagon, pixout)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}
		break;
	}


#ifdef USE_MPI
	if(rank==0){
		pixon_tot = new long long[sky_size];
	}
	MPI_Reduce(pixon,pixon_tot,sky_size,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#else
	pixon_tot=pixon;
#endif


	npix = 0;
	if(rank==0){

		indpix = new long long[sky_size];
		for(long long ii=0; ii< sky_size; ii++)
			if (pixon_tot[ii] != 0)
				indpix[ii] = npix++;
			else
				indpix[ii] = -1;
		/*!
		 * Write indpix to a binary file : ind_size = factdupl*nn*nn+2 + addnpix;
		 * npix : total number of filled pixels,
		 * flagon : if some pixels are apodized or outside the map
		 */
		if(write_indpix(sky_size, npix, indpix, dir.tmp_dir, flagon)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}
		if(write_indpsrc((long long) NAXIS1*NAXIS2, npixsrc, indpsrc,  dir.tmp_dir)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}


		delete [] indpsrc;
		delete [] indpix;
	}


	if (pixout)
		printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
	if(rank==0){
		printf("Total number of Scans : %d \n", (int) samples_struct.ntotscan);
		printf("Size of the map : %ld x %ld (using %lld pixels)\n", NAXIS1, NAXIS2, sky_size);
		printf("Total Number of filled pixels : %lld\n", npix);
	}

#ifdef DEBUG_PRINT
	t3=time(NULL);
#endif


	if(rank==0){
#ifdef DEBUG_PRINT
		printf("\nProcessing time : %d sec\n",(int)(t3-t2));
#endif
		printf("\nCleaning up\n");
	}

	// clean up
	delete [] mask;
	delete [] coordscorner;
	delete [] samples_struct.nsamples;

	delete [] pixon;

	delete [] samples_struct.fits_table;
	delete [] samples_struct.noise_table;
	delete [] samples_struct.index_table;

	int nwcs=1;
	wcsvfree(&nwcs, &wcs);

	if(rank==0)
		printf("End of sanePos\n");


#ifdef USE_MPI
	if(rank==0)
		delete [] pixon_tot;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	return 0;
}
