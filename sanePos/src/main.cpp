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

#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"

extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
}


#include "sanePos_map_making.h"
#include "blastSpecific.h"
#include "parsePos.h"
#include "sanePos_preprocess.h"


#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;



template<class T> void vector2array(std::vector<T> l, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}




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
{



	int size/*, size_det*/; /*! size = number of processor used for this step*/
	int rank/*, rank_det*/; /*! rank = processor MPI rank*/

#ifdef USE_MPI
	// int tag = 10;
	//MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	size = 1;
	rank = 0;
	cout << "Mpi is not used for this step" << endl;
#endif


	//struct user_options_sanepos u_opt;
	struct samples samples_struct;
	struct input_commons com;
	struct directories dir;
	struct detectors det;


	//default value of the data to pointing shift
	//u_opt.shift_data_to_point = 0; /*! default value = 0 */


	//DEFAULT PARAMETERS
	com.napod = 0; /*! number of samples to apodize, =0 -> no apodisation */
	//double errarcsec = 15.0; /*! source error, rejection criteria : scerr[ii] > errarcsec, sample is rejected */


	long iframe_min, iframe_max; /*! frame number min and max each processor has to deal with */

	int flagon = 0; /*! if rejectsample [ii]==3, flagon=1*/
	//u_opt.bfixc = 0; /*! indicates that 4 corners are given for the cross corelation removal box */
	bool pixout = 0; /*! indicates that at least one pixel has been flagged and is out */
	com.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	com.flgdupl = 0; /*! 1 if flagged data are put in a separate map */


	//set coordinate system
	double *coordscorner; /* srccoord = source coordinates, coordscorner = map corners coordinates*/
	//u_opt.srccoord = new double[2]; // RA/DEC source
	coordscorner = new double[4]; // map min/max RA/DEC coords (-N,-t,-T absents)
	//u_opt.srccoord[0] = -1000; // RA tangent point/source
	//u_opt.srccoord[1] = -1000; // DEC tangent point/source
	//u_opt.radius = -1.0; /*! map radius (half a side) in degrees */


	// data parameters
	//long *nsamples ; /*! number of samples table array */


	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels used*/
	//int nnf; /*! number of noise file */
	long long addnpix=0; /*!add a number 'n' of pixels to the map */



	// map making parameters
	//double pixdeg; /*! size of pixels (degree) */

	long long npix; /*! npix = number of filled pixels */
	long long npixsrc; /*! number of pixels included in Crossing Constraint Removal */

	struct wcsprm wcs;    // wcs structure of the image
	long NAXIS1, NAXIS2;  // size of the image

	double ra_min, ra_max, dec_min, dec_max; /*! ra/dec min/max coordinates of the map*/
	//float *scoffsets; /*! source offsets depending on wavelength */
	//scoffsets = new float[6];  useless now !

	//	int nfoff; /*! number of offsets */
	//	foffset *foffsets; /*! tableau d'offsets */

	//	double *tancoord; /*! tangent point coordinates RA/dec */
	//	double *tanpix; /*! tangent pixel coordinates in the map */
	double gra_min, gra_max, gdec_min, gdec_max; /*! global ra/dec min and max (to get the min and max of all ra/dec min/max computed by different processors) */

	//internal data params
	//	long ns; /*! number of samples for this scan */


	string fname; /*! parallel scheme file name */


	unsigned short *mask;
	long long *indpix, *indpsrc; /*! pixels indices, CCR mask pixels indices */

	long long *pixon; /*! this array is used to store the rules for pixels : they are seen or not */
	long long *pixon_tot;



	string field; /*! actual boloname in the bolo loop */
	string bolofield; /*! bolofield = boloname + bextension */
	string flagfield; /*! flagfield = field+fextension;*/
	//string dirfile; /*! data directory*/
	//string tmp_dir; /*! output directory*/
	//string poutdir; /*! current path (pPath) or output dir (outdir)*/

	/* parser inputs */
	//std::vector<string> bolonames/*, extentnoiseSP*/; /*! bolometer list, noise file prefix */
	std::vector<struct box> boxFile; /*! box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y) */
	//std::vector<string> fitsvect;
	//std::vector<string> noisevect;
	//std::vector<long> scans_index;

	time_t t2, t3;//, t3, t4, t5, dt;




	com.pixdeg = -1.0; /*! "Size of pixels (deg)"*/

	// -----------------------------------------------------------------------------//

	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		int parsed=1;
		/*parsed=parse_sanePos_ini_file(argv[1],u_opt,ntotscan,ndet,
				bolonames,nsamples,boxFile,fitsvect,scans_index);*/

		parsed=parse_sanePos_ini_file(argv[1],com,dir,
				det, samples_struct, boxFile);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}



	//fname = tmp_dir + parallel_scheme_filename;

	///////////////: debug ///////////////////////////////
	cout << "ntotscan : " << samples_struct.ntotscan << endl;
	/*

	std::vector<long>::iterator it;

	cout << "frames" << endl;
	for(it=fframes_vec.begin();it<fframes_vec.end();it++)
		cout << *it << " ";

	cout << "\nnsamples" << endl;
	for(it=nsamples_vec.begin();it<nsamples_vec.end();it++)
		cout << *it << " ";
	cout << endl;*/
	///////////////: debug ///////////////////////////////

	// -----------------------------------------------------------------------------//
	t2=time(NULL);



	// convert lists to regular arrays (MPI_BCas works only on array...
	//fframes       = new long[ntotscan];
	//nsamples      = new long[ntotscan];
	//string *fits_table, *noise_table;
	//long *index_table;

	samples_struct.fits_table  = new string[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];
	samples_struct.index_table = new long[samples_struct.ntotscan];

	//vector2array(fitsvect, fits_table);
	//vector2array(scans_index,  index_table);
	//cout << fframes[0] << fframes[1] << fframes[2] << endl;
	//cout << nsamples[0] << nsamples[1] << nsamples[2] << endl;
	//	cout << nsamples[0] << endl;

	/*
	if (NORMLIN)
		printf("NO BASELINE REMOVED\n");*/


	/*if (projgaps)
		printf("Flaged data are binned. iterative solution to fill gaps with noise only.\n");
	 */


#ifdef USE_MPI
	/********************* Define parallelization scheme   *******/

	//long *frnum;
	//frnum = new long[ntotscan+1];

	//	if (rank == 0){

	int test=0;
	fname = dir.outdir + parallel_scheme_filename;
	cout << fname << endl;
	test=define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct.ntotscan,size,samples_struct.nsamples,samples_struct.fitsvect,samples_struct.noisevect,samples_struct.fits_table, samples_struct.noise_table,samples_struct.index_table);

	if(test==-1){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}

	cout << "Et ca donne ca !" << endl;

	cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
	cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
	cout << samples_struct.index_table[0] << " " << samples_struct.index_table[1] << " " << samples_struct.index_table[2] << " " << samples_struct.index_table[3] << endl;
	cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;

	//	}

	//	MPI_Barrier(MPI_COMM_WORLD);
	//if(rank==0){
	//MPI_Bcast(nsamples,ntotscan,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
	// MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	//MPI_Bcast(index_table,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	//}

	cout << "mon rank : " << rank << endl;

	iframe_min = -1;
	//iframe_max = -1;

	for(long ii=0;ii<samples_struct.ntotscan;ii++){
		if((samples_struct.index_table[ii]==rank)&&(iframe_min == -1)){
			iframe_min=ii;
			break;
		}
	}

	iframe_max=iframe_min;
	for(iframe_max=iframe_min;iframe_max<samples_struct.ntotscan-1;iframe_max++)
		if(samples_struct.index_table[iframe_max]!=rank){
			iframe_max--;
			break;
		}

	iframe_max++;

	cout << rank << " iframe_min : " << iframe_min << endl;
	cout << rank << " iframe_max : " << iframe_max << endl;

	for(long ii=0;ii<samples_struct.ntotscan;ii++)
		samples_struct.fits_table[ii] = dir.dirfile + samples_struct.fits_table[ii];

	//exit(0);
	/*long *frnum ;

	  if (rank == 0){

		long *ruleorder ;
		long *fframesorder ;
		long *nsamplesorder ;
		//string *extentnoiseSp_allorder;

	  check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, &frnum);
	  // reorder nsamples
	  //find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
	  //cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";


	  fframesorder  = new long[ntotscan];
	  //extentnoiseSp_allorder = new string[ntotscan];
	  nsamplesorder = new long[ntotscan];

	  for (long ii=0;ii<ntotscan;ii++){
	  nsamplesorder[ii] = nsamples[ruleorder[ii]];
	  fframesorder[ii] = fframes[ruleorder[ii]];
	  //extentnoiseSp_allorder[ii] = extentnoiseSp_all[ruleorder[ii]];
	  }
	  for (long ii=0;ii<ntotscan;ii++){
	  nsamples[ii] = nsamplesorder[ii];
	  fframes[ii] = fframesorder[ii];
	  //extentnoiseSp_all[ii] = extentnoiseSp_allorder[ii];
	  //printf("frnum[%d] = %d\n",ii,frnum[ii]);
	  }

	  delete [] fframesorder;
	  delete [] nsamplesorder;
	  //delete [] extentnoiseSp_allorder;

	  delete [] ruleorder;


	  } else {
	  frnum = new long[ntotscan+1];
	  }
	 */
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//if(rank==0){
	// MPI_Bcast(nsamples,ntotscan,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
	// MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	//MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);
	//}

	//	iframe_min = frnum[rank];
	//iframe_max = frnum[rank+1];
	//rank_det = 0;
	//size_det = 1;
	//delete [] frnum;



#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

	//rank_det = rank;
	//size_det = size;

#endif




	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	/************************ Look for distriBoxution failure *******************************/
	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > samples_struct.ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}

	//	exit(0);

	/* END PARAMETER PROCESSING */




	/********** Allocate memory ***********/
	printf("[%2.2i] Allocating Memory\n",rank);

	// seek maximum number of samples
	//	ns = *max_elements(nsamples, nsamples+ntotscan);

	//	ra = new double[2*ns]; // RA bolo de ref
	//	dec = new double[2*ns]; // DEc du bolo de ref
	//	phi = new double[2*ns]; // (du bolo de ref) angle de la matrice de detecteur par rapport a RA/dec
	//scerr = new double[2*ns]; // BLAST SPECIFIC : mesure l'erreur de pointage, si trop grande on flag la donnée
	//	xx = new int[2*ns]; // sample column coordinates in the map
	//	yy = new int[2*ns]; // sample row coordinates in the map
	//	samptopix = new long[2*ns]; // sample to pixel conversion index
	//flag = new unsigned char[2*ns]; // flag data => =1
	//	flag = new short[2*ns];
	//	rejectsamp = new unsigned char[2*ns]; // rejected samples after flag conditions
	//flpoint = new unsigned char[2*ns]; // flpoint est un flag du pointage/time. Savoir au temps t, si tu prends ces données là, ou non.
	//	flpoint = new short[2*ns];
	//	tancoord = new double[2]; // coordinates in ra/dec of the tangent point
	//	tanpix = new double[2]; // coordinates in the map of the tangent point

	//	froffsets = new double[2]; //
	//	offsets = new double[2];

	//	offmap = new double[2]; // map offsets


	// default value for map variables
	//	ra_min  = 1000.0;
	//	ra_max  = -1000.0;
	//	dec_min = 1000.0;
	//	dec_max = -1000.0;

	/*
	offmap[0] = 0.0;
	offmap[1] = 0.0;
	 */

	//********************************************************************************
	//*************  find coordinates of pixels in the map
	//********************************************************************************

	printf("[%2.2i] Finding coordinates of pixels in the map\n",rank);

	//	bool default_projection = 1;

	// TODO: Different ways of computing the map parameters :
	// 1 - find minmax of the pointings on the sky -> define map parameters from that
	// 2 - defined minmax of the map -> define map parameters from that
	// (3 - define center of the map and radius -> define map parameters from that)


	computeMapMinima(det.boloname,samples_struct.fits_table,
			iframe_min,iframe_max,samples_struct.nsamples,com.pixdeg,
			ra_min,ra_max,dec_min,dec_max);

#ifdef USE_MPI
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
		printf("[%2.2i] ra  = [ %7.3f, %7.3f ] \n",rank, gra_min, gra_max );
		printf("[%2.2i] dec = [ %7.3f, %7.3f ] \n",rank, gdec_min, gdec_max);
	}


	computeMapHeader(com.pixdeg, (char *) "EQ", (char *) "TAN", coordscorner, wcs, NAXIS1, NAXIS2);

	if (rank == 0) {
		printf("[%2.2i] %ld x %ld pixels\n",rank, NAXIS1, NAXIS2);
		save_MapHeader(dir.tmp_dir,wcs, NAXIS1, NAXIS2);
	}

	//************************************* Deal with masking the point sources
	mask    = new unsigned short[NAXIS1*NAXIS2];
	indpsrc = new long long[NAXIS1*NAXIS2];

	// Initialize the masks
	addnpix=0;
	npixsrc=0;
	for (long long ii=0; ii<NAXIS1*NAXIS2; ii++){
		mask[ii]    =  1;
		indpsrc[ii] = -1;
	}

	// TODO : untested....
	// if a box for crossing constraint removal is given in ini file
	// TODO : save mask in fits file
	// TODO : being able to read a mask in fits file format
	for (long long iBox = 0; iBox < boxFile.size(); iBox++){
		for (long ii=boxFile[iBox].blc.x; ii<boxFile[iBox].trc.x ; ii++)
			for (long jj=boxFile[iBox].blc.y; jj<boxFile[iBox].trc.y; jj++){
				mask[jj*NAXIS1 + ii] = 0;
				indpsrc[jj*NAXIS1 + ii] = npixsrc++;
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
	(com.flgdupl) ? factdupl = 2: factdupl = 1; //  default 1 : if flagged data are put in a duplicated map

	// pixon indicates pixels that are seen
	// factdupl if flagged data are to be projected onto a separete map
	// 1 more pixel for flagged data
	// 1 more pixel for all data outside the map
	long long sky_size = factdupl*NAXIS1*NAXIS2 + 1 + 1 + addnpix;

	pixon = new long long[sky_size];
	pixon_tot = new long long[sky_size];
	fill(pixon,pixon+(sky_size),0);
	fill(pixon_tot,pixon_tot+(sky_size),0);


	//**********************************************************************************
	// Compute pixels indices
	//**********************************************************************************

	printf("[%2.2i] Compute Pixels Indices\n",rank);

	computePixelIndex(samples_struct.ntotscan,dir.tmp_dir, det.boloname,
			samples_struct.fits_table, iframe_min, iframe_max, samples_struct.nsamples,
			wcs, NAXIS1, NAXIS2,
			mask,
			com.napod, com.NOFILLGAP, com.flgdupl,factdupl,
			addnpix, pixon, rank,
			indpsrc, npixsrc, flagon, pixout);


#ifdef USE_MPI

	MPI_Reduce(pixon,pixon_tot,sky_size,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#else
	// TODO : pointer assignement only ?
	for(long long ii=0;ii<sky_size;ii++){
		pixon_tot[ii]=pixon[ii];
	}
#endif

	delete [] pixon;


	npix = 0;
	if(rank==0){

		indpix = new long long[sky_size];
		fill(indpix, indpix+(sky_size),-1);
		for(long long ii=0; ii< sky_size; ii++)
			if (pixon_tot[ii] != 0)
				indpix[ii] = npix++;

		/*!
		 * Write indpix to a binary file : ind_size = factdupl*nn*nn+2 + addnpix;
		 * npix : total number of filled pixels,
		 * flagon : if some pixels are apodized or outside the map
		 */
		write_indpix(sky_size, npix, indpix, dir.tmp_dir, flagon);
	}


	if (pixout)
		printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
	if(rank==0){
		printf("[%2.2i] Total number of detectors : %d\t Total number of Scans : %d \n",rank,(int)det.ndet, (int) samples_struct.ntotscan);
		printf("[%2.2i] Size of the map : %ld x %ld (using %lld pixels)\n",rank, NAXIS1, NAXIS2, sky_size);
		printf("[%2.2i] Total Number of filled pixels : %lld\n",rank, npix);
	}

	t3=time(NULL);
	cout << "temps de traitement : " << t3-t2 << " sec" << endl;




#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	/* ---------------------------------------------------------------------------------------------*/


	// Close MPI process


#ifdef USE_MPI
	MPI_Finalize();
#endif

	printf("[%2.2i] Cleaning up\n",rank);

	// TODO : Check all variable declaration/free
	// clean up
	delete [] mask;
	delete [] pixon_tot;
	delete [] coordscorner;
	delete [] samples_struct.nsamples;
	delete [] indpix;
	delete [] indpsrc;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.noise_table;
	delete [] samples_struct.index_table;

	printf("[%2.2i] End of sanePos\n",rank);

	return 0;
}
