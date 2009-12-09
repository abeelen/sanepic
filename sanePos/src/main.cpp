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
//#include <algorithm>
//#include <fstream>
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

	//TODO : Why the defaults are here ????
	//DEFAULT PARAMETERS
	com.napod  = 0; /*! number of samples to apodize, =0 -> no apodisation */
	com.pixdeg = -1.0; /*! "Size of pixels (deg)"*/

	long iframe_min=0, iframe_max=0; /*! frame number min and max each processor has to deal with */

	int flagon = 0; /*! if rejectsample [ii]==3, flagon=1*/
	bool pixout = 0; /*! indicates that at least one pixel has been flagged and is out */
	com.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	com.flgdupl = 0; /*! 1 if flagged data are put in a separate map */


	//set coordinate system
	double *coordscorner; /* srccoord = source coordinates, coordscorner = map corners coordinates*/
	coordscorner = new double[4]; // map min/max RA/DEC coords (-N,-t,-T absents)

	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels used*/

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
	long long *pixon_tot;

	string field; /*! actual boloname in the bolo loop */
	string bolofield; /*! bolofield = boloname + bextension */
	string flagfield; /*! flagfield = field+fextension;*/

	//TODO : Get rid of this
	std::vector<struct box> boxFile; /*! box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y) */

	//TODO : Debug this...
	time_t t2, t3;//, t3, t4, t5, dt;

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
				det, samples_struct, boxFile, rank);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}

	// -----------------------------------------------------------------------------//
	t2=time(NULL);

	samples_struct.fits_table  = new string[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];
	samples_struct.index_table = new long[samples_struct.ntotscan];



#ifdef USE_MPI

	ofstream file;


	// User has not given a processor order
	if(samples_struct.scans_index.size()==0){

		int test=0;
		fname = dir.outdir + parallel_scheme_filename;
		cout << fname << endl;
		//test=define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct.ntotscan,size,samples_struct.nsamples,samples_struct.fitsvect,samples_struct.noisevect,samples_struct.fits_table, samples_struct.noise_table,samples_struct.index_table);
		test = define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
		// user has given a processor order
	}else{
		int test=0;
		test = verify_parallelization_scheme(rank,dir.outdir,samples_struct, size, iframe_min, iframe_max);


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&test,1,MPI_INT,0,MPI_COMM_WORLD);

		if(test>0){
			MPI_Finalize();
			exit(0);

		}

	}

	if(rank==0){
		//				file.close();
		cout << "on aura : \n";
		cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
		cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
		cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;
		//cout << samples_struct.filename << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // test
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";

	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(long ii=0;ii<size;ii++){
		if(rank==ii)
			cout << "[ " << rank << " ]. iframemin : " << iframe_min << " iframemax : " << iframe_max << endl;
		else
			MPI_Barrier(MPI_COMM_WORLD);
	}

#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);


	ofstream file;
	string outfile = dir.outdir + parallel_scheme_filename;
	cout << "outfile : " << outfile << endl;
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


//	if (0){
		if(iframe_min!=iframe_max)
			printf("[%2.2i] Determining the size of the map\n",rank);

		// TODO: Different ways of computing the map parameters :
		// 1 - find minmax of the pointings on the sky -> define map parameters from that
		// 2 - defined minmax of the map -> define map parameters from that
		// (3 - define center of the map and radius -> define map parameters from that)



		//TODO : Should not be needed ! computeMapMinima should skip the loop
		if(iframe_min!=iframe_max)
//			computeMapMinima_HIPE(det.boloname,samples_struct,
//					iframe_min,iframe_max,com.pixdeg,
//					ra_min,ra_max,dec_min,dec_max);
			computeMapMinima(det.boloname,samples_struct,
					iframe_min,iframe_max,com.pixdeg,
					ra_min,ra_max,dec_min,dec_max);

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
			printf("[%2.2i] ra  = [ %7.3f, %7.3f ] \n",rank, gra_min, gra_max );
			printf("[%2.2i] dec = [ %7.3f, %7.3f ] \n",rank, gdec_min, gdec_max);
		}

		computeMapHeader(com.pixdeg, (char *) "EQ", (char *) "TAN", coordscorner, wcs, NAXIS1, NAXIS2);

		npixsrc = 0;
		// Initialize the masks
		mask    = new short[NAXIS1*NAXIS2];
		indpsrc = new long long[NAXIS1*NAXIS2];

		for (long long ii=0; ii<NAXIS1*NAXIS2; ii++){
			mask[ii]    =  0;
			indpsrc[ii] = -1;
		}

//	} else {
//		// Map header is determined from the mask file
//
//		if(iframe_min!=iframe_max)
//			printf("[%2.2i] Reading Mask map\n",rank);
//
//		string fname="mask_RCW120.fits";
//		string extname="mask";
//
//		if (read_mask_wcs(fname, extname, /*(char) 's',*/ wcs, NAXIS1, NAXIS2, mask ))
//			cerr << "Error Reading Mask file" << endl;
//
//		npixsrc = 0;
//		indpsrc = new long long[NAXIS1*NAXIS2];
//		long long ll;
//
//		for (long jj=0; jj<NAXIS2; jj++) {
//			for (long ii=0; ii<NAXIS1; ii++) {
//				ll = NAXIS1*jj+ii;
//				if (mask[ll] != 0)
//					indpsrc[ll] = npixsrc++;
//				else
//					indpsrc[ll] = -1;
//			}
//		}
//	}


	if (rank == 0) {
		printf("[%2.2i] %ld x %ld pixels\n",rank, NAXIS1, NAXIS2);
		save_MapHeader(dir.tmp_dir,wcs, NAXIS1, NAXIS2);
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
	fill(pixon,pixon+(sky_size),0);

	//**********************************************************************************
	// Compute pixels indices
	//**********************************************************************************

	if(iframe_min!=iframe_max)
		printf("[%2.2i] Compute Pixels Indices\n",rank);


	computePixelIndex_HIPE(dir.tmp_dir, det.boloname,samples_struct,
			com,iframe_min, iframe_max,
			wcs, NAXIS1, NAXIS2,
			mask,factdupl,
			addnpix, pixon, rank,
			indpsrc, npixsrc, flagon, pixout);



#ifdef USE_MPI
	if(rank==0){
		pixon_tot = new long long[sky_size];
		fill(pixon_tot,pixon_tot+(sky_size),0);
	}
	MPI_Reduce(pixon,pixon_tot,sky_size,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	delete [] pixon;
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
		write_indpix(sky_size, npix, indpix, dir.tmp_dir, flagon);
		write_indpsrc((long long) NAXIS1*NAXIS2, npixsrc, indpsrc,  dir.tmp_dir);

		delete [] indpsrc;
		delete [] indpix;
	}


	if (pixout)
		printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
	if(rank==0){
		printf("[%2.2i] Total number of detectors : %d\t Total number of Scans : %d \n",rank,(int)det.ndet, (int) samples_struct.ntotscan);
		printf("[%2.2i] Size of the map : %ld x %ld (using %lld pixels)\n",rank, NAXIS1, NAXIS2, sky_size);
		printf("[%2.2i] Total Number of filled pixels : %lld\n",rank, npix);
	}

	t3=time(NULL);


	if(iframe_min!=iframe_max){
		printf("[%2.2i] Temps de traitement : %d sec\n",rank,(int)(t3-t2));
		printf("[%2.2i] Cleaning up\n",rank);
	}

	// clean up
	delete [] mask;
	delete [] coordscorner;
	delete [] samples_struct.nsamples;

	//TODO : NO !
#ifdef USE_MPI
	if(rank==0)
		delete [] pixon_tot;
#else
	delete [] pixon;
#endif

	delete [] samples_struct.fits_table;
	delete [] samples_struct.noise_table;
	delete [] samples_struct.index_table;

	wcsfree(wcs);

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	if(iframe_min!=iframe_max)
		printf("[%2.2i] End of sanePos\n",rank);

	return 0;
}
