/*
 * sanepic_Main_Loop.cpp
 *
 *  Created on: 29 mai 2009
 *      Author: matthieu
 */





// liste des variables a donner :


// relancer le mpi et faire un best frame order

/*includes*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include "todprocess.h"
#include "map_making.h"

//#include "sane_io.h"
#include "binaryFileIO.h"
#include "boloIO.h"
#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"


#include "parseSanepic.h"
#include "sanepic_preprocess.h"
#include "conjugate_gradient.h"

#include "estimPS_sanepic.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "mpi_architecture_builder.h"
#include <time.h>
#include <fftw3.h>
//#include <fcntl.h>
//#include <unistd.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPI
#include "mpi.h"
#include "mpi_architecture_builder.h"
#endif

using namespace std;



/*
template<class T> void list2array(list<T> l, T* a)
{
	// copy list of type T to array of type T
	typename list<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}*/

/*
template<class T> void vector2array(std::vector<T> l, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}*/


//**********************************************************************************//
//**********************************************************************************//
//***************** Beginning of conjugate gradient program ************************//
//**********************************************************************************//
//**********************************************************************************//




int main(int argc, char *argv[])
{

	//

	int size, size_det;
	int rank, rank_det;
#ifdef USE_MPI
	// int tag = 10;
	MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	cout << size << endl;
	cout << rank << endl;

#else
	size = 1;
	rank = 0;
	cout << "Mpi will not be used for the main loop" << endl;
#endif

	/*char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);*/





	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//


	bool projgaps = 0; /*!1: data flagged are put in a single pixel  (assume no signal in this pixel),
	0: data flagged are not reprojected */


	int shift_data_to_point = 0; /*! default value of the data to pointing shift */

	int samples_per_frames = 20; /*! blast Specific : Each frame has 20 samples */

	//DEFAULT PARAMETERS
	long napod = 0; /*!  number of samples to apodize */
	double fsamp = 0.0; //25.0; /*!  sampling frequency : BLAST Specific */


	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if one sample is rejected, flagon=1 */
	int iterw = 10; /*!  period in iterations to which the data are written to disk, 0 = no intermediate map to be written*/
	bool NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	bool NOFILLGAP = 0; /*!  fill the gap ? default is YES (debug parameter) */
	//bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	bool flgdupl = 0; /*!  1 if flagged data are put in a separate map */
	bool remove_polynomia = 1; /*! Remove a fitted polynomia from the data ? */
	bool CORRon = 1; /*!  correlation included in the analysis (=1), else 0, default 0 */
	//bool parallel_frames = 0; // a parallel scheme is used : mpi has been launched
	int factdupl = 1; /*! map duplication factor */
	long addnpix=0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long npixsrc = 0; /*! number of pix in box constraint */
	//bool bfixc = 0; // indicates that 4 corners are given for the cross corelation removal box


	// data parameters
	long *fframes  ; /*!  first frames table  */
	long *nsamples ; /*!  number of samples (for each frame) table */

#ifdef USE_MPI
	long *fframesorder ; /*! reordered frames */
	long *nsamplesorder ; /*! reordered number of samples */
	long *ruleorder ; /*! which frame goes to which processor */
	long *frnum ;
	string *extentnoiseSp_allorder; /*! reordered noise file names */
#endif

	long ntotscan; /*! total number of scans */
	long ndet; /*! number of channels */
	int nnf; /*! extentnoiseSp_list number of elements */


	// map making parameters
	double pixdeg; /*! size of pixels (degree) */

	int npix2; /*! used to check PNd reading was correct */
	long ind_size; /*! indpix read size */
	int nn, npix; /*! nn = side of the map, npix = number of filled pixels */
	double *tancoord; /*! tangent point coordinates */
	double *tanpix; /*! tangent pixel */

	//internal data params
	long ns, ff; /*! number of samples for this scan, first frame number of this scan */
	double f_lp, f_lp_Nk; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples */




	double *PNdtot; /*! to deal with mpi parallelization : Projected noised data */
	long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */



	string field; /*! actual boloname in the bolo loop */
	string *extentnoiseSp_all; /*! noise file name */
	string bolofield; /*! bolofield = boloname + bextension */
	string dirfile; /*! data directory */
	string outdir; /*! output directory */
	string tmp_dir; /*! temporary directory */
	//string poutdir; // current path (pPath) or output dir (outdir)
	string bextension; /*! bolometer field extension */
	string fextension = "NOFLAG"; /*! flag field extension */
	string pextension; /*! pointing extension */
	string termin; /*! output file suffix */
	string noiseSppreffile; /*! noise file suffix */
	string extentnoiseSp; /*! noise file */
	string prefixe; /*! prefix used for temporary name file creation */
	string termin_internal = "internal_data"; /*! internal data suffix */

	// utilisé lors de la lecture des coord de la map en pixel (dans la f° read_data)
	string ra_field; /*! RA data file suffix */
	string dec_field;/*! DEC data file suffix */
	string phi_field;/*! PHI data file suffix */
	string scerr_field = "ERR"+pextension; /*! Pointing error file suffix */
	string flpoint_field = "FLPOINTING"; /*! pointing data file suffix */


	string MixMatfile = "NOFILE"; /*! mixing matrix filename */
	bool doInitPS = 0; /*! Do we rat a PS estimation from the elaborated map */

	/* DEFAULT PARAMETERS */
	int coordsyst = 1; /*! coordinatesystem :  Default is RA/DEC = 1 */
	int coordsyst2 = -1; /*! used to check binary reading of InfoPointing file */


	std::vector<long> fframes_vec,nsamples_vec; /*! first frames number vector, number of samples vector */
	std::vector<long>  xxi, xxf, yyi, yyf; /*! box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)*/
	std::vector<double> fcut; /*! noise cutting frequency vector */
	std::vector<string> extentnoiseSP; /*! noise filenames vector of string */
	std::vector<string> bolonames; /*! bolonames vector */


	//time t2, t3, t4, t5, dt;



	f_lp = 0.0; // low pass filter frequency
	f_lp_Nk = 0.0; // noise PS frequency threshold
	pixdeg = -1.0; // "Size of pixels (deg)"


	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	string fname; /*! parallel scheme filename */



	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		//parse_sanePos_ini_file(argv[1]);
		int parsed=1;
		parsed=parse_sanePic_ini_file(argv[1],pixdeg,shift_data_to_point,napod,fsamp,NOFILLGAP,NORMLIN,projgaps,remove_polynomia,flgdupl,
				CORRon,iterw,doInitPS, ntotscan,ndet,nnf,f_lp,f_lp_Nk,dirfile,outdir,tmp_dir,bextension,fextension,
				pextension,termin,noiseSppreffile,coordsyst,MixMatfile,bolonames,fframes_vec,nsamples_vec,fname,xxi,xxf,yyi,yyf,fcut,extentnoiseSP);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Finalize();
#endif
			exit(1);
		}
		//exit(0);

	}
	////////////////////////////////////////////////////////////////

	if (CORRon) printf("[%2.2i] CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);
	if (!CORRon) printf("[%2.2i] NO CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);


	if (f_lp_Nk == 0.0)
		f_lp_Nk = f_lp;

	if (napod){
		printf("[%2.2i] Data are apodized\n", rank);
	} else {
		printf("[%2.2i] Data are not apodized\n", rank);
	}


	if (pixdeg < 0){
		cerr << "ERROR: enter pixel size -p keyword\n";
		exit(1);
	}

	if (fsamp<=0.0){
		cerr << "ERROR: enter a correct sampling frequency -R keyword\n";
		exit(1);
	}


	fframes  = new long[ntotscan];
	nsamples = new long[ntotscan];
	extentnoiseSp_all = new string[ntotscan];




	// convert vectors to regular arrays
	vector2array(nsamples_vec, nsamples);
	vector2array(fframes_vec,  fframes);
	vector2array(extentnoiseSP,  extentnoiseSp_all);



	if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}


	if (coordsyst == 2){
		ra_field = "L"+pextension;
		dec_field = "B"+pextension;
		phi_field = "PHIG"+pextension;
		printf("[%2.2i] Coordinate system: Galactic\n",rank );
	}else{
		ra_field = "RA"+pextension;
		dec_field = "DEC"+pextension;
		phi_field = "PHI"+pextension;
		if (coordsyst == 3){
			printf("[%2.2i] Map in Telescope coordinates. Reference coordinate system is RA/DEC (J2000)\n", rank);
		} else {
			printf("[%2.2i] Coordinate system: RA/DEC (J2000)\n", rank);
		}
	}


	if (NORMLIN)
		printf("NO BASELINE REMOVED\n");


	if (projgaps)
		printf("Flaged data are binned. iterative solution to fill gaps with noise only.\n");


	// path in which data are written
	/*if (pPath != NULL){
		poutdir = pPath;
	} else {
		poutdir = outdir;
	}*/

	///////////////////////////////////////////////////////////////////



	/********************* Define parallelization scheme   *******/

#ifdef USE_MPI
	cout << "parallel_frames : 1 " << "   size : " << size  << endl;

	if (rank == 0){

		fframesorder = new long[ntotscan];
		nsamplesorder = new long[ntotscan];
		ruleorder = new long[ntotscan];
		frnum = new long[ntotscan+1];
		extentnoiseSp_allorder = new string[ntotscan];

		check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, &frnum);
		// reorder nsamples
		//find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";
		for (long ii=0;ii<ntotscan;ii++){
			nsamplesorder[ii] = nsamples[ruleorder[ii]];
			fframesorder[ii] = fframes[ruleorder[ii]];
			extentnoiseSp_allorder[ii] = extentnoiseSp_all[ruleorder[ii]];
		}
		for (long ii=0;ii<ntotscan;ii++){
			nsamples[ii] = nsamplesorder[ii];
			fframes[ii] = fframesorder[ii];
			extentnoiseSp_all[ii] = extentnoiseSp_allorder[ii];
			//printf("frnum[%d] = %d\n",ii,frnum[ii]);
		}

	}

	delete [] fframesorder;
	delete [] nsamplesorder;
	delete [] ruleorder;
	delete [] extentnoiseSp_allorder;


	MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);

	iframe_min = frnum[rank];
	iframe_max = frnum[rank+1];
	rank_det = 0;
	size_det = 1;
	delete [] frnum;

#else
	iframe_min = 0;
	iframe_max = ntotscan;
	rank_det = rank;
	size_det = size;
#endif


	// allocate memory
	tancoord = new double[2];
	tanpix = new double[2];


	// read nn, coordsyst, tanpix, tancoord
	read_info_pointing(nn, tmp_dir, termin_internal, coordsyst2, tanpix, tancoord);
	//cout << tanpix[0] << " " << tanpix[1] << endl;
	//cout << tancoord[0] << " " << tancoord[1] << endl;



	cout << "Map size :" << nn << "x" << nn << endl;

	if (coordsyst!=coordsyst2){
		cerr << "Error : coordinates systems must be the same for preprocessing and mapmaking" << endl;
		exit(0);
	}



	//******************************** some preprocess again  ****************/

	// MALLOC
	indpsrc = new long[nn*nn];

	// compute indpsrc and addnpix
	sanepic_preprocess(nn, xxi, xxf, yyi, yyf, indpsrc, npixsrc, ntotscan, addnpix);


	if (flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read npix, PNdtot from file
	read_PNd(PNdtot, npix, termin_internal, tmp_dir);
	/*for (ii=0;ii<20;ii++)
			cout << PNdtot[ii] << " ";
		cout << endl << "avant read indpix\n";
		exit(0);*/



	// read indpix
	read_indpix(ind_size, npix2, indpix, termin_internal, tmp_dir, flagon);

	if(ind_size!=(factdupl*nn*nn+2 + addnpix)){
		cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}
	if (npix!=npix2){
		cout << "Warning ! Indpix_for_conj_grad.bi and PNdCorr_*.bi are not compatible, npix!=npix2" << endl;
		exit(0);
	}


	/*************************************************************/

	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}

	/* END PARAMETER PROCESSING */




	printf("[%2.2i] Main Conjugate gradient loop\n",rank);

	//MALLOC
	S = new double[npix];
	fill(S,S+npix,0.0);

	// conjugate GRADIENT LOOP
	sanepic_conjugate_gradient(flgdupl, npix, S, iframe_min, iframe_max,
			nsamples, fframes, fcut,f_lp, fsamp, indpix, nn, factdupl, tmp_dir, termin, termin_internal,
			ndet,extentnoiseSp_all,noiseSppreffile, bolonames, size_det, rank_det, iterw,
			pixdeg,tancoord, tanpix,coordsyst,indpsrc, npixsrc,flagon, projgaps, rank, CORRon,
			dirfile, PNdtot, ntotscan,addnpix,NORMLIN,NOFILLGAP,napod,shift_data_to_point,
			remove_polynomia,fextension,bextension,flpoint_field,scerr_field, outdir);





	//*******************************************************************//
	//******************  Update noise power spectra  *******************//

	if (doInitPS){
		//printf("EstimPS will be run  with this mixing matrix file : %s\n",MixMatfile.c_str());

		if (MixMatfile != "NOFILE"){
			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = nsamples[iframe];
				ff = fframes[iframe];
				extentnoiseSp = extentnoiseSp_all[iframe];

				EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,S,
						/*MixMatfile,*/bolonames,dirfile,bextension,fextension,shift_data_to_point,
						tmp_dir,termin,termin_internal,NORMLIN,NOFILLGAP,remove_polynomia,noiseSppreffile,extentnoiseSp,outdir);

			}
		}
	}



	//******************************************************************//
	//******************************************************************//
	//*********************** End of program *************************//
	//******************************************************************//
	//******************************************************************//





	if (rank == 0){
		//write infos for second part
		write_info_for_second_part(outdir, termin, nn, npix,pixdeg, tancoord, tanpix, coordsyst, flagon, indpix);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	// clean up
	delete [] S;
	delete [] fframes;
	delete [] nsamples;
	delete [] extentnoiseSp_all;
	delete [] tanpix;
	delete [] tancoord;
	delete [] indpsrc;
	delete [] indpix;
	delete [] PNdtot;





#ifdef USE_MPI
	MPI_Finalize();
#endif



	return 0;
}



