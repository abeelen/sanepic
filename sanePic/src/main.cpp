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
#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"


#include "parseSanepic.h"
#include "sanepic_preprocess.h"
#include "conjugate_gradient.h"

//#include "estimPS_sanepic.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "mpi_architecture_builder.h"
#include <time.h>
//#include <fcntl.h>
//#include <unistd.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include <fftw3.h>
#include "wcslib/wcs.h"
}

#ifdef USE_MPI
#include "mpi.h"
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

	int size;//, size_det;
	int rank;//, rank_det;
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
	struct user_options u_opt;
	struct samples samples_struct;
	struct input_commons com;
	struct directories dir;
	struct detectors det;

	u_opt.projgaps = 0; /*!1: data flagged are put in a single pixel  (assume no signal in this pixel),
	0: data flagged are not reprojected */


	//u_opt.shift_data_to_point = 0; /*! default value of the data to pointing shift */

	//int samples_per_frames = 20; /*! blast Specific : Each frame has 20 samples */

	//DEFAULT PARAMETERS
	com.napod = 0; /*!  number of samples to apodize */
	u_opt.fsamp = 0.0; //25.0; /*!  sampling frequency : BLAST Specific */


	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if one sample is rejected, flagon=1 */
	int iterw = 10; /*!  period in iterations to which the data are written to disk, 0 = no intermediate map to be written*/
	u_opt.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	com.NOFILLGAP = 0; /*!  fill the gap ? default is YES (debug parameter) */
	//bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	com.flgdupl = 0; /*!  1 if flagged data are put in a separate map */
	u_opt.remove_polynomia = 1; /*! Remove a fitted polynomia from the data ? */
	u_opt.CORRon = 1; /*!  correlation included in the analysis (=1), else 0, default 0 */
	//bool parallel_frames = 0; // a parallel scheme is used : mpi has been launched
	int factdupl = 1; /*! map duplication factor */
	long long addnpix=0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /*! number of pix in box constraint */
	//bool bfixc = 0; // indicates that 4 corners are given for the cross corelation removal box


	// data parameters
	//	long *fframes  ; /*!  first frames table  */
	//long *nsamples ; /*!  number of samples (for each frame) table */

	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels */
	//int nnf; /*! extentnoiseSp_list number of elements */


	// map making parameters
	com.pixdeg=-1.0; /*! size of pixels (degree) */

	long long npix2; /*! used to check PNd reading was correct */
	long long ind_size; /*! indpix read size */
	long NAXIS1, NAXIS2;
	long long npix; /*! nn = side of the map, npix = number of filled pixels */
	//	double *tancoord; /*! tangent point coordinates */
	//	double *tanpix; /*! tangent pixel */

	//internal data params
	//long ns, ff; /*! number of samples for this scan, first frame number of this scan */
	u_opt.f_lp=0.0; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples */




	double *PNdtot; /*! to deal with mpi parallelization : Projected noised data */
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */



	string field; /*! actual boloname in the bolo loop */
	string *extentnoiseSp_all; /*! noise file name */
	//string bolofield; /*! bolofield = boloname + bextension */
	//string dirfile; /*! data directory */
	//string outdir; /*! output directory */
	//string tmp_dir; /*! temporary directory */
	//string poutdir; // current path (pPath) or output dir (outdir)
	//string bextension; /*! bolometer field extension */
	//string fextension = "NOFLAG"; /*! flag field extension */
	//string pextension; /*! pointing extension */
	//	string termin; /*! output file suffix */
	//string noiseSppreffile; /*! noise file suffix */
	string extentnoiseSp; /*! noise file */
	string prefixe; /*! prefix used for temporary name file creation */
	//	string termin_internal = "internal_data"; /*! internal data suffix */

	// utilisé lors de la lecture des coord de la map en pixel (dans la f° read_data)
	//string ra_field; /*! RA data file suffix */
	//string dec_field;/*! DEC data file suffix */
	//string phi_field;/*! PHI data file suffix */
	//string scerr_field = "ERR"+pextension; /*! Pointing error file suffix */
	//string flpoint_field = "FLPOINTING"; /*! pointing data file suffix */

	//string signame; /*! name of the signal file for sanePS */
	//string MixMatfile = "NOFILE"; /*! mixing matrix filename */
	//bool doInitPS = 0; /*! Do we rat a PS estimation from the elaborated map */

	/* DEFAULT PARAMETERS */
	//	int coordsyst = 1; /*! coordinatesystem :  Default is RA/DEC = 1 */
	//	int coordsyst2 = -1; /*! used to check binary reading of InfoPointing file */


	//std::vector<long> fframes_vec,nsamples_vec; /*! first frames number vector, number of samples vector */
	//std::vector<long>  xxi, xxf, yyi, yyf; /*! box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)*/
	std::vector<double> fcut; /*! noise cutting frequency vector */
	std::vector<string> extentnoiseSP; /*! noise filenames vector of string */
	//std::vector<string> bolonames; /*! bolonames vector */
	//std::vector<string> fitsvect, noisevect;
	//std::vector<long> scans_index;
	std::vector<struct box> boxFile;

	//time t2, t3, t4, t5, dt;



	//f_lp = 0.0; // low pass filter frequency
	//f_lp_Nk = 0.0; // noise PS frequency threshold
	//pixdeg = -1.0; // "Size of pixels (deg)"


	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	//string fname; /*! parallel scheme filename */



	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		//parse_sanePos_ini_file(argv[1]);
		int parsed=1;

		/*parsed=parse_sanePic_ini_file(argv[1],pixdeg,shift_data_to_point,napod,fsamp,NOFILLGAP,NORMLIN,projgaps,remove_polynomia,flgdupl,
				CORRon,iterw, ntotscan,ndet,f_lp,dirfile,outdir,tmp_dir,
				termin,MixMatfile,bolonames,fframes,nsamples,fname,xxi,xxf,yyi,yyf,fcut,extentnoiseSP, fitsvect, noisevect, scans_index);*/
		parsed=parse_sanePic_ini_file(argv[1],u_opt,iterw, dir, samples_struct,com,
				det,boxFile,extentnoiseSP, fcut);
		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}
		//exit(0);

	}
	////////////////////////////////////////////////////////////////

	if (u_opt.CORRon) printf("[%2.2i] CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);
	if (!(u_opt.CORRon)) printf("[%2.2i] NO CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);


	/*if (f_lp_Nk == 0.0)
		f_lp_Nk = f_lp;*/

	if (com.napod){
		printf("[%2.2i] Data are apodized\n", rank);
	} else {
		printf("[%2.2i] Data are not apodized\n", rank);
	}


	if (com.pixdeg < 0){
		cerr << "ERROR: enter pixel size -p keyword\n";
		exit(1);
	}

	if (u_opt.fsamp<=0.0){
		cerr << "ERROR: enter a correct sampling frequency -R keyword\n";
		exit(1);
	}



	//fframes  = new long[ntotscan];
	//nsamples = new long[ntotscan];
	extentnoiseSp_all = new string[samples_struct.ntotscan];


	//string *fits_table;
	//long *index_table;

	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];

	// convert vectors to regular arrays
	//vector2array(nsamples_vec, nsamples);
	//vector2array(fframes_vec,  fframes);

	vector2array(extentnoiseSP,  extentnoiseSp_all);

	//	cout << fframes[0] << endl;
	cout << samples_struct.nsamples[0] << endl;



	if (u_opt.NORMLIN)
		printf("NO BASELINE REMOVED\n");


	if (u_opt.projgaps)
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
	/********************* Define parallelization scheme   *******/

	//long *frnum;
	//frnum = new long[ntotscan+1];

	//	if (rank == 0){

	int test=0;
	string fname;
	fname = dir.outdir + parallel_scheme_filename;
	cout << fname << endl;
	test=define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct.ntotscan,size,samples_struct.nsamples,samples_struct.fitsvect,
			samples_struct.noisevect,samples_struct.fits_table, samples_struct.noise_table,samples_struct.index_table);

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

	for(int ii=0;ii<samples_struct.ntotscan;ii++){
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

	for(int ii=0;ii<samples_struct.ntotscan;ii++)
		samples_struct.fits_table[ii] = dir.dirfile + samples_struct.fits_table[ii];
#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

#endif

	//
	//	// allocate memory
	//	tancoord = new double[2];
	//	tanpix = new double[2];

	//	read_info_pointing(NAXIS1, NAXIS2, u_opt.outdir, tanpix, tancoord);
	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs, &NAXIS1, &NAXIS2);

	// read nn, coordsyst, tanpix, tancoord
	//	read_info_pointing(NAXIS1, NAXIS2, u_opt.tmp_dir, tanpix, tancoord);
	//cout << tanpix[0] << " " << tanpix[1] << endl;
	//cout << tancoord[0] << " " << tancoord[1] << endl;

	cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;


	// TODO : The mask should not be remade everytime...
	//******************************** some preprocess again  ****************/

	unsigned short *mask;
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
	for (unsigned long iBox = 0; iBox < boxFile.size(); iBox++){
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

	// compute indpsrc and addnpix
	//sanepic_preprocess(NAXIS1,NAXIS2, xxi, xxf, yyi, yyf, indpsrc, npixsrc, ntotscan, addnpix);


	if (com.flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read npix, PNdtot from file
	read_PNd(PNdtot, npix,  dir.tmp_dir);
	/*for (int ii=0;ii<20;ii++)
			cout << PNdtot[ii] << " ";
		cout << endl << "avant read indpix\n";
		exit(0);*/



	// read indpix
	read_indpix(ind_size, npix2, indpix,  dir.tmp_dir, flagon);

	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}
	if (npix!=npix2){
		cout << "Warning ! Indpix_for_conj_grad.bi and PNdCorr_*.bi are not compatible, npix!=npix2" << endl;
		exit(0);
	}


	/*************************************************************/

	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > samples_struct.ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}

	/* END PARAMETER PROCESSING */




	printf("[%2.2i] Main Conjugate gradient loop\n",rank);

	//MALLOC
	S = new double[npix];
	fill(S,S+npix,0.0);

	// conjugate GRADIENT LOOP
	sanepic_conjugate_gradient(com.flgdupl, npix, S, iframe_min, iframe_max,
			samples_struct.nsamples, fcut,u_opt.f_lp, u_opt.fsamp,
			indpix,
			wcs, NAXIS1, NAXIS2,
			factdupl, dir.tmp_dir, det.ndet,
			extentnoiseSp_all,dir.tmp_dir, det.boloname, iterw,
			indpsrc, npixsrc,flagon, u_opt.projgaps, rank, u_opt.CORRon,
			dir.dirfile, PNdtot, samples_struct.ntotscan,addnpix,u_opt.NORMLIN,com.NOFILLGAP,
			com.napod, u_opt.remove_polynomia, dir.outdir,samples_struct.fits_table);

	//
	//
	//string signame;
	//signame = tmp_dir + "Signal_internal_data.bin";
	//write_signal(npix, S, signame);



	//	//*******************************************************************//
	//	//******************  Update noise power spectra  *******************//
	//
	//	if (doInitPS){
	//		//printf("EstimPS will be run  with this mixing matrix file : %s\n",MixMatfile.c_str());
	//
	//		if (MixMatfile != "NOFILE"){
	//			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
	//				ns = nsamples[iframe];
	//				ff = fframes[iframe];
	//				extentnoiseSp = extentnoiseSp_all[iframe];
	//
	//				EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,S,
	//						/*MixMatfile,*/bolonames,dirfile,bextension,fextension,shift_data_to_point,
	//						tmp_dir,termin,termin_internal,NORMLIN,NOFILLGAP,remove_polynomia,tmp_dir,extentnoiseSp,outdir);
	//
	//			}
	//		}
	//	}



	//******************************************************************//
	//******************************************************************//
	//*********************** End of program *************************//
	//******************************************************************//
	//******************************************************************//




	// TODO : This will be rewrite differently
	//	if (rank == 0){
	//		//write infos for second part
	//		write_info_for_second_part(u_opt.outdir, NAXIS1, NAXIS2, npix,u_opt.pixdeg, tancoord, tanpix, coordsyst, flagon, indpix);
	//	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	// clean up
	delete [] S;
	//	delete [] fframes;
	delete [] samples_struct.nsamples;


	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;

	delete [] extentnoiseSp_all;
	//	delete [] tanpix;
	//	delete [] tancoord;
	delete [] indpsrc;
	delete [] indpix;
	delete [] PNdtot;





#ifdef USE_MPI
	MPI_Finalize();
#endif



	return 0;
}



