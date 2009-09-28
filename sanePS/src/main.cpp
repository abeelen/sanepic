/*
 * main.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

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




#include "positionsIO.h"
#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"


#include "parsePS.h"
#include "estimPS.h"
#include "todprocess.h"

#ifdef USE_MPI
#include "mpi.h"
#include "mpi_architecture_builder.h"
#endif

using namespace std;



int main(int argc, char *argv[])
{


	// read framesorder

	int size;//,size_det;
	int rank;//,rank_det;
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


	//default value of the data to pointing shift
	int shift_data_to_point = 0;

	//DEFAULT PARAMETERS
	long napod = 0; // number of samples to apodize
	double fsamp = 0.0; //25.0; // sampling frequency : BLAST Specific
	bool NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	bool NOFILLGAP = 0; // fill the gap ? default is YES (debug parameter)
	bool remove_polynomia=1;
	int factdupl = 1;
	bool flgdupl = 0; // 1 if flagged data are put in a separate map
	int flagon;
	long ind_size;

	//int samples_per_frames=20;
	long *indpix;

	long ntotscan; // total number of scans
	long ndet; // number of channels


	// map making parameters
	//double pixdeg; // size of pixels (degree)


	int NAXIS1, NAXIS2, npix; // nn = side of the map, npix = number of filled pixels

	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan


	string field; // actual boloname in the bolo loop
	//string bolofield; // bolofield = boloname + bextension
	string dirfile; // data directory
	string outdir; // output directory
	string tmp_dir; // current path (pPath) or output dir (outdir)
	//string bextension; // bolometer field extension
	//string fextension = "NOFLAG"; // flag field extension
	//string pextension; // pointing extension
	//string termin; // output file suffix
	string noiseSppreffile; // noise file suffix
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation
	//string termin_internal = "internal_data";

	string MixMatfile = "NOFILE";


	std::vector<long> /*fframes_vec,nsamples_vec,*/xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)
	std::vector<string> extentnoiseSP;
	std::vector<string> bolonames;
	std::vector<string> fitsvect, noisevect;
		std::vector<long> scans_index;

	// data parameters
	long *fframes  ; // first frames table ff_in list -> fframes
	long *nsamples ; // number of samples table nf_in list -> nsamples


	string *extentnoiseSp_all;
	//time t2, t3, t4, t5, dt;


	long iframe_min, iframe_max;


	//pixdeg = -1.0; // "Size of pixels (deg)"


	// main loop variables
	double *S;

	// parallel scheme file
	string fname;
	string signame;


	int parsed=0;
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		parsed=parse_sanePS_ini_file(argv[1], shift_data_to_point, napod,fsamp, NOFILLGAP, NORMLIN,remove_polynomia, flgdupl,
				ntotscan, ndet, dirfile, outdir, tmp_dir, /*bextension,
				fextension,*/ termin, noiseSppreffile,
				bolonames, fframes,  nsamples, extentnoiseSP, MixMatfile,signame,fitsvect,noisevect,scans_index);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}



	if (napod){
		printf("[%2.2i] Data are apodized\n", rank);
	} else {
		printf("[%2.2i] Data are not apodized\n", rank);
	}

	if (flgdupl)
		factdupl=2;


	if (fsamp<=0.0){
		cerr << "ERROR: enter a correct sampling frequency -R keyword\n";
		exit(1);
	}


	//	fframes  = new long[ntotscan];
	//nsamples = new long[ntotscan];
	extentnoiseSp_all = new string[ntotscan];

	string *fits_table;
	long *index_table;

	fits_table = new string[ntotscan];
	index_table= new long[ntotscan];
	//convert vector to standard C array to speed up memory accesses
	//vector2array(nsamples_vec, nsamples);
	//vector2array(fframes_vec,  fframes);
	vector2array(fitsvect, fits_table);
	//vector2array(noisevect, );
	vector2array(scans_index,  index_table);
	vector2array(extentnoiseSP,  extentnoiseSp_all);

	cout << fframes[0] << endl;
	cout << nsamples[0] << endl;



	/*if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}*/






	//indpix=new long[factdupl*nn*nn+2 + addnpix];

	//int npix2;
	read_indpix(ind_size, npix, indpix, tmp_dir, flagon);

	S = new double[npix];


	//signame = "optimMap_sanepic_flux.fits";
	//signame = "NOSIGFILE";

	//First time run S=0, after sanepic, S = Pure signal
	if(signame != "NOSIGFILE"){
		// if second launch of estimPS, read S and nn in the previously generated fits map

		//int npix2;
		//long nn2;
		//read_signal(npix2,S,signame);
		read_fits_signal(signame, S, indpix, NAXIS1, NAXIS2, npix);
		FILE * fp;
		fp = fopen("test_signal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);

		fclose(fp);

		cout << setprecision(10) << S[0] << endl;
		cout <<  setprecision(10) << S[1] << endl;
	}else{
		// read nn in InfoPoiting


		int coordsyst2;
		//double *tancoord;
		//double *tanpix;
		// allocate memory
		//tancoord = new double[2];
		//tanpix = new double[2];


		// read nn, coordsyst, tanpix, tancoord
		read_info_pointing(NAXIS1, NAXIS2, tmp_dir, coordsyst2, NULL, NULL); //juste to read nn
		//cout << tanpix[0] << " " << tanpix[1] << endl;
		//cout << tancoord[0] << " " << tancoord[1] << endl;

		//delete [] tancoord;
		//delete [] tanpix;

		for(long ii=0;ii<npix;ii++)
			S[ii]=0.0;
	}


	//cout << "Map size :" << nn << "x" << nn << endl;
	//getchar();
	cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;
	getchar();

#ifdef USE_MPI
	/********************* Define parallelization scheme   *******/

	fname = tmp_dir + parallel_scheme_filename;

	int test=0;
	long *frnum;
	test=define_parallelization_scheme(rank,fname,&frnum,ntotscan,size,nsamples,fframes);

	if(test==-1){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}

	/*if (rank == 0){

		long *ruleorder ;
		long *fframesorder ;
		long *nsamplesorder ;
		string *extentnoiseSp_allorder;

		check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, &frnum);
		// reorder nsamples
		//find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";


		fframesorder  = new long[ntotscan];
		extentnoiseSp_allorder = new string[ntotscan];
		nsamplesorder = new long[ntotscan];

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

		delete [] fframesorder;
		delete [] nsamplesorder;
		delete [] extentnoiseSp_allorder;

		delete [] ruleorder;

	}*/



	MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);

	iframe_min = frnum[rank];
	iframe_max = frnum[rank+1];
	//rank_det = 0;
	//size_det = 1;

	delete [] frnum;

#else
	iframe_min = 0;
	iframe_max = ntotscan;
	//rank_det = rank;
	//size_det = size;
#endif

	string fits_filename;

	if (MixMatfile != "NOFILE"){
		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
			extentnoiseSp = extentnoiseSp_all[iframe];
			fits_filename=fits_table[iframe];

			cout << ff ;
			cout << " " ;
			cout << ns ;
			cout << " " ;
			cout << extentnoiseSp ;
			cout << endl;

			EstimPowerSpectra(fsamp, ns, ff, ndet, NAXIS1,NAXIS2, npix, napod,	iframe, flgdupl, factdupl, indpix,	S, MixMatfile, bolonames, dirfile, /*bextension,
					fextension,*/ shift_data_to_point, tmp_dir, NORMLIN, NOFILLGAP,remove_polynomia, noiseSppreffile,extentnoiseSp, outdir,fits_filename);
			// fsamp = bolometers sampling freq
			// ns = number of samples in the "iframe" scan
			// ff = first sample number
			// ndet = total number of detectors
			// nn = side of the map
			// npix = total number of filled pixels
			// napod = number of border pixels used to apodize data
			// iframe == scan number
			// flgdupl = flagged data map duplication indicator
			// factdupl = duplication factor (1 or 2)
			// indpix = pixels index
			// S = Pnd
			// MixMatfile = this file contains the number of components that interviene in the common-mode component of the noise
			// and the value of alpha, the amplitude factor which depends on detectors but not on time (see formulae (3) in "Sanepic:[...], Patanchon et al.")
			// bolonames = detectors names
			// dirfile = data directory
			// bextension = -B option : "_data" for example
			// fextension = "NOFLAG" or -G option ("_flag" for example)
			// cextension = "NOCALP" or -R option ("_calp" for example)
			// shift_data_to_point (default 0), for subtracting a time offset to the data to match the pointing
			// poutdir = outpout dir or current path (default)
			// termin = output file suffix
			// NORMLIN = baseline is remove from the data, default =0, option -L
			// NOFILLGAP = fill the gap ? default yes => 0
			// noiseSppreffile = noise power spectrum file suffix = path
			// extentnoiseSp = noise file
			// outdir = output directory
		}
	}




#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	//clean up
	delete [] fframes;
	delete [] nsamples;
	delete [] extentnoiseSp_all;
	delete [] S;

	return 0;
}
