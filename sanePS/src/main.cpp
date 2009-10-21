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

#include "mpi_architecture_builder.h"
#include "parsePS.h"
#include "estimPS.h"
#include "todprocess.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;



int main(int argc, char *argv[])
{


	// read framesorder

	int size;//,size_det;
	int rank;//,rank_det;
#ifdef USE_MPI
	// int tag = 10;
	//MPI_Status status;

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

	struct user_options u_opt;
	//default value of the data to pointing shift
	u_opt.shift_data_to_point = 0;

	//DEFAULT PARAMETERS
	u_opt.napod = 0; // number of samples to apodize
	u_opt.fsamp = 0.0; //25.0; // sampling frequency : BLAST Specific
	u_opt.NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	u_opt.NOFILLGAP = 0; // fill the gap ? default is YES (debug parameter)
	u_opt.remove_polynomia=1;
	int factdupl = 1;
	u_opt.flgdupl = 0; // 1 if flagged data are put in a separate map
	int flagon;
	long long ind_size;

	//int samples_per_frames=20;
	long long *indpix;

	long ntotscan; // total number of scans
	long ndet; // number of channels


	// map making parameters
	//double pixdeg; // size of pixels (degree)


	long NAXIS1, NAXIS2;
	long long npix; // nn = side of the map, npix = number of filled pixels

	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan


	string field; // actual boloname in the bolo loop
	//string bolofield; // bolofield = boloname + bextension
	//string dirfile; // data directory
	//string outdir; // output directory
	//string tmp_dir; // current path (pPath) or output dir (outdir)
	string ellFile; // file containing the ells
	//string bextension; // bolometer field extension
	//string fextension = "NOFLAG"; // flag field extension
	//string pextension; // pointing extension
	//string termin; // output file suffix
	//string noiseSppreffile; // noise file suffix
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation
	//string termin_internal = "internal_data";

	string MixMatfile = "NOFILE";


	//std::vector<long> /*fframes_vec,nsamples_vec,*/xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)
	std::vector<string> extentnoiseSP;
	std::vector<string> bolonames;
	std::vector<string> fitsvect, noisevect;
	std::vector<long> scans_index;

	// data parameters
//	long *fframes  ; // first frames table ff_in list -> fframes
	long *nsamples ; // number of samples table nf_in list -> nsamples


	string *extentnoiseSp_all;
	//time t2, t3, t4, t5, dt;


	long iframe_min, iframe_max;


	//pixdeg = -1.0; // "Size of pixels (deg)"


	// main loop variables
	double *S;

	// parallel scheme file
//	string fname;
	string signame;


	int parsed=0;
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		parsed=parse_sanePS_ini_file(argv[1], u_opt,
				ntotscan, ndet,
				bolonames, nsamples, extentnoiseSP, MixMatfile, ellFile, signame,
				fitsvect,noisevect, scans_index);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}

//TODO : Check in the parser if we have more than one bolometer....


	if (u_opt.napod){
		printf("[%2.2i] Data are apodized\n", rank);
	} else {
		printf("[%2.2i] Data are not apodized\n", rank);
	}

	if (u_opt.flgdupl)
		factdupl=2;


	if (u_opt.fsamp<=0.0){
		cerr << "ERROR: enter a correct sampling frequency -R keyword\n";
		exit(1);
	}


	//	fframes  = new long[ntotscan];
	//nsamples = new long[ntotscan];
	extentnoiseSp_all = new string[ntotscan];

	string *fits_table,*noise_table;
	long *index_table;

	fits_table = new string[ntotscan];
	index_table= new long[ntotscan];
	noise_table = new string[ntotscan];
	//convert vector to standard C array to speed up memory accesses
	//vector2array(nsamples_vec, nsamples);
	//vector2array(fframes_vec,  fframes);
	//vector2array(fitsvect, fits_table);
	//vector2array(noisevect, );
	//vector2array(scans_index,  index_table);
	vector2array(extentnoiseSP,  extentnoiseSp_all);

	//cout << fframes[0] << endl;
	cout << nsamples[0] << endl;



	/*if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}*/






	//indpix=new long[factdupl*nn*nn+2 + addnpix];

	//int npix2;

	// TODO : Should not be here
	// TODO : Ugly fix... remove this for the moment
	read_indpix(ind_size, npix, indpix, u_opt.tmp_dir, flagon);

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


		//		int coordsyst2;
		//double *tancoord;
		//double *tanpix;
		// allocate memory
		//tancoord = new double[2];
		//tanpix = new double[2];

		// TODO : Should change to use the wcs structure
		// read nn, coordsyst, tanpix, tancoord
//		read_info_pointing(NAXIS1, NAXIS2, u_opt.tmp_dir, NULL, NULL); //juste to read nn
		//cout << tanpix[0] << " " << tanpix[1] << endl;
		//cout << tancoord[0] << " " << tancoord[1] << endl;
		//	read_info_pointing(NAXIS1, NAXIS2, u_opt.outdir, tanpix, tancoord);
			struct wcsprm * wcs;
			read_MapHeader(u_opt.outdir,wcs, &NAXIS1, &NAXIS2);

		//delete [] tancoord;
		//delete [] tanpix;

		for(long ii=0;ii<npix;ii++)
			S[ii]=0.0;
	}


	//cout << "Map size :" << nn << "x" << nn << endl;
	//getchar();
	cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;
	//	getchar();

#ifdef USE_MPI
	/********************* Define parallelization scheme   *******/
	int test=0;
	string fname;
	fname = u_opt.outdir + parallel_scheme_filename;
	cout << fname << endl;
	test=define_parallelization_scheme(rank,fname,u_opt.dirfile,ntotscan,size,nsamples,fitsvect,noisevect,fits_table, noise_table,index_table);

	if(test==-1){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}

	cout << "Et ca donne ca !" << endl;

	cout << fits_table[0] << " " << fits_table[1] << " " << fits_table[2] << " " << fits_table[3] << endl;
	cout << noise_table[0] << " " << noise_table[1] << " " << noise_table[2] << " " << noise_table[3] << endl;
	cout << index_table[0] << " " << index_table[1] << " " << index_table[2] << " " << index_table[3] << endl;
	cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << nsamples[3] << endl;

	//	}

	//	MPI_Barrier(MPI_COMM_WORLD);
	//if(rank==0){
	//MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	// MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	//MPI_Bcast(index_table,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	//}

	cout << "mon rank : " << rank << endl;

	iframe_min = -1;
	//iframe_max = -1;

	for(int ii=0;ii<ntotscan;ii++){
		if((index_table[ii]==rank)&&(iframe_min == -1)){
			iframe_min=ii;
			break;
		}
	}

	iframe_max=iframe_min;
	for(iframe_max=iframe_min;iframe_max<ntotscan-1;iframe_max++)
		if(index_table[iframe_max]!=rank){
			iframe_max--;
			break;
		}

	iframe_max++;

	cout << rank << " iframe_min : " << iframe_min << endl;
	cout << rank << " iframe_max : " << iframe_max << endl;

	for(int ii=0;ii<ntotscan;ii++)
		fits_table[ii] = u_opt.dirfile + fits_table[ii];

#else
	iframe_min = 0;
	iframe_max = ntotscan;
	vector2array(fitsvect, fits_table);
	vector2array(scans_index,  index_table);

#endif

	string fits_filename;

	// TODO: useless test
	if (MixMatfile != "NOFILE"){
		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
//			ff = fframes[iframe];
			ff = iframe;
			extentnoiseSp = extentnoiseSp_all[iframe];
			fits_filename=fits_table[iframe];

			cout << ff ;
			cout << " " ;
			cout << ns ;
			cout << " " ;
			cout << extentnoiseSp ;
			cout << endl;

			EstimPowerSpectra(u_opt.fsamp, ns, ff, ndet, NAXIS1,NAXIS2, npix, u_opt.napod,
					iframe, u_opt.flgdupl, factdupl, indpix,
					S, MixMatfile, bolonames, u_opt.dirfile, ellFile,
					 u_opt.tmp_dir,
					u_opt.NORMLIN, u_opt.NOFILLGAP,u_opt.remove_polynomia, u_opt.tmp_dir,
					extentnoiseSp, u_opt.outdir,fits_filename);
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
			// tmp_dir = noise power spectrum file suffix = path
			// extentnoiseSp = noise file
			// outdir = output directory
		}
	}




#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	//clean up
//	delete [] fframes;
	delete [] nsamples;
	delete [] extentnoiseSp_all;
	delete [] S;

	delete [] fits_table;
	delete [] index_table;
	delete [] noise_table;

	return 0;
}
