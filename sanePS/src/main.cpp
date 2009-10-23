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


	struct user_options u_opt;
	struct samples samples_struct;
	struct input_commons com;
	struct directories dir;
	struct detectors det;


	//DEFAULT PARAMETERS
	com.napod = 0; // number of samples to apodize
	u_opt.fsamp = 0.0; //25.0; // sampling frequency : BLAST Specific
	u_opt.NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	com.NOFILLGAP = 0; // fill the gap ? default is YES (debug parameter)
	u_opt.remove_polynomia=1;
	int factdupl = 1;
	com.flgdupl = 0; // 1 if flagged data are put in a separate map
	int flagon;
	long long ind_size;

	//int samples_per_frames=20;
	long long *indpix;

	samples_struct.ntotscan=0; // total number of scans
	det.ndet=0; // number of channels


	// map making parameters
	com.pixdeg=-1.0; // size of pixels (degree)


	long NAXIS1, NAXIS2;
	long long npix; // nn = side of the map, npix = number of filled pixels

	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan


	string field; // actual boloname in the bolo loop
	string ellFile; // file containing the ells
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation


	string MixMatfile = "NOFILE";
	string signame;

	//std::vector<long> /*fframes_vec,nsamples_vec,*/xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)
	std::vector<string> extentnoiseSP;
	//std::vector<string> bolonames;
	//std::vector<string> fitsvect, noisevect;
	//std::vector<long> scans_index;

	// data parameters
	//	long *fframes  ; // first frames table ff_in list -> fframes
	//long *nsamples ; // number of samples table nf_in list -> nsamples


	string *extentnoiseSp_all;
	//time t2, t3, t4, t5, dt;


	long iframe_min, iframe_max;


	//pixdeg = -1.0; // "Size of pixels (deg)"


	// main loop variables
	double *S;





	int parsed=0;
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		parsed=parse_sanePS_ini_file(argv[1], u_opt, dir, samples_struct,com,det,
				extentnoiseSP, MixMatfile, ellFile, signame);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}

	//TODO : Check in the parser if we have more than one bolometer....

	extentnoiseSp_all = new string[samples_struct.ntotscan];

	//string *fits_table,*noise_table;
	//long *index_table;

	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];

	//convert vector to standard C array to speed up memory accesses
	vector2array(extentnoiseSP,  extentnoiseSp_all);


	cout << samples_struct.nsamples[0] << endl;



	// TODO : Should not be here
	// TODO : Ugly fix... remove this for the moment
	read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon);

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

		// TODO : Should change to use the wcs structure
		// read nn, coordsyst, tanpix, tancoord
		//		read_info_pointing(NAXIS1, NAXIS2, u_opt.tmp_dir, NULL, NULL); //juste to read nn
		//cout << tanpix[0] << " " << tanpix[1] << endl;
		//cout << tancoord[0] << " " << tancoord[1] << endl;
		//	read_info_pointing(NAXIS1, NAXIS2, u_opt.outdir, tanpix, tancoord);
		struct wcsprm * wcs;
		read_MapHeader(dir.outdir,wcs, &NAXIS1, &NAXIS2);

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

	string fits_filename;

	// TODO: useless test
	if (MixMatfile != "NOFILE"){
		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = samples_struct.nsamples[iframe];
			//			ff = fframes[iframe];
			ff = iframe;
			extentnoiseSp = extentnoiseSp_all[iframe];
			fits_filename=samples_struct.fits_table[iframe];

			cout << ff ;
			cout << " " ;
			cout << ns ;
			cout << " " ;
			cout << extentnoiseSp ;
			cout << endl;

			EstimPowerSpectra(u_opt.fsamp, ns, ff, det.ndet, NAXIS1,NAXIS2, npix, com.napod,
					iframe, com.flgdupl, factdupl, indpix,
					S, MixMatfile, det.boloname, dir.dirfile, ellFile,
					dir.tmp_dir, u_opt.NORMLIN, com.NOFILLGAP,u_opt.remove_polynomia, dir.tmp_dir,
					extentnoiseSp, dir.outdir,fits_filename);
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
	delete [] samples_struct.nsamples;
	delete [] extentnoiseSp_all;
	delete [] S;

	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;

	return 0;
}
