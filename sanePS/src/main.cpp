/*
 * main.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>






#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"
#include "parsePS.h"
#include "estimPS.h"


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
	com.flgdupl = 0; // 1 if flagged data are put in a separate map
	int flagon;
	long long ind_size;
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


	long iframe_min=0, iframe_max=0;

	// main loop variables
	double *S;

	string fname;


	int parsed=0;
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		parsed=parse_sanePS_ini_file(argv[1], u_opt, dir, samples_struct,com,det,
				MixMatfile, ellFile, signame);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}

	//TODO : Check in the parser if we have more than one bolometer....
	long *frames_index;

	frames_index = new long [samples_struct.ntotscan];

	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];

	cout << samples_struct.nsamples[0] << endl;



	// TODO : Should not be here
	// TODO : Ugly fix... remove this for the moment
	read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon);

	S = new double[npix];



	//First time run S=0, after sanepic, S = Pure signal
	if(signame != "NOSIGFILE"){
		// if second launch of estimPS, read S and nn in the previously generated fits map

		read_fits_signal(signame, S, indpix, NAXIS1, NAXIS2, npix);
		FILE * fp;
		fp = fopen("test_signal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);

		fclose(fp);

		cout << setprecision(10) << S[0] << endl;
		cout <<  setprecision(10) << S[1] << endl;
	}else{

		// TODO : Should change to use the wcs structure
		struct wcsprm * wcs;
		read_MapHeader(dir.tmp_dir,wcs, &NAXIS1, &NAXIS2);


		for(long ii=0;ii<npix;ii++)
			S[ii]=0.0;
	}


	cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;

#ifdef USE_MPI


	int test=0;
	fname = dir.outdir + parallel_scheme_filename;
	cout << fname << endl;
	test = define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct,size, iframe_min, iframe_max);

	if(test==-1){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}

#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;
	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

	for(long ii=0; ii<samples_struct.ntotscan;ii++)
		frames_index[ii] = ii;

#endif

	string fits_filename;

	// TODO: useless test
	if (MixMatfile != "NOFILE"){
		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = samples_struct.nsamples[iframe];
			ff = iframe;
			extentnoiseSp = samples_struct.noise_table[iframe];
			fits_filename=samples_struct.fits_table[iframe];

			cout << ff ;
			cout << " " ;
			cout << ns ;
			cout << " " ;
			cout << extentnoiseSp ;
			cout << endl;

			EstimPowerSpectra(u_opt,det,dir,com, ns, ff, NAXIS1,NAXIS2, npix,
								iframe, indpix,	S, MixMatfile, ellFile,
								extentnoiseSp, fits_filename);
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
	delete [] samples_struct.nsamples;
	delete [] S;

	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;

	delete [] frames_index;

	return 0;
}
