#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

#include "imageIO.h"
#include "inline_IO2.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "parser_functions.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "write_maps_to_disk.h"

extern "C" {
#include "wcslib/wcshdr.h"
}



#if defined(PARA_BOLO) && ! defined(USE_MPI)
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
	if(rank==0)
		printf("\nsanepic_conjugate_gradient:\n\n");

#else
	size = 1;
	rank = 0;
	printf("\nsanepic_conjugate_gradient:\n\n");
	cout << "Mpi will not be used for the main loop" << endl;
#endif



	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//

	struct param_process proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct;  /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_positions pos_param; /*! A structure that contains user options about map projection and properties */
	struct common dir; /*! structure that contains output input temp directories */
	struct detectors det; /*! A structure that contains everything about the detectors names and number */

	int nwcs=1; // number of wcs : 1
	int iterw; // sanePic writes a temporary fits file (map) to disk each iterw iterations (conjugate gradient)
	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if one sample is rejected, flagon=1 */
	int factdupl = 1; /*! map duplication factor */
	long long addnpix=0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /*! number of pix in box constraint */


	// map making parameters
	long long npix2; /*! used to check PNd reading was correct */
	long long ind_size; /*! indpix read size */
	long NAXIS1, NAXIS2; // map dimensions
	long long npix; /*! nn = side of the map, npix = number of filled pixels */


	double *PNdtot; /*! to deal with mpi parallelization : Projected noised data */
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */


	string field; /*! actual boloname in the bolo loop */
	string prefixe; /*! prefix used for temporary name file creation */

	std::vector<double> fcut; /*! noise cutting frequency vector */



	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	string fname; /*! parallel scheme filename */


	int parsed=0;

	if (argc<2) // no enough arguments
		parsed=1;
	else{
		// Parse ini file

		// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
		double fcut_sanePS=0.0;
		string MixMatfile, ellFile, signame;
		long ncomp=1;

		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], dir, det, samples_struct, pos_param, proc_param, fcut,
				fcut_sanePS, MixMatfile, ellFile, signame, ncomp, iterw, rank, size);
	}

	if (parsed>0){ // error during parser phase
		if (rank==0)
			switch (parsed){

			case 1: printf("Please run %s using a *.ini file\n",argv[0]);
			break;

			case 2 : printf("Wrong program options or argument. Exiting !\n");
			break;

			case 3 : cerr << "You are using too many processors : " << size << " processors for only " << det.ndet << " detectors! Exiting...\n";
			break;

			default :;
			}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}

	string name_rank = dir.output_dir + "debug_sanePic.txt"; // log file name


#ifdef DEBUG
	std::ostringstream oss;
	oss << dir.output_dir + "debug_sanePre_" << rank << ".txt";
	name_rank = oss.str();

	ofstream file_rank;
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	file_rank.open(name_rank.c_str(), ios::out | ios::trunc); //& ios::trunc
	if(!file_rank.is_open()){
		cerr << "File [" << file_rank << "] Invalid." << endl;
		return -1;
	}
	file_rank << "Opening file for writing debug at " << asctime (timeinfo)  << endl;
	file_rank.close();

#endif

	////////////////////////////////////////////////////////////////

	// Memory allocation
	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new int[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	/********************* Define parallelization scheme   *******/

#if defined(USE_MPI) && !defined(PARA_BOLO)
	ofstream file;

	if(samples_struct.scans_index.size()==0){

		int test=0;
		fname = dir.output_dir + parallel_scheme_filename;
		cout << fname << endl;

		test = define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
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
#if defined(USE_MPI) && defined(PARA_BOLO)

	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

#else
	fname = dir.output_dir + parallel_scheme_filename;
	int test=0;
	test=check_ParallelizationScheme(fname,dir.dirfile,samples_struct,size);
	if (test==-1){
		if(rank==0)
			cerr << "erreur dans check_parallelizationScheme non-MPI " << endl;
		exit(0);
	}
#endif
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

#endif

	//	read_info_pointing(NAXIS1, NAXIS2, proc_param.output_dir, tanpix, tancoord);
	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs, &NAXIS1, &NAXIS2); // read keyrec file

	if(rank==0)
		cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl << endl; // print map size

	long long test_size;
	read_indpsrc( test_size, npixsrc, indpsrc,  dir.tmp_dir); // read mask index
	if(test_size != NAXIS1*NAXIS2){ // check size compatibility
		if(rank==0)
			cout << "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos" << endl;
		exit(0);
	}

	// each frame contains npixsrc pixels with index indsprc[] for which
	// crossing constraint are removed
	// thus
	// addnpix = number of pix to add in pixon
	//         = number of scans * number of pix in box crossing constraint removal
	addnpix = samples_struct.ntotscan*npixsrc;

	if (pos_param.flgdupl) factdupl = 2; // default 0 : if flagged data are put in a duplicated map

	// read indpix
	read_indpix(ind_size, npix, indpix,  dir.tmp_dir, flagon); // read map index

	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){ // check size compatibility
		if(rank==0)
			cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}

	// read (At N-1 d) from file
	read_PNd(PNdtot, npix2,  dir.tmp_dir);

	if (npix!=npix2){ // check size compatibility
		if(rank==0)
			cout << "Warning ! Indpix_for_conj_grad.bi and PNdCorr_*.bi are not compatible, npix!=npix2" << endl;
		exit(0);
	}


	/*************************************************************/

	//	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	if (iframe_min < 0 || iframe_min > iframe_max || iframe_max > samples_struct.ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}

	/* END PARAMETER PROCESSING */



	if(rank==0)
		printf("\nMain Conjugate gradient loop\n");

	//MALLOC
	S = new double[npix];
	fill(S,S+npix,0.0);

	//////////////////////////////////// Computing of sanePic starts here


	FILE *fp;
	string testfile; // log file to follow evolution of both criteria
	ostringstream temp_stream; // fits files filename string stream

	// inititialisation of the Conjugate gradient with preconditioner
	// see (for a complete description of the following variables) : http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
	double *PtNPmatS,  *PtNPmatStot=NULL, *r, *q, *qtot=NULL, *d, *Mp, *Mptot=NULL, *s; // =NULL to avoid warnings
	// Mp = M in the paper = preconditioner
	long *hits, *hitstot=NULL; // coverage map and global coverage for MPI
	double *PNd; // (At N-1 d)

	double var0 = 0.0, var_n = 0.0, delta0 = 0.0, delta_n = 0.0, alpha = 0.0; // conjugate gradient convergence criteria
	double delta_o, rtq, beta; // conjugate gradient needed parameters (see paper for a complete description)

	long mi;
	double *map1d; // buffer used to store maps before exporting to fits

	int iter; // conjugate gradient loop counter
	double  f_lppix_Nk,f_lppix; // noise cut-off frequency (in terms of samples number), filter cut-off freq (samples)
	long ns; // number of samples for the considered scan
	long long npixeff; // number of filled pixels

	// memory allocs
	r           = new double[npix];
	q           = new double[npix];
	d           = new double[npix];
	Mp          = new double[npix];
	s           = new double[npix];
	PtNPmatS    = new double[npix];
	hits        = new long[npix];
	PNd         = new double[npix];
	map1d       = new double[NAXIS1*NAXIS2];

	// in case flagged pixels are put in a duplicated map
	for (int idupl = 0;idupl<=pos_param.flgdupl;idupl++){

		if (pos_param.projgaps || !flagon){
			npixeff = npix;
		} else {
			npixeff = npix-1;
		}

		fill(PtNPmatS,PtNPmatS+npix,0.0);
		fill(Mp,Mp+npix,0.0);
		fill(hits,hits+npix,0);
		fill(r,r+npix,0.0);
		fill(d,d+npix,0.0);
		fill(s,s+npix,0.0);
		fill(PNd,PNd+npix,0.0);



		for (long iframe=iframe_min;iframe<iframe_max;iframe++){

			ns = samples_struct.nsamples[iframe];
			f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

			// preconditioner computation : Mp
			if (proc_param.CORRon){

				//TODO : What in the normal case ?

#if defined(USE_MPI) && !defined(PARA_BOLO)
				write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, 0, 1);
				// read pointing + deproject + fourier transform
#else
				write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
#endif


#ifdef DEBUG
				time ( &rawtime );
				timeinfo = localtime ( &rawtime );
				file_rank.open(name_rank.c_str(), ios::out | ios::app);
				file_rank << "rank " << rank << " a fini write_tfAS et attend a " << asctime (timeinfo) << " \n";
				file_rank.close();
#endif
#ifdef PARA_BOLO
				MPI_Barrier(MPI_COMM_WORLD);
#endif
#if defined(USE_MPI) && !defined(PARA_BOLO)
				do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
						proc_param.fsamp,ns, 0,1,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank);
				// return Pnd = At N-1 d
#else
				do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
						proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank);

#endif

			} else {


				do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
						proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
						iframe, PtNPmatS, Mp, hits,rank,size);
			}

		} // end of iframe loop

#ifdef USE_MPI

		if(rank==0){ // malloc only for first processor that reduces the data
			PtNPmatStot = new double[npix];
			hitstot=new long[npix];
			Mptot = new double[npix];
			qtot = new double[npix];

			fill(PtNPmatStot,PtNPmatStot+npix,0.0);
			fill(hitstot,hitstot+npix,0);
			fill(Mptot,Mptot+npix,0.0);
			fill(qtot,qtot+npix,0.0);
		}

		MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else

		PtNPmatStot=PtNPmatS; // in case non-MPI : pointer equality
		hitstot=hits;
		Mptot=Mp;

#endif

		// inititialisation of the Conjugate gradient with preconditioner
		// see : http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
		if (rank == 0) {

			for (long ii=0;ii<npixeff;ii++)
				if (Mptot[ii] == 0)
					printf("ERROR: Mp[%ld] has elements = 0\n",ii);
				else
					Mptot[ii] = 1.0/Mptot[ii]; // M : preconditioner


			for (long ii=0;ii<npixeff;ii++)
				r[ii] = PNdtot[ii] - PtNPmatStot[ii]; // r = b - Ax


			for (long ii=0;ii<npixeff;ii++)
				d[ii] =  Mptot[ii] * r[ii]; // d = M-1 * r


			delta_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				delta_n += r[ii]*d[ii]; // delta_new = rT * d

			var_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				var_n += r[ii]*r[ii];


			delta0 = delta_n; // delta_0 <= delta_new
			var0 = var_n;
			printf("var0 = %lf\n",var0);

		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&var0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(d,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

		//		printf("[%2.2i] Main Conjugate gradient loop started\n",rank);


		//start loop
		iter = 0; // max iter = 2000, but ~100 iterations are required to achieve convergence

		// while i<imax and var_new > epsilon² * var_0 : epsilon² = 1e-10 => epsilon = 1e-5
		while( ( (iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !pos_param.flgdupl) )
				|| (!idupl && var_n/var0 > 1e-6) ){

			fill(q,q+npixeff,0.0); // q <= A*d

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = samples_struct.nsamples[iframe];
				f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

				if (proc_param.CORRon){

#if defined(USE_MPI) && !defined(PARA_BOLO)
					write_tfAS(d,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, 0, 1);
					// read pointing + deproject + fourier transform
#else
					write_tfAS(d,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);

#endif

#ifdef DEBUG
					time ( &rawtime );
					timeinfo = localtime ( &rawtime );
					file_rank.open(name_rank.c_str(), ios::out | ios::app);
					file_rank << "rank " << rank << " a fini write_tfAS et attend a " << asctime (timeinfo) << " \n";
					file_rank.close();
#endif
#ifdef PARA_BOLO
					MPI_Barrier(MPI_COMM_WORLD);
#endif

#if defined(USE_MPI) && !defined(PARA_BOLO)
					do_PtNd(q, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
							proc_param.fsamp,ns, 0,1,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
#else
					do_PtNd(q, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);

#endif

				} else {


					do_PtNPS_nocorr(d, samples_struct.noise_table, dir, det,f_lppix_Nk,
							proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
							iframe, q, NULL, NULL,rank,size);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			qtot=q;
#endif



			if (rank == 0){
				rtq= 0.0;
				for (long ii=0;ii<npixeff;ii++)
					rtq += qtot[ii] * d[ii]; // rtq = (dT * q)

				alpha = delta_n/rtq; // alpha <= delta_new / (dT * q)


				for (long ii=0;ii<npixeff;ii++)
					S[ii] += alpha*d[ii]; // x = x + alpha * d, x = S = signal
			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			//cout << rank << " S bcast\n";
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif


			// every 10 iterations do ....
			if ((iter % 10) == 0){ // if iter is divisible by 10, recompute PtNPmatStot

				fill(PtNPmatS,PtNPmatS+npixeff,0.0);

				for (long iframe=iframe_min;iframe<iframe_max;iframe++){
					ns = samples_struct.nsamples[iframe];
					f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;


					if (proc_param.CORRon){


#if defined(USE_MPI) && !defined(PARA_BOLO)
						write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, 0, 1);
						// read pointing + deproject + fourier transform
#else
						write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank,size);
#endif
#ifdef DEBUG
						time ( &rawtime );
						timeinfo = localtime ( &rawtime );
						file_rank.open(name_rank.c_str(), ios::out | ios::app);
						file_rank << "rank " << rank << " a fini write_tfAS et attend a " << asctime (timeinfo) << " \n";
						file_rank.close();
#endif
#ifdef PARA_BOLO
						MPI_Barrier(MPI_COMM_WORLD);
#endif
#if defined(USE_MPI) && !defined(PARA_BOLO)
						do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
								proc_param.fsamp,ns,0,1,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
						// return Pnd = At N-1 d
#else
						do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
								proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
#endif
					} else {


						do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
								proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
								iframe, PtNPmatS, NULL, NULL,rank,size);
					}
				} // end of iframe loop


#ifdef USE_MPI
				//				cout << rank << " PtNPmatS reduction\n";
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
				PtNPmatStot=PtNPmatS;
#endif


				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] = PNdtot[ii] - PtNPmatStot[ii]; //r = b - Ax
				}


			} else { // if iter is not divisible by 10 ...

				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] -= alpha*qtot[ii]; // r = r - alpha * q
				}
			}

			if (rank == 0){

				for (long ii=0;ii<npixeff;ii++)
					s[ii] = Mptot[ii]*r[ii]; // s = M-1 * r


				delta_o = delta_n; // delta_0 <= delta_new

				delta_n = 0.0;
				for (long ii=0;ii<npixeff;ii++)
					delta_n += r[ii]*s[ii]; // delta_new = rT * s

				var_n = 0.0;
				for (long ii=0;ii<npixeff;ii++)
					var_n += r[ii]*r[ii];



				beta = delta_n/delta_o; // beta = delta_new / delta_0
				for (long ii=0;ii<npixeff;ii++)
					d[ii] = s[ii] + beta*d[ii]; // d = s + beta * d



				if (iterw && (iter % iterw) == 0){ // saving iterated maps

					// Every iterw iteration compute the map and save it

					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = S[indpix[mi]];
							} else {
								map1d[mi] = NAN;
							}
						}
					}

					temp_stream << dir.output_dir + "optimMap_" << iter << "b.fits";
					fname= temp_stream.str();
					temp_stream.str("");
					write_fits_wcs("!" + fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,(char *) "Image",0);

					if (pos_param.flgdupl){
						for (long ii=0; ii<NAXIS1; ii++) {
							for (long jj=0; jj<NAXIS2; jj++) {
								mi = jj*NAXIS1 + ii;
								if (indpix[mi] >= 0){ // pixel observed in first map
									if (indpix[mi+NAXIS1*NAXIS2] >=0)
										map1d[mi] = S[indpix[mi + NAXIS1*NAXIS2]]; //-finalmap[ii][jj];
									else
										map1d[mi] = INFINITY;
								} else {
									map1d[mi] = NAN;
								}
							}
						}


						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Duplicated Image",1);
					}



					if (addnpix){
						// initialize the container
						for (long jj=0; jj<NAXIS2 ; jj++){
							for (long ii=0; ii<NAXIS1 ; ii++){
								mi = jj*NAXIS1 + ii;
								map1d[mi] = 0.0;
							}
						}
						// loop thru frame to coadd all pixels
						for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
									if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
										map1d[mi] += S[indpix[ll]]/Mptot[mi];;
								}
							}
						}
						// replace the non observed pixels by NAN
						for (long ii=0; ii<NAXIS1; ii++) {
							for (long jj=0; jj<NAXIS2; jj++) {
								mi = jj*NAXIS1 + ii;
								if (map1d[mi] == 0.0)
									map1d[mi] = NAN;
							}
						}

						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, "Flagged pixels temporary map", 1);

					}
					write_fits_hitory(fname , NAXIS1, NAXIS2, dir.dirfile, proc_param, pos_param , fcut, det, samples_struct);
				} // end of saving iterated maps


				cout << "iter = " << iter;
				cout << ", crit  = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << var_n/var0;
				cout << ", crit2 = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << delta_n/delta0;
				cout << "\r " << flush;

				temp_stream << dir.output_dir + "ConvFile.txt";

				// Transform into string
				testfile= temp_stream.str();

				// Clear ostringstream buffer
				temp_stream.str("");
				fp = fopen(testfile.c_str(),"a");
				fprintf(fp,"iter = %d, crit = %10.15g, crit2 = %10.15g\n",iter,var_n/var0, delta_n/delta0);
				fclose(fp);

			} // end of if (rank == 0)


#ifdef USE_MPI
			// BCast updated criteria and map d
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++; // i = i +1

		} // end of while loop
		printf("\n");








		// If the gaps are projected on a second map, redo the preprocessing/writing of fdata,
		// which will now replace the flagged data by data from the signal map

		if  ((pos_param.projgaps || (pos_param.flgdupl)) && !idupl){


			fill(PNd,PNd+npix,0.0);

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = samples_struct.nsamples[iframe];
				f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp;
				f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

				if (proc_param.CORRon){


#if defined(USE_MPI) && !defined(PARA_BOLO)
					write_ftrProcesdata(S,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
							npixsrc,addnpix,f_lppix,ns,	iframe, 0, 1, name_rank);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)
#else
					write_ftrProcesdata(S,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
							npixsrc,addnpix,f_lppix,ns,	iframe, rank, size, name_rank);
#endif

#ifdef DEBUG
					time ( &rawtime );
					timeinfo = localtime ( &rawtime );
					file_rank.open(name_rank.c_str(), ios::out | ios::app);
					file_rank << "rank " << rank << " a fini et attend a " << asctime (timeinfo) << " \n";
					file_rank.close();
#endif
#ifdef PARA_BOLO
					MPI_Barrier(MPI_COMM_WORLD);
#endif

#if defined(USE_MPI) && !defined(PARA_BOLO)
					do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,"fdata_",det,f_lppix_Nk,
							proc_param.fsamp,ns, 0, 1,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
					// return Pnd = At N-1 d
#else
					do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,"fdata_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
#endif
				} else {

					do_PtNd_nocorr(PNd,dir.tmp_dir,proc_param, pos_param, samples_struct,det, f_lppix, f_lppix_Nk,
							addnpix, ns,indpix, indpsrc, NAXIS1, NAXIS2, npix, npixsrc, iframe, S,rank,size);
				}

			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			PNdtot=PNd; // ajout Mat 02/07
#endif
		}



	}// end of idupl loop



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	//******************************  write final map in file ********************************

	if (rank == 0){
		printf(" after CC INVERSION %lld\n",npix*(npix+1)/2);

#ifdef DEBUG
		FILE * fp;
		fp = fopen("output_signal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);
		fclose(fp);
#endif

		write_maps_to_disk(S, NAXIS1, NAXIS2, dir.output_dir, indpix, indpsrc,
				Mptot, addnpix, npixsrc, factdupl, samples_struct.ntotscan,
				proc_param, pos_param, det, samples_struct, fcut,
				wcs, pos_param.maskfile);
	}// end of rank==0



#ifdef USE_MPI
	delete [] qtot;
	delete [] Mptot;
	delete [] PtNPmatStot;
	delete [] hitstot;
#endif



	// clean up of conjugate gradient variables
	delete [] r;
	delete [] q;
	delete [] d;
	delete [] Mp;
	delete [] s;
	delete [] PtNPmatS;
	delete [] hits;
	delete [] map1d;
	delete [] PNd;




	//******************************************************************//
	//******************************************************************//
	//*********************** End of sanePic ***************************//
	//******************************************************************//
	//******************************************************************//



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


#ifdef DEBUG
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	file_rank.open(name_rank.c_str(), ios::out | ios::app);
	if(!file_rank.is_open()){
		cerr << "File [" << file_rank << "] Invalid." << endl;
		return -1;
	}

	file_rank << "[ " << rank << " ] Finish Time : " << asctime (timeinfo) << endl;
	file_rank.close();
#endif


	// clean up
	delete [] S;

	delete [] samples_struct.nsamples;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;


	delete [] indpsrc;
	delete [] indpix;
	delete [] PNdtot;

	wcsvfree(&nwcs, &wcs);




#ifdef USE_MPI
	MPI_Finalize();
#endif

	cout << "\nEnd of sanePic" << endl;


	return 0;
}



