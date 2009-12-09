/*
 * sanepic_Main_Loop.cpp
 *
 *  Created on: 29 mai 2009
 *      Author: matthieu
 */



#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "imageIO.h"
#include "inline_IO2.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "parseSanepic.h"
#include "mpi_architecture_builder.h"

extern "C" {
#include "wcslib/wcshdr.h"
}



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

#else
	size = 1;
	rank = 0;
	cout << "Mpi will not be used for the main loop" << endl;
#endif



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

	//TODO : why do we need to define the default value here... should be in the parser...
	//DEFAULT PARAMETERS
	com.napod = 0; /*!  number of samples to apodize */
	u_opt.fsamp = 0.0; //25.0; /*!  sampling frequency : BLAST Specific */
	u_opt.projgaps = 0; /*!1: data flagged are put in a single pixel  (assume no signal in this pixel),
	0: data flagged are not reprojected */

	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if one sample is rejected, flagon=1 */
	int iterw = 10; /*!  period in iterations to which the data are written to disk, 0 = no intermediate map to be written*/
	u_opt.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	com.NOFILLGAP = 0; /*!  fill the gap ? default is YES (debug parameter) */
	com.flgdupl = 0; /*!  1 if flagged data are put in a separate map */
	u_opt.remove_polynomia = 1; /*! Remove a fitted polynomia from the data ? */
	u_opt.CORRon = 1; /*!  correlation included in the analysis (=1), else 0, default 0 */
	int factdupl = 1; /*! map duplication factor */
	long long addnpix=0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /*! number of pix in box constraint */



	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels */


	// map making parameters
	com.pixdeg=-1.0; /*! size of pixels (degree) */

	long long npix2; /*! used to check PNd reading was correct */
	long long ind_size; /*! indpix read size */
	long NAXIS1, NAXIS2;
	long long npix; /*! nn = side of the map, npix = number of filled pixels */


	//internal data params
	u_opt.f_lp=0.0; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples */




	double *PNdtot; /*! to deal with mpi parallelization : Projected noised data */
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */



	string field; /*! actual boloname in the bolo loop */
	string prefixe; /*! prefix used for temporary name file creation */

	std::vector<double> fcut; /*! noise cutting frequency vector */
	//std::vector<string> extentnoiseSP; /*! noise filenames vector of string */
	std::vector<struct box> boxFile;



	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	string fname; /*! parallel scheme filename */



	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		int parsed=1;

		parsed=parse_sanePic_ini_file(argv[1],u_opt,iterw, dir, samples_struct,com,
				det,boxFile, fcut, rank);
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


	//frames_index = new long [samples_struct.ntotscan];
	//extentnoiseSp_all = new string[samples_struct.ntotscan];


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	//vector2array(extentnoiseSP,  extentnoiseSp_all);
	cout << samples_struct.nsamples[0] << endl;


	/********************* Define parallelization scheme   *******/

#ifdef USE_MPI
	ofstream file;

	if(samples_struct.scans_index.size()==0){

		int test=0;
		fname = dir.outdir + parallel_scheme_filename;
		cout << fname << endl;

		test = define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
	}else{
		long size_tmp = 0;
		int return_error = 0;
		int num_frame = 0;
		char c;
		vector2array(samples_struct.scans_index,  samples_struct.index_table); // TODO : passer index_table en int plutot que long

		if(rank==0){
			//check the processor order given is correct
			//			size_tmp = *max_element(samples_struct.index_table, samples_struct.index_table+samples_struct.ntotscan);

			struct sortclass_long sortobject;
			sort(samples_struct.scans_index.begin(), samples_struct.scans_index.end(), sortobject);

			std::vector<long>::iterator it;
			//			int size_tmp=0;

			// using default comparison:
			it = unique(samples_struct.scans_index.begin(), samples_struct.scans_index.end());
			size_tmp = it - samples_struct.scans_index.begin();

			cout << "size unique : " << size_tmp << endl;

			cout << size << " vs size : " <<  size_tmp << endl;

			if((size_tmp)>size){
				cerr << "Number of processors are different between MPI and parallel scheme. Exiting\n";
				return_error =1;
			}else{

				samples_struct.scans_index.resize( size_tmp );

				cout << "trié + unique : " << samples_struct.scans_index[0] <<  " " << samples_struct.scans_index[1] << endl;


				if((size_tmp)<size){
					cout << "Warning. The number of processors used in fits_filelist is < to the number of processor used by MPI !\n";
					cout << "Do you wish to continue ? (y/n)\n";
					c=getchar();
					switch (c){
					case('y') :
						cout << "Let's continue with only " << (size_tmp) << " processor(s) !\n";
					break;
					default:
						cout << "Exiting ! Please modify fits filelist to use the correct number of processors\n";
						return_error =1;
						break;
					}

					for(long ii=0;ii<size_tmp;ii++)
						if(samples_struct.scans_index[ii]==0)
							num_frame++;

					if(num_frame==0){
						cout << "Exiting ! Please modify fits filelist to use at least processor 0 \n";
						return_error =1;
					}


				}else{


					for(long ii=0;ii<size_tmp;ii++)
						if(samples_struct.scans_index[ii]!=ii){
							cerr << "There is a problem in the fits filelist : you have forgot a processor to use. Exiting" << endl;
							return_error =1;
						}
				}
			}
		}





		if(rank==0){

			string outfile = dir.outdir + samples_struct.filename + "_sanepic.txt";
			cout << "outfile : " << outfile;
			file.open(outfile.c_str(), ios::out);
			if(!file.is_open()){
				cerr << "File [" << fname << "] Invalid." << endl;
				return_error = 1;
			}
		}


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&return_error,1,MPI_INT,0,MPI_COMM_WORLD);

		if(return_error>0){
			MPI_Finalize();
			exit(0);

		}

		string temp;
		size_t found;

		num_frame=0;
		iframe_min=0;
		iframe_max=0;

		long * nsamples_temp;
		nsamples_temp = new long[samples_struct.ntotscan];

		for(long jj = 0; jj<samples_struct.ntotscan; jj++)
			nsamples_temp[jj]= samples_struct.nsamples[jj];


		for(long ii = 0; ii<size; ii++){
			if(rank==ii)
				iframe_min=num_frame;
			for(long jj = 0; jj<samples_struct.ntotscan; jj++){
				if(samples_struct.index_table[jj]==ii){

					samples_struct.fits_table[num_frame]=samples_struct.fitsvect[jj];
					samples_struct.noise_table[num_frame]=samples_struct.noisevect[jj];
					samples_struct.nsamples[num_frame]=nsamples_temp[jj];
					if(rank==0){
						temp = samples_struct.fits_table[num_frame];
						found=temp.find_last_of('/');
						file << temp.substr(found+1) << " " << samples_struct.noise_table[num_frame] << " " << ii << endl;

					}
					num_frame++;
				}
			}
			if(rank==ii)
				iframe_max=num_frame;
		}

		delete [] nsamples_temp;

	}

	if(rank==0){
		file.close();
		cout << "on aura : \n";
		cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
		cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
		cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;
		//cout << samples_struct.filename << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // test
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
		//		MPI_Finalize();
		//		exit(0);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(long ii=0;ii<size;ii++){
		if(rank==ii)
			cout << "[ " << rank << " ]. iframe min : " << iframe_min << " iframemax : " << iframe_max << endl;
		else
			MPI_Barrier(MPI_COMM_WORLD);
	}

	//////// temp
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	MPI_Finalize();
	//	exit(0);

#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
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

	long long test_size;
	read_indpsrc( test_size, npixsrc, indpsrc,  dir.tmp_dir);
	if(test_size != NAXIS1*NAXIS2){
		cout << "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos" << endl;
		exit(0);
	}
	// each frame contains npixsrc pixels with index indsprc[] for which
	// crossing constraint are removed
	// thus
	// addnpix = number of pix to add in pixon
	//         = number of scans * number of pix in box crossing constraint removal
	addnpix = samples_struct.ntotscan*npixsrc;

	if (com.flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read indpix
	read_indpix(ind_size, npix, indpix,  dir.tmp_dir, flagon);

	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}

	// read npix, PNdtot from file
	read_PNd(PNdtot, npix2,  dir.tmp_dir);
	//	for (int ii=0;ii<20;ii++)
	//			cout << PNdtot[ii] << " ";
	//		cout << endl << "avant read indpix\n";
	//		getchar();


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

	FILE *fp;

	string testfile;
	ostringstream temp_stream;

	double *PtNPmatS,  *PtNPmatStot=NULL, *r, *q, *qtot=NULL, *d, *Mp, *Mptot=NULL, *s; // =NULL to avoid warnings
	// Mp = M in the paper = preconditioner
	long *hits, *hitstot=NULL;
	double *PNd;

	double var0 = 0.0, var_n = 0.0, delta0 = 0.0, delta_n = 0.0, alpha = 0.0;
	double delta_o, rtq, beta;

	long mi;
	double *map1d;

	int iter;
	double  f_lppix_Nk,f_lppix;
	long ns;
	long long npixeff;

	// memory allocs
	r           = new double[npix];
	q           = new double[npix];
	//	qtot        = new double[npix];
	d           = new double[npix];
	Mp          = new double[npix];
	//	Mptot       = new double[npix];
	s           = new double[npix];
	PtNPmatS    = new double[npix];
	//	PtNPmatStot = new double[npix];
	hits        = new long[npix];
	//	hitstot     = new long[npix];
	PNd         = new double[npix];

	map1d       = new double[NAXIS1*NAXIS2];

	for (int idupl = 0;idupl<=com.flgdupl;idupl++){



		//Conjugate gradien Inversion
		if (u_opt.projgaps || !flagon){
			npixeff = npix;
		} else {
			npixeff = npix-1;
		}



		printf("[%2.2i] npix = %lld, npixeff = %lld\n", rank, npix, npixeff);


		//t1 = time(0);

		fill(PtNPmatS,PtNPmatS+npix,0.0);
		fill(Mp,Mp+npix,0.0);
		fill(hits,hits+npix,0);
		fill(r,r+npix,0.0);
		fill(d,d+npix,0.0);
		fill(s,s+npix,0.0);
		fill(PNd,PNd+npix,0.0);



		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = samples_struct.nsamples[iframe];
			f_lppix_Nk = fcut[iframe]*double(ns)/u_opt.fsamp;


			// preconditioner computation : Mp
			if (u_opt.CORRon){
				write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,com.flgdupl, dir.tmp_dir,ns,iframe);
				// read pointing + deproject + fourier transform

				do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
						u_opt.fsamp,ns,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits);

				// return Pnd = At N-1 d
			} else {
				do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
						u_opt.fsamp, com.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
						iframe, PtNPmatS, Mp, hits);
			}

		} // end of iframe loop

#ifdef USE_MPI

		if(rank==0){
			PtNPmatStot = new double[npix];
			hitstot=new long[npix];
			Mptot = new double[npix];

			// TODO : Is it really necessary to initialize those variables ?
			fill(PtNPmatStot,PtNPmatStot+npix,0.0);
			fill(hitstot,hitstot+npix,0);
			fill(Mptot,Mptot+npix,0.0);
		}

		MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
		PtNPmatStot=PtNPmatS;
		hitstot=hits;
		Mptot=Mp;
#endif

		// intitialisation of the Conjugate gradient with preconditioner
		// see : http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
		if (rank == 0) {

			for (long ii=0;ii<npixeff;ii++)
				if (Mptot[ii] == 0)
					printf("ERROR: Mp[%ld] has elements = 0\n",ii);


			for (long ii=0;ii<npixeff;ii++)
				Mptot[ii] = 1.0/Mptot[ii]; // M : preconditioner


			for (long ii=0;ii<npixeff;ii++)
				r[ii] = PNdtot[ii] - PtNPmatStot[ii]; // r = b - Ax

			//			cout << r[0] << " " << r[1] << " " << r[2] << " " << endl;

			for (long ii=0;ii<npixeff;ii++)
				d[ii] =  Mptot[ii] * r[ii]; // d = M-1 * r


			delta_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				delta_n += r[ii]*d[ii]; // delta_new = rT * d

			//			cout << "delta_n : " << delta_n << endl;

			var_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				var_n += r[ii]*r[ii];


			delta0 = delta_n; // delta_0 <= delta_new
			var0 = var_n;
			printf("[%2.2i] var0 = %lf\n",rank, var0);

		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&var0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(d,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

		printf("[%2.2i] Main Conjugate gradient loop started\n",rank);


		//start loop
		iter = 0; // max iter = 2000, but ~100 iterations are required to achieve convergence

		// while i<imax and var_new > epsilon² * var_0 : epsilon² = 1e-10 => epsilon = 1e-5
		while(((iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !com.flgdupl)) || (!idupl && var_n/var0 > 1e-4)){ // 2000

			fill(q,q+npixeff,0.0); // q <= A*d

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = samples_struct.nsamples[iframe];
				//				ff = fframes[iframe];
				f_lppix_Nk = fcut[iframe]*double(ns)/u_opt.fsamp;

				if (u_opt.CORRon){
					write_tfAS(d,det,indpix,NAXIS1, NAXIS2,npix,com.flgdupl, dir.tmp_dir,ns,iframe);
					// read pointing + deproject + fourier transform

					do_PtNd(q, samples_struct.noise_table,dir.tmp_dir, "fPs_" ,det,f_lppix_Nk,
							u_opt.fsamp,ns,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNPS_nocorr(d, samples_struct.noise_table, dir, det,f_lppix_Nk,
							u_opt.fsamp, com.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
							iframe, q, NULL, NULL);
				}
			} // end of iframe loop




#ifdef USE_MPI
			if(rank==0){
				qtot = new double[npix];
				// TODO : Is it really necessary to initialize those variables ?
				fill(qtot,qtot+npix,0.0);
			}

			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			qtot=q;
#endif



			if (rank == 0){
				rtq= 0.0;
				for (long ii=0;ii<npixeff;ii++)
					rtq += qtot[ii] * d[ii]; // rtq = (dT * q)

				alpha = delta_n/rtq; // alpha <= delta_new / (dT * q)
				//				cout << "rtq : " << rtq << endl;
				//				cout << "alpha : " << alpha << endl;


				for (long ii=0;ii<npixeff;ii++)
					S[ii] += alpha*d[ii]; // x = x + alpha * d, x = S = signal
			}

#ifdef USE_MPI
			// Distribute the map to everyone
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif


			// every 10 iterations do ....
			if ((iter % 10) == 0){ // if iter is divisible by 10, recompute PtNPmatStot

				fill(PtNPmatS,PtNPmatS+npixeff,0.0);
				fill(PtNPmatStot,PtNPmatStot+npixeff,0.0);

				for (long iframe=iframe_min;iframe<iframe_max;iframe++){

					ns = samples_struct.nsamples[iframe];
					//					ff = fframes[iframe];
					f_lppix_Nk = fcut[iframe]*double(ns)/u_opt.fsamp;

					if (u_opt.CORRon){
						write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,com.flgdupl, dir.tmp_dir,ns,iframe);
						// read pointing + deproject + fourier transform

						do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir, "fPs_",det,f_lppix_Nk,
								u_opt.fsamp,ns,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
						// return Pnd = At N-1 d
					} else {

						do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
								u_opt.fsamp, com.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
								iframe, PtNPmatS, NULL, NULL);
					}
				} // end of iframe loop


#ifdef USE_MPI
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
				//				for(long ii=0;ii<npixeff;ii++)
				PtNPmatStot=PtNPmatS;
#endif


				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] = PNdtot[ii] - PtNPmatStot[ii]; //r = b - Ax
				}


			} else {

				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] -= alpha*qtot[ii]; // else r = r - alpha * q
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


				cout << "iter = " << iter;
				cout << ", crit  = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << var_n/var0;
				cout << ", crit2 = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << delta_n/delta0;
				cout << "\r " << flush;

				/*sprintf(testfile,"%s%s%s%s",outdir.c_str(),"Iteration_file_",termin.c_str(),".txt");
						      fp = fopen(testfile,"w");
						      if(fp!=NULL){
							for(ii=0;ii<npix;ii++)
							  fprintf(fp,"[%2.2i] iter = %d, crit = %10.15g, crit2 = %10.15g  \n",rank, iter,var_n/var0,delta_n/delta0);
							fclose(fp);}*/

				//      printf("[%2.2i] iter = %d, crit = %10.15g, crit2 = %10.15g     \n",rank, iter,var_n/var0,delta_n/delta0);


				//TODO : Do we really need that ? Should be covered by sanePre
				if (iter == 0){
					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]];
							} else {
								map1d[mi] = NAN;
							}
						}
					}


					fname = '!' + dir.outdir + "binMap_noisevar.fits"; // write preconditioner
					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					// TODO: Should not be to much different from a naive map...
					//					for (long ii=0; ii<NAXIS1 ; ii++){
					//						for (long jj=0; jj<NAXIS2; jj++){
					//							mi = jj*NAXIS1 + ii;
					//							if (indpix[mi] >= 0){
					//								map1d[mi] = Mptot[indpix[mi]] * PNdtot[indpix[mi]];
					//							} else {
					//								map1d[mi] = NAN;
					//							}
					//						}
					//					}
					//					fname = '!' + dir.outdir + "binMap_" + "_flux.fits";
					//					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					//					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					// TODO: Remove this, done in sanePre
					//
					//					for (long ii=0; ii<NAXIS1; ii++) {
					//						for (long jj=0; jj<NAXIS2; jj++) {
					//							mi = jj*NAXIS1 + ii;
					//							if (indpix[mi] >= 0){
					//								map1d[mi] = hitstot[indpix[mi]];
					//							} else {
					//								map1d[mi] = 0.0;
					//							}
					//						}
					//					}
					//					//TODO : Treat the 2 pixels properly
					//
					//					if (addnpix){
					//						for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
					//							for (long ii=0; ii<NAXIS1; ii++) {
					//								for (long jj=0; jj<NAXIS2; jj++) {
					//									mi = jj*NAXIS1 + ii;
					//									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc] != -1))
					//										map1d[mi] += hitstot[indpix[indpsrc[mi]+factdupl*NAXIS2*NAXIS2+iframe*npixsrc]];
					//								}
					//							}
					//						}
					//					}
					//
					//					fname = '!' + dir.outdir + "optimMap_" + "_hits.fits";
					//					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					//					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							map1d[mi] = 0.0;
						}
					}

					//TODO : Treat the 2 pixels properly
					if (addnpix){
						for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
									if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
										map1d[mi] += 1.0/Mptot[indpix[ll]];
								}
							}
						}

						fname = '!' + dir.outdir + "optimMap_" + "_invnoisevaruncpix.fits";
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
					}

				} // end of if (iter == 0)

				if (iterw && (iter % iterw) == 0){

					// make the map
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

					/*sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + "_flux" + iterstr + "b.fits";*/
					temp_stream << "!" + dir.outdir + "optimMap_" + "flux" << iter << "b.fits";


					// récupérer une chaîne de caractères
					fname= temp_stream.str();
					// Clear ostringstream buffer
					temp_stream.str("");
					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					if (com.flgdupl){
						for (long ii=0; ii<NAXIS1; ii++) {
							for (long jj=0; jj<NAXIS2; jj++) {
								mi = jj*NAXIS1 + ii;
								if (indpix[mi] >= 0){
									map1d[mi] = S[indpix[mi+NAXIS1*NAXIS2]]; //-finalmap[ii][jj];
								} else {
									map1d[mi] = 0.0;
								}
							}
						}


						/*sprintf(iterchar,"%d",iter);
							iterstr = iterchar;
							fname = '!' + outdir + "optimMap_" + "_fluxflags" + iterstr + "b.fits";*/
						temp_stream << "!" + dir.outdir + "optimMap_fluxflags_" << iter << "b.fits";

						// récupérer une chaîne de caractères
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
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


						/*sprintf(iterchar,"%d",iter);
							iterstr = iterchar;
							fname = '!' + outdir + "optimMap_" + termin + "_fluxuncpix_" + iterstr + "b.fits";*/
						temp_stream << "!" + dir.outdir + "optimMap_fluxuncpix_" << iter << "b.fits";

						// récupérer une chaîne de caractères
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					}
				} // end of if (iterw && (iter % iterw) == 0)


				//	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"ConvFile_",termin.c_str(),".txt");
				temp_stream << dir.outdir + "ConvFile.txt";

				// récupérer une chaîne de caractères
				testfile= temp_stream.str();
				// Clear ostringstream buffer
				temp_stream.str("");
				fp = fopen(testfile.c_str(),"a");
				fprintf(fp,"iter = %d, crit = %10.15g, crit2 = %10.15g\n",iter,var_n/var0, delta_n/delta0);
				fclose(fp);

			} // end of if (rank == 0)


#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++; // i = i +1

		} // end of while loop
		printf("\n");




		/*	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
			fp = fopen(testfile,"a");
			fprintf(fp,"test en sortie de la boucle while \n");
			fclose(fp);*/





		if  ((u_opt.projgaps || (com.flgdupl)) && !idupl){


			fill(PNd,PNd+npix,0.0);
			fill(PNdtot,PNdtot+npix,0.0);

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = samples_struct.nsamples[iframe];
				//				ff = fframes[iframe];
				f_lppix = u_opt.f_lp*double(ns)/u_opt.fsamp;
				f_lppix_Nk = fcut[iframe]*double(ns)/u_opt.fsamp;

				if (u_opt.CORRon){
					write_ftrProcesdata(S,u_opt,samples_struct,com,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
							npixsrc,addnpix,f_lppix,ns,	iframe);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)

					do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,"fdata_",det,f_lppix_Nk,
							u_opt.fsamp,ns,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);

					// return Pnd = At N-1 d
				} else {
					do_PtNd_nocorr(PNd,dir.tmp_dir,u_opt,samples_struct,com,det, f_lppix, f_lppix_Nk,
							addnpix, ns,indpix, indpsrc, NAXIS1, NAXIS2, npix, npixsrc, iframe, S);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			PNdtot=PNd;
#endif
		}



	}// end of idupl loop



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
	temp_stream << dir.outdir + "testfile.txt";

	// récupérer une chaîne de caractères
	testfile= temp_stream.str();
	// Clear ostringstream buffer
	temp_stream.str("");
	fp = fopen(testfile.c_str(),"a");
	fprintf(fp,"test avant ecriture \n");
	fclose(fp);


	//TODO : Should be outside of this function
	// write map function : NAXIS1, NAXIS2, S, indpix, outdir, termin, addnpix, ntotscan,
	// pixdeg, tancoord, tanpix, coordsyst, Mptot, indpsrc, npixsrc, factdupl,
	//

	//******************************  write final map in file ********************************



	if (rank == 0){

		printf(" after CC INVERSION %lld\n",npix*(npix+1)/2);

		//TODO : writing the maps should NOT be here...
		//TODO : In general separate the reading/writing from the computation

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

		fname = '!' + dir.outdir + "optimMap_flux.fits";
		write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = Mptot[indpix[mi]];
				} else {
					map1d[mi] = NAN;
				}
			}
		}


		fname = '!' + dir.outdir + "optimMap_noisevar.fits";
		//		write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
		write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

		if (addnpix){
			for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
				for (long ii=0; ii<NAXIS1; ii++) {
					for (long jj=0; jj<NAXIS2; jj++) {
						mi = jj*NAXIS1 + ii;
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1)  && (indpix[ll] != -1)){
							map1d[mi] = S[ll];
						} else {
							map1d[mi] = 0.0;
						}
					}
				}



				temp_stream << "!" + dir.outdir + "optimMap_flux_fr" << iframe << ".fits";
				// récupérer une chaîne de caractères
				fname= temp_stream.str();
				// Clear ostringstream buffer
				temp_stream.str("");
				write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

				for (long ii=0; ii<NAXIS1; ii++) {
					for (long jj=0; jj<NAXIS2; jj++) {
						mi = jj*NAXIS1 + ii;
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1)  && (indpix[ll] != -1)){
							map1d[mi] = Mptot[indpix[ll]];
						} else {
							map1d[mi] = 0.0;
						}
					}
				}

				//fname = '!' + outdir + "optimMap_" + termin + "_noisevar_fr" + iframestr + ".fits";
				temp_stream << "!" + dir.outdir + "optimMap_noisevar_fr" << iframe << ".fits";

				// récupérer une chaîne de caractères
				fname= temp_stream.str();
				// Clear ostringstream buffer
				temp_stream.str("");
				//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
				write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
			}
		}
	} // end of write map function

#ifdef USE_MPI
	delete [] qtot;
	delete [] Mptot;
	delete [] PtNPmatStot;
	delete [] hitstot;
#endif

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
	//*********************** End of program ***************************//
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

	delete [] samples_struct.nsamples;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;


	delete [] indpsrc;
	delete [] indpix;
	delete [] PNdtot;

	//delete [] frames_index;

	//	wcsfree(wcs);
	int nwcs = 1;
	wcsvfree(&nwcs, &wcs);




#ifdef USE_MPI
	MPI_Finalize();
#endif



	return 0;
}



