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
#include "parseSanepic.h"
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

	struct param_process proc_param;
	struct samples samples_struct;
	struct param_positions pos_param;
	struct directories dir;
	struct detectors det;

	int nwcs=1;
	int iterw;
	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if one sample is rejected, flagon=1 */
	int factdupl = 1; /*! map duplication factor */
	long long addnpix=0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /*! number of pix in box constraint */


	// map making parameters
	long long npix2; /*! used to check PNd reading was correct */
	long long ind_size; /*! indpix read size */
	long NAXIS1, NAXIS2;
	long long npix; /*! nn = side of the map, npix = number of filled pixels */



	double *PNdtot; /*! to deal with mpi parallelization : Projected noised data */
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */


	string field; /*! actual boloname in the bolo loop */
	string prefixe; /*! prefix used for temporary name file creation */

	std::vector<double> fcut; /*! noise cutting frequency vector */
	//std::vector<string> extentnoiseSP; /*! noise filenames vector of string */



	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	string fname; /*! parallel scheme filename */


	int parsed=0;

	if (argc<2)
		parsed=1;
	else{
		// Parse ini file
		parsed=parse_sanePic_ini_file(argv[1],proc_param, pos_param, iterw, dir, samples_struct,
				det, fcut, rank);

		if(size>det.ndet) parsed=3;
	}

	if (parsed>0){
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


#ifdef DEBUG
	std::ostringstream oss;
	string name_rank;
	oss << dir.outdir + "debug_sanePre_" << rank << ".txt";
	name_rank = oss.str();
#else
	string name_rank = dir.outdir + "debug_sanePre.txt";

#endif

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

	////////////////////////////////////////////////////////////////


	//frames_index = new long [samples_struct.ntotscan];
	//extentnoiseSp_all = new string[samples_struct.ntotscan];


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new int[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	//vector2array(extentnoiseSP,  extentnoiseSp_all);
	//	cout << samples_struct.nsamples[0] << endl;


	/********************* Define parallelization scheme   *******/

#if defined(USE_MPI) && !defined(PARA_BOLO)
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
		vector2array(samples_struct.scans_index,  samples_struct.index_table);

		if(rank==0){
			//check the processor order given is correct
			//			size_tmp = *max_element(samples_struct.index_table, samples_struct.index_table+samples_struct.ntotscan);

			struct sortclass_int sortobject;
			sort(samples_struct.scans_index.begin(), samples_struct.scans_index.end(), sortobject);

			std::vector<int>::iterator it;
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


		for(int ii = 0; ii<size; ii++){
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
#if defined(USE_MPI) && defined(PARA_BOLO)
	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

	//	cout <<" final list : " << endl;
	//	for(int ii = 0; ii< samples_struct.ntotscan;ii++){
	//		cout << samples_struct.fits_table[ii] << " " << samples_struct.noise_table[ii] << " " << samples_struct.index_table[ii] << endl;
	//		//		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	//	}
#else
	fname = dir.outdir + parallel_scheme_filename;
	int test=0;
	test=check_ParallelizationScheme(fname,dir.dirfile,samples_struct,size);
	if (test==-1){
		if(rank==0)
			cerr << "erreur dans check_parallelizationScheme non-MPI " << endl;
		exit(0);
	}

	//	cout <<" final list : " << endl;
	//	for(int ii = 0; ii< samples_struct.ntotscan;ii++){
	//		cout << samples_struct.fits_table[ii] << " " << samples_struct.noise_table[ii] << " " << samples_struct.index_table[ii] << endl;
	//		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	//	}

#endif
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

#endif

	//	// allocate memory
	//	tancoord = new double[2];
	//	tanpix = new double[2];

	//	read_info_pointing(NAXIS1, NAXIS2, proc_param.outdir, tanpix, tancoord);
	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs, &NAXIS1, &NAXIS2);

	// read nn, coordsyst, tanpix, tancoord
	// read_info_pointing(NAXIS1, NAXIS2, proc_param.tmp_dir, tanpix, tancoord);
	// cout << tanpix[0] << " " << tanpix[1] << endl;
	// cout << tancoord[0] << " " << tancoord[1] << endl;

	if(rank==0)
		cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl << endl;

	long long test_size;
	read_indpsrc( test_size, npixsrc, indpsrc,  dir.tmp_dir);
	if(test_size != NAXIS1*NAXIS2){
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

	if (pos_param.flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read indpix
	read_indpix(ind_size, npix, indpix,  dir.tmp_dir, flagon);

	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		if(rank==0)
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
		if(rank==0)
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



	if(rank==0)
		printf("\nMain Conjugate gradient loop\n");

	//MALLOC
	S = new double[npix];
	fill(S,S+npix,0.0);


	//	sanepic_conjugate_gradient(samples_struct,pos_param,det,dir,proc_param, npix, S, iframe_min, iframe_max,
	//			fcut,indpix,wcs, NAXIS1, NAXIS2, iterw,
	//			indpsrc, npixsrc,flagon, rank,size, PNdtot, addnpix);

	//////////////////////////////////////////////////////// conjugate gradient


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

	//		int factdupl = 1;
	if(pos_param.flgdupl==1) factdupl = 2;

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

	for (int idupl = 0;idupl<=pos_param.flgdupl;idupl++){



		//Conjugate gradien Inversion
		if (pos_param.projgaps || !flagon){
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

#ifdef LARGE_MEMORY
			fftw_complex  *fdata_buffer;
			fftw_complex *fdata_buffer_tot=NULL;
#endif

			ns = samples_struct.nsamples[iframe];
			f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;


			//    cout << "[" << rank << "] " << iframe << "/" << iframe_max << endl;


			// preconditioner computation : Mp
			if (proc_param.CORRon){
#ifdef LARGE_MEMORY
				// A fdata buffer will be used to avoid binary writing
				fdata_buffer = new fftw_complex[det.ndet*(ns/2+1)];


				if (rank==0){
					fdata_buffer_tot = new fftw_complex[det.ndet*(ns/2+1)];
					for (long ii=0;ii<det.ndet*(ns/2+1);ii++){
						fdata_buffer_tot[ii][0] = 0.0;
						fdata_buffer_tot[ii][1] = 0.0;
					}


				}
				write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size,fdata_buffer);
#else

				write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
				// read pointing + deproject + fourier transform
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
#ifdef LARGE_MEMORY
				MPI_Reduce(fdata_buffer,fdata_buffer_tot,(ns/2+1)*2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
				MPI_Bcast(fdata_buffer,(ns/2+1)*2,MPI_DOUBLE,0,MPI_COMM_WORLD);
				do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
						proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits,name_rank,fdata_buffer);
#else

				do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
						proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank);
				// return Pnd = At N-1 d
#endif

#ifdef LARGE_MEMORY
				delete [] fdata_buffer;
				if(rank==0)
					delete [] fdata_buffer_tot;
#endif
			} else {


				do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
						proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
						iframe, PtNPmatS, Mp, hits,rank,size);
			}

		} // end of iframe loop

		//cout << "do_ptnd" << endl;


#ifdef USE_MPI

		if(rank==0){
			PtNPmatStot = new double[npix];
			hitstot=new long[npix];
			Mptot = new double[npix];
			qtot = new double[npix];

			//				fill(PtNPmatStot,PtNPmatStot+npix,0.0);
			//				fill(hitstot,hitstot+npix,0);
			//				fill(Mptot,Mptot+npix,0.0);
			//				fill(qtot,qtot+npix,0.0);
		}

		MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else

		PtNPmatStot=PtNPmatS; // ajout Mat 02/06
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
		while(((iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !pos_param.flgdupl)) || (!idupl && var_n/var0 > 1e-4)){ // 2000
			// added brackets in order to avoid warning, mat-27/05
			//			if(iter==2)
			//				break;
			fill(q,q+npixeff,0.0); // q <= A*d



			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

#ifdef LARGE_MEMORY
				fftw_complex  *fdata_buffer;
				//if(rank==0)
				fftw_complex *fdata_buffer_tot=NULL;
#endif

				ns = samples_struct.nsamples[iframe];
				//				ff = fframes[iframe];
				f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

				if (proc_param.CORRon){

#ifdef LARGE_MEMORY
					// A fdata buffer will be used to avoid binary writing
					fdata_buffer = new fftw_complex[det.ndet*(ns/2+1)];


					if (rank==0){
						fdata_buffer_tot = new fftw_complex[det.ndet*(ns/2+1)];
						for (long ii=0;ii<det.ndet*(ns/2+1);ii++){
							fdata_buffer_tot[ii][0] = 0.0;
							fdata_buffer_tot[ii][1] = 0.0;
						}


					}
					write_tfAS(d,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size, fdata_buffer);
					// read pointing + deproject + fourier transform
#else
					write_tfAS(d,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
					// read pointing + deproject + fourier transform
#endif
					//					do_PtNd(q,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,f_lppix_Nk,
					//							fsamp,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
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
#ifdef LARGE_MEMORY
					MPI_Reduce(fdata_buffer,fdata_buffer_tot,(ns/2+1)*2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					MPI_Bcast(fdata_buffer,(ns/2+1)*2,MPI_DOUBLE,0,MPI_COMM_WORLD);
					do_PtNd(q, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank, fdata_buffer);
#else


					do_PtNd(q, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
#endif
					// return Pnd = At N-1 d

#ifdef LARGE_MEMORY
					delete [] fdata_buffer;
					if(rank==0)
						delete [] fdata_buffer_tot;
#endif
				} else {

					//					do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
					//							f_lppix_Nk,fsamp,flgdupl,factdupl,ns,ndet,/*size_det,rank_det,*/indpix,
					//							NAXIS1, NAXIS2,npix,iframe,q,NULL,NULL);

					do_PtNPS_nocorr(d, samples_struct.noise_table, dir, det,f_lppix_Nk,
							proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
							iframe, q, NULL, NULL,rank,size);
				}
			} // end of iframe loop




#ifdef USE_MPI
			//cout << " q reduction\n";
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			qtot=q; // ajout mat 02/06
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
			MPI_Barrier(MPI_COMM_WORLD);
			//cout << rank << " S bcast\n";
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif


			// every 10 iterations do ....
			if ((iter % 10) == 0){ // if iter is divisible by 10, recompute PtNPmatStot

				fill(PtNPmatS,PtNPmatS+npixeff,0.0);


				for (long iframe=iframe_min;iframe<iframe_max;iframe++){
#ifdef LARGE_MEMORY
					fftw_complex  *fdata_buffer;
					//if(rank==0)
					fftw_complex *fdata_buffer_tot=NULL;
#endif
					ns = samples_struct.nsamples[iframe];
					f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

					if (proc_param.CORRon){

#ifdef LARGE_MEMORY
						// A fdata buffer will be used to avoid binary writing
						fdata_buffer = new fftw_complex[det.ndet*(ns/2+1)];


						if (rank==0){
							fdata_buffer_tot = new fftw_complex[det.ndet*(ns/2+1)];
							for (long ii=0;ii<det.ndet*(ns/2+1);ii++){
								fdata_buffer_tot[ii][0] = 0.0;
								fdata_buffer_tot[ii][1] = 0.0;
							}


						}
						write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size,fdata_buffer);
						// read pointing + deproject + fourier transform
#else
						write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
						// read pointing + deproject + fourier transform
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
#ifdef LARGE_MEMORY
						MPI_Reduce(fdata_buffer,fdata_buffer_tot,(ns/2+1)*2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
						MPI_Bcast(fdata_buffer,(ns/2+1)*2,MPI_DOUBLE,0,MPI_COMM_WORLD);
						do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
								proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL,name_rank, fdata_buffer);
						// return Pnd = At N-1 d
#else
						do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
								proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);
						// return Pnd = At N-1 d
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

				cout << endl;
				cout << "iter = " << iter;
				cout << ", crit  = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << var_n/var0;
				cout << ", crit2 = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << delta_n/delta0;
				cout << "\r " << flush;


				// This part is now done in sanePre
				//					if (iter == 0){
				//						for (long ii=0; ii<NAXIS1; ii++) {
				//							for (long jj=0; jj<NAXIS2; jj++) {
				//								mi = jj*NAXIS1 + ii;
				//								if (indpix[mi] >= 0){
				//									map1d[mi] = Mptot[indpix[mi]];
				//								} else {
				//									map1d[mi] = NAN;
				//								}
				//							}
				//						}
				//
				//
				//						fname = '!' + dir.outdir + "binMap_noisevar.fits"; // write preconditioner
				//						//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
				//						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
				//
				//								for (long ii=0; ii<NAXIS1; ii++) {
				//							for (long jj=0; jj<NAXIS2; jj++) {
				//								mi = jj*NAXIS1 + ii;
				//								map1d[mi] = 0.0;
				//							}
				//						}
				//
				//						if (addnpix){
				//							for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
				//								for (long ii=0; ii<NAXIS1; ii++) {
				//									for (long jj=0; jj<NAXIS2; jj++) {
				//										mi = jj*NAXIS1 + ii;
				//										long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
				//										if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
				//											map1d[mi] += 1.0/Mptot[indpix[ll]];
				//									}
				//								}
				//							}
				//
				//							fname = '!' + dir.outdir + "optimMap_" + "_invnoisevaruncpix.fits";
				//							//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
				//							write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
				//						}
				//
				//					} // end of if (iter == 0)

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

					temp_stream << "!" + dir.outdir + "optimMap_" + "flux" << iter << "b.fits";


					// récupérer une chaîne de caractères
					fname= temp_stream.str();
					// Clear ostringstream buffer
					temp_stream.str("");
					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,(char *) "Iterative Map",0);

					if (pos_param.flgdupl){
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


						temp_stream << "!" + dir.outdir + "optimMap_fluxflags_" << iter << "b.fits";

						// récupérer une chaîne de caractères
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, (char *)"Duplicated temporary map",0);
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

						temp_stream << "!" + dir.outdir + "optimMap_fluxuncpix_" << iter << "b.fits";

						// récupérer une chaîne de caractères
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d, "Flagged pixels temporary map", 0);

					}
				} // end of if (iterw && (iter % iterw) == 0)


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









		if  ((pos_param.projgaps || (pos_param.flgdupl)) && !idupl){


			fill(PNd,PNd+npix,0.0);

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
#ifdef LARGE_MEMORY
				fftw_complex  *fdata_buffer;
				//if(rank==0)
				fftw_complex *fdata_buffer_tot=NULL;
#endif

				ns = samples_struct.nsamples[iframe];
				f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp;
				f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

				if (proc_param.CORRon){

#ifdef LARGE_MEMORY
					// A fdata buffer will be used to avoid binary writing
					fdata_buffer = new fftw_complex[det.ndet*(ns/2+1)];


					if (rank==0){
						fdata_buffer_tot = new fftw_complex[det.ndet*(ns/2+1)];
						for (long ii=0;ii<det.ndet*(ns/2+1);ii++){
							fdata_buffer_tot[ii][0] = 0.0;
							fdata_buffer_tot[ii][1] = 0.0;
						}


					}
					write_ftrProcesdata(S,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
							npixsrc,addnpix,f_lppix,ns,	iframe, rank, size,name_rank,fdata_buffer);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)

#else
					write_ftrProcesdata(S,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
							npixsrc,addnpix,f_lppix,ns,	iframe, rank, size, name_rank);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)
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
#ifdef LARGE_MEMORY
					MPI_Reduce(fdata_buffer,fdata_buffer_tot,(ns/2+1)*2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					MPI_Bcast(fdata_buffer,(ns/2+1)*2,MPI_DOUBLE,0,MPI_COMM_WORLD);
					do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,"fdata_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank, size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank, fdata_buffer);
#else
					do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,"fdata_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank, size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL, name_rank);

					// return Pnd = At N-1 d

#endif
				} else {

					do_PtNd_nocorr(PNd,dir.tmp_dir,proc_param, pos_param, samples_struct,det, f_lppix, f_lppix_Nk,
							addnpix, ns,indpix, indpsrc, NAXIS1, NAXIS2, npix, npixsrc, iframe, S,rank,size);
				}
#ifdef LARGE_MEMORY
				delete [] fdata_buffer;
				if(rank==0)
					delete [] fdata_buffer_tot;
#endif
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

		write_maps_to_disk(S, NAXIS1, NAXIS2, dir.outdir, indpix, indpsrc,
				Mptot, addnpix, npixsrc, factdupl, samples_struct.ntotscan, wcs);
	}// end of rank==0



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
	//		write_info_for_second_part(proc_param.outdir, NAXIS1, NAXIS2, npix,proc_param.pixdeg, tancoord, tanpix, coordsyst, flagon, indpix);
	//	}

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

	//delete [] frames_index;

	//	wcsfree(wcs);
	wcsvfree(&nwcs, &wcs);




#ifdef USE_MPI
	MPI_Finalize();
#endif

	cout << "\nEnd of sanePic" << endl;


	return 0;
}



