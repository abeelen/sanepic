#include <iostream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>


#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "parsePre.h"
#include "sanePre_preprocess.h"
#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}

#ifdef PARA_BOLO
#define USE_MPI
#endif

//temp
#include <fstream>


#ifdef PARA_BOLO
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
/*! Sanepic preprocess main function */
{




	int size;/*!< number of processors */
	int rank;

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
	cout << "Mpi is not used for this step" << endl;
#endif


	struct param_process proc_param;
	struct samples samples_struct;
	struct param_positions pos_param;
	struct directories dir;
	struct detectors det;


	//DEFAULT PARAMETERS
	proc_param.napod = 0; /*! number of samples to apodize*/
	proc_param.fsamp = 0.0;// 25.0; /*! sampling frequency : BLAST Specific*/

	int nwcs=1;

	//Parser parameter (Program options)
	long iframe_min=0, iframe_max=0; /*!  min and max number of frame (used with mpi) */
	proc_param.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	proc_param.NOFILLGAP = 0; /*! fill the gap ? default is YES*/
	int factdupl=1; /*! map duplication factor */
	int flagon = 0; /*! if a data is flagged */
	proc_param.CORRon = 1; /*! correlation included in the analysis (=1), else 0, default 0*/
	proc_param.remove_polynomia = 1; /*! remove a polynomia fitted to the data*/


	samples_struct.ntotscan=0; /*! total number of scans*/
	det.ndet=0; /*! number of channels*/


	long NAXIS1, NAXIS2;
	long long npix; /*! npix = number of filled pixels*/
	long long npixsrc = 0; /*! number of pixels contained in box constraint removal */
	long long addnpix=0; /* number of pixels to add to the final map */
	long long ind_size; /* indpix readed size */

	//internal data params
	long ns; /*! number of samples for this scan, first frame number of this scan*/
	double f_lppix, f_lppix_Nk; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples*/
	proc_param.f_lp = 0.0; // low pass filter frequency

	//fftw_complex *fdata_buffer; /*! buffer used to store all the fdata arrays instead of writing on disk */


	double *PNd, *PNdtot=NULL; /*!  projected noised data, and global Pnd for mpi utilization */
	double *Mp, *Mptot=NULL;
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices*/
	long *hits, *hitstot=NULL; /*! naivmap parameters : hits count */

	string field; /*! actual boloname in the bolo loop*/
	string fname;
	string prefixe; /*! prefix used for temporary name file creation*/


	/* Parser inputs */
	std::vector<double> fcut; /* noise cutting frequency */

	// Processing time estimation
	time_t t2, t3;// t4, t5, dt;

	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		int parsed=1;
		parsed=parse_sanePre_ini_file(argv[1],proc_param, pos_param, dir, samples_struct,
				det, fcut, rank);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}
	}



	// processing begins here
	t2=time(NULL);
	//
	//	long *frames_index;
	//
	//	frames_index = new long [samples_struct.ntotscan];


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	if (pos_param.flgdupl) factdupl = 2;

	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs,&NAXIS1, &NAXIS2);

	//	cout << "nwcs : " << nwcs << endl;

	if(rank==0)
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

	if (pos_param.flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	//read projection vector from a file
	read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon);

	// Check indpix readed size = expected size
	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}

	//TODO: check the noise matrix file at the beginning of each frame/start of the program...


#if defined(USE_MPI) && ! defined(PARA_BOLO)

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

			string outfile = dir.outdir + samples_struct.filename + "_sanepre.txt";
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
						cout << temp.substr(found+1) << " " << samples_struct.noise_table[num_frame] << " " << ii << endl;
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
		//		cout << "on aura : \n";
		//		cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
		//		cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
		//		cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;
		//		//cout << samples_struct.filename << endl;
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

	//	MPI_Barrier(MPI_COMM_WORLD);
#else
#if defined(USE_MPI) && defined(PARA_BOLO)
	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);

	cout <<" final list : " << endl;
	for(int ii = 0; ii< samples_struct.ntotscan;ii++){
		cout << samples_struct.fits_table[ii] << " " << samples_struct.noise_table[ii] << " " << samples_struct.index_table[ii] << endl;
		//		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	}
#else
	fname = dir.outdir + parallel_scheme_filename;
	int test=0;
	test=check_ParallelizationScheme(fname,dir.dirfile,samples_struct,size);
	if (test==-1){
		cerr << "erreur dans check_parallelizationScheme non-MPI " << endl;
		exit(0);
	}

	cout <<" final list : " << endl;
	for(int ii = 0; ii< samples_struct.ntotscan;ii++){
		cout << samples_struct.fits_table[ii] << " " << samples_struct.noise_table[ii] << " " << samples_struct.index_table[ii] << endl;
		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	}

#endif
	// add reading parallel scheme procedure !

	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	//convert vector to standard C array to speed up memory accesses
	//	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	//	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	//	vector2array(samples_struct.scans_index,  samples_struct.index_table);

	//	for(long ii=0; ii<samples_struct.ntotscan;ii++)
	//		frames_index[ii] = ii;


#endif

	//exit(0);

	//At N-1 D memory allocation
	PNd = new double[npix];

	// global At N-1 D malloc for mpi
	hits=new long[npix];
	Mp = new double[npix];

	// initialisation to 0.0
	fill(PNd,PNd+npix,0.0);
	fill(hits,hits+npix,0);
	fill(Mp,Mp+npix,0.0);


#ifdef USE_MPI

	if(rank==0){
		PNdtot = new double[npix];
		hitstot=new long[npix];
		Mptot = new double[npix];

		fill(PNdtot,PNdtot+npix,0.0);
		fill(hitstot,hitstot+npix,0);
		fill(Mptot,Mptot+npix,0.0);
	}

#endif


	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//
	if(iframe_min!=iframe_max)
		printf("[%2.2i] Pre-processing of the data\n",rank);

	/*fftw_complex **fdatas;
	long nsamp_max=0;
	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		ns = nsamples[iframe];
		if(ns>nsamp_max)
			nsamp_max=ns;
	}*/

	/*fdatas=new fftw_complex*[ndet];
	for (long ii=0;ii<ndet;ii++)
		fdatas[ii]=new fftw_complex[nsamp_max/2+1];*/

	prefixe = "fdata_"; // Fourier transform of the data file prefixe

	// loop over the scans
	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

		ns = samples_struct.nsamples[iframe]; // number of samples for this scan
		f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp; // knee freq of the filter in terms of samples in order to compute fft
		f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp; // noise PS threshold freq, in terms of samples

		if(iframe_min!=iframe_max)
			printf("[%2.2i] iframe : %ld/%ld\n",rank,iframe+1,iframe_max);

		// if there is correlation between detectors
		if (proc_param.CORRon){
			// var double *S =  NULL
			// indpix = pixel indice
			// indpsrc = box crossing constraint removal pixel indice
			// nn = size of map
			// npix = number of pixels that are seen
			// npixsrc = number of pixels in box CCRemoval
			// ntotscan = total number of scans
			// addnpix = number of added pixels in the map
			// flgdupl = flaggued pixels are in a duplicate map : 1/0
			// factdupl = duplication de la map : 2/1
			// fillg =2 ????
			// poutdir = outpout dir or current path (default)
			// termin = output file suffix
			// errarcsec = pointing error threshold
			// dirfile = data directory
			// scerr_field = "ERR" + pextension (_def for example)
			// flpoint_field = "FLPOINTING"
			// bolonames = bolo names array
			// bextension = -B option : "_data" for example
			// fextension = "NOFLAG" or -G option ("_flag" for example)
			// cextension = "NOCALP" or -R option ("_calp" for example)
			// shift_data_to_point (default 0), for subtracting a time offset to the data to match the pointing
			// f_lppix = filter freq in term of sample
			// ff = first frame number of this scan
			// ns = number of sample for this scan
			// napod = number of samples to apodize -A option
			// ndet = bolo total number
			// NORMLIN = baseline is remove from the data, default =0, option -L
			// NOFILLGAP = fill the gap ? default yes => 0
			// iframe = scan number : 0=> ntotscan

			// A fdata buffer will be used to avoid binary writing
			//fdata_buffer = new fftw_complex[ndet*(ns/2+1)];

			write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,rank,size);

			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)

			//			cout << "avant time ! \n";
			//Processing stops here
			t3=time(NULL);

			//debug : computation time
			if(iframe_min!=iframe_max)
				cout << " [ " << rank << " ] temps : " << t3-t2 << " sec\n";
#ifdef PARA_BOLO
			//cout << "rank " << rank << " a fini et attend ! \n";
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			// PNd = npix dimension, initialised to 0.0
			// extentnoiseSp_all = list of power spectrum file names (for each scan or same PS for all the scans)
			// noiseSppreffile = noise power spectrum file suffix = path
			// poutdir = outpout dir or current path (default)
			// prefixe = "fdata"; => prefixe de lecture/sauvegarde des données
			// termin = output file suffix
			// bolonames = bolo names array
			// f_lppix_Nk = freq threshold noise en terme de sample
			// fsamp = freq echantillonage des data
			// ff = n° premier sample du scan
			// ns = nombre de sample ds le scan
			// ndet = nombre de bolo
			// size = 1 // cf mpi
			// rank = 0 // cf mpi
			// indpix = pixel indice double[nn*nn]
			// nn = taille de la carte (1 coté)
			// npix = total number of filled pixels (pixel dont on a les data correspondantes)
			// iframe = indice du scan
			// *Mp = Null : la map ???
			// *Hits = Null
			do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits/*,fdata_buffer*/);
			// Returns Pnd = (At N-1 d)

			// delete fdata buffer
			//delete [] fdata_buffer;

			/*do_PtNd2(PNd,NULL,extentnoiseSp_all,noiseSppreffile,outdir,prefixe,termin_internal,bolonames,f_lppix_Nk,
								fsamp,ff,ns,ndet,size_det,rank_det,indpix,indpsrc,npixsrc,ntotscan,addnpix,flgdupl,factdupl,
								2,errarcsec,dirfile,scerr_field,flpoint_field,bextension,fextension,shift_data_to_point,
								napod,NORMLIN,NOFILLGAP,remove_polynomia,nn,npix,iframe,NULL,NULL);*/

		} else {


			do_PtNd_nocorr(PNd, dir.tmp_dir,proc_param,pos_param,samples_struct,
					det,f_lppix,f_lppix_Nk,addnpix,
					ns/*,size_det,rank_det*/,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,iframe,NULL);
			// fillgaps + butterworth filter + fourier transform and PNd generation

		}


	} // end of iframe loop




#ifdef USE_MPI
	if(iframe_min!=iframe_max)
		printf("[%2.2i] End of Pre-Processing\n",rank);

	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

#else
	printf("End of Pre-Processing\n");

	hitstot=hits;
	PNdtot=PNd;
	Mptot=Mp;

#endif



	if (rank == 0){
		// write (At N-1 d) in a file
		write_PNd(PNdtot,npix,dir.tmp_dir);

		//temp
		//		ofstream filee;
		//		string outfile = dir.outdir + "test_Pnd_Mp.txt";
		//		filee.open(outfile.c_str(), ios::out);
		//		if(!filee.is_open()){
		//			cerr << "File [" << fname << "] Invalid." << endl;
		//			exit(0);
		//		}


		cout << "naive step" << endl;
		string fnaivname;
		double *map1d;
		long long mi;
		map1d = new double[NAXIS1*NAXIS2];

		//TODO: replace the NAN by 0 + addnpix pixel + replace map[hits[mi] = 0] = NAN
		//TODO: save everythings in a single fits file
		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = PNdtot[indpix[mi]]/Mptot[indpix[mi]];
					//					filee << PNdtot[indpix[mi]] << " " << Mptot[indpix[mi]] << endl;
					//					getchar();
				} else {
					map1d[mi] = NAN;
				}
			}
		}

		//		filee.close();
		fnaivname = '!' + dir.outdir + "naivMap.fits";
		cout << fnaivname << endl;
		//write_fits(fnaivname, 0, NAXIS1, NAXIS2, tanpix, tancoord, 1, 'd', (void *)map1d);
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = hitstot[indpix[mi]];
				} else {
					map1d[mi] = 0;
				}
			}
		}

		if (addnpix){
			for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++){
				for (long jj=0; jj<NAXIS2; jj++) {
					for (long ii=0; ii<NAXIS1; ii++) {
						mi = jj*NAXIS1 + ii;
						long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
						if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
							map1d[mi] += hitstot[indpix[ll]];
					}
				}
			}
		}

		fnaivname = '!' + dir.outdir + "hits.fits";
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


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


		fnaivname = '!' + dir.outdir + "binMap_noisevar.fits"; // write preconditioner
		//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

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

			fnaivname = '!' + dir.outdir + "optimMap_" + "_invnoisevaruncpix.fits";
			//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
			write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
		}


		delete [] map1d;
		printf("End of saneNaiv\n");

	}
	/* ---------------------------------------------------------------------------------------------*/

	//Processing stops here
	t3=time(NULL);

	//debug : computation time
#ifdef USE_MPI
	if(iframe_min!=iframe_max)
		printf("[%2.2i] Time : %d sec\n",rank, (int)(t3-t2));

	if(rank==0){
		// clean up
		delete [] PNdtot;
		delete [] Mptot;
		delete [] hitstot;

	}
#else
	cout << "Time : " << t3-t2 << " sec\n";
#endif

	// clean up
	delete [] PNd;
	delete [] Mp;
	delete [] hits;


	delete [] samples_struct.nsamples;
	delete [] indpix;
	delete [] indpsrc;

	delete [] samples_struct.noise_table;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;

	//	wcsfree(wcs);
	wcsvfree(&nwcs, &wcs);


#ifdef USE_MPI
	printf("[%2.2i] End of sanePre \n",rank);
	MPI_Finalize();
#else
	printf("End of sanePre \n");
#endif


	return 0;
}

//******************************************************************//
//******************************************************************//
//**********************  End of init loop *************************//
//******************************************************************//
//******************************************************************//



