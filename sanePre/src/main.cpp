#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>


#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "parsePre.h"
#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}

// not sure this is needed
#ifdef PARA_BOLO
#define USE_MPI
#endif

//temp
#include <fstream>


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
	//	cout << size << endl;
	//	cout << rank << endl;
	//	cout << "rank " << rank << " size : " << size << endl;

	if(rank==0)
		printf("\nsanepic_preprocess\n");

#else
	size = 1;
	rank = 0;
	printf("\n sanepic_preprocess\n");
	cout << "Mpi is not used for this step" << endl;
#endif


	struct param_process proc_param;
	struct samples samples_struct;
	struct param_positions pos_param;
	struct common dir;
	struct detectors det;


	// default parameters
	int nwcs=1;
	long iframe_min=0, iframe_max=0; /*!  min and max number of frame (used with mpi) */
	int factdupl=1; /*! map duplication factor */
	int flagon = 0; /*! if a data is flagged */




	long NAXIS1, NAXIS2;
	long long npix; /*! npix = number of filled pixels*/
	long long npixsrc = 0; /*! number of pixels contained in box constraint removal */
	long long addnpix=0; /* number of pixels to add to the final map */
	long long ind_size; /* indpix readed size */

	//internal data params
	long ns; /*! number of samples for this scan, first frame number of this scan*/
	double f_lppix, f_lppix_Nk; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples*/


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

	int parsed=0;


	if (argc<2)
		parsed=1;
	else {
		// Parse ini file
		parsed=parse_sanePre_ini_file(argv[1],proc_param, pos_param, dir, samples_struct,
				det, fcut, rank, size);
	}



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

	if ((parsed>0)||(!compute_dirfile_format_fdata(dir.tmp_dir, det, samples_struct.ntotscan, rank))){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}


#ifdef DEBUG
	std::ostringstream oss;
	string name_rank;
	oss << dir.output_dir + "debug_sanePre_" << rank << ".txt";
	name_rank = oss.str();
#else
	string name_rank = dir.output_dir + "debug_sanePre.txt";

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

	// processing begins here
	t2=time(NULL);

	//	long *frames_index;
	//
	//	frames_index = new long [samples_struct.ntotscan];


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new int[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	if (pos_param.flgdupl) factdupl = 2;// -M =1, default 0 : if flagged data are put in a duplicated map

	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs,&NAXIS1, &NAXIS2);

	//	cout << "nwcs : " << nwcs << endl;

	if(rank==0)
		cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;


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


	//read projection vector from a file
	read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon);

	// Check indpix readed size = expected size
	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		if(rank==0){
			cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
			cout << ind_size << " != "  << (factdupl*NAXIS1*NAXIS2+2 + addnpix) << " " << factdupl << " " << addnpix << endl;
		}
		exit(0);
	}

	//TODO: check the noise matrix file at the beginning of each frame/start of the program...


#if defined(USE_MPI) && ! defined(PARA_BOLO)

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

	//	cout <<" final list : " << endl;
	//	for(int ii = 0; ii< samples_struct.ntotscan;ii++){
	//		cout << samples_struct.fits_table[ii] << " " << samples_struct.noise_table[ii] << " " << samples_struct.index_table[ii] << endl;
	//		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	//}
#else
	fname = dir.output_dir + parallel_scheme_filename;
	int test=0;
	test=check_ParallelizationScheme(fname,dir.dirfile,samples_struct,size);
	if (test==-1){
		if(rank==0)
			cerr << "erreur dans check_parallelizationScheme non-MPI " << endl;
		exit(0);
	}

	//	cout <<" final list : " << endl;
	for(int ii = 0; ii< samples_struct.ntotscan;ii++){
		//		cout << samples_struct.fits_table[ii] << " " << samples_struct.noise_table[ii] << " " << samples_struct.index_table[ii] << endl;
		samples_struct.fits_table[ii]=dir.dirfile + samples_struct.fits_table[ii];
	}

#endif
	// add reading parallel scheme procedure !

	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;


#endif


	//	cout << "rank " << rank << " frame : " << iframe_min << " " << iframe_max << endl;
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
	if(rank==0)
		printf("\nPre-processing of the data\n");

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

#ifdef LARGE_MEMORY
		fftw_complex  *fdata_buffer;
		//if(rank==0)
		fftw_complex *fdata_buffer_tot=NULL;
#endif

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

			write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,rank,size,name_rank,fdata_buffer);
#else

#if defined(USE_MPI) && ! defined(PARA_BOLO)
			write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,0,1, name_rank);
#else
			write_ftrProcesdata(NULL,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe,rank,size, name_rank);
#endif
#endif
			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)

			//			cout << "avant time ! \n";
			//Processing stops here
			t3=time(NULL);

			//debug : computation time
			//			if(rank==0)
			//				cout << " [ " << rank << " ] temps : " << t3-t2 << " sec\n";
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
#ifdef LARGE_MEMORY
			do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank,fdata_buffer);
			// Returns Pnd = (At N-1 d)
#else

#if defined(USE_MPI) && ! defined(PARA_BOLO)
			do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,0,1,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank);
#else
			do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
					proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits, name_rank);

#endif
#endif
			// delete fdata buffer
			//delete [] fdata_buffer;

			/*do_PtNd2(PNd,NULL,extentnoiseSp_all,noiseSppreffile,outdir,prefixe,termin_internal,bolonames,f_lppix_Nk,
								fsamp,ff,ns,ndet,size_det,rank_det,indpix,indpsrc,npixsrc,ntotscan,addnpix,flgdupl,factdupl,
								2,errarcsec,dirfile,scerr_field,flpoint_field,bextension,fextension,shift_data_to_point,
								napod,NORMLIN,NOFILLGAP,remove_polynomia,nn,npix,iframe,NULL,NULL);*/
#ifdef LARGE_MEMORY
			delete [] fdata_buffer;
			if(rank==0)
				delete [] fdata_buffer_tot;
#endif
		} else {


			do_PtNd_nocorr(PNd, dir.tmp_dir,proc_param,pos_param,samples_struct,
					det,f_lppix,f_lppix_Nk,addnpix,
					ns/*,size_det,rank_det*/,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,iframe,NULL,rank,size);
			// fillgaps + butterworth filter + fourier transform and PNd generation

		}


	} // end of iframe loop



#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

#else
	hitstot=hits;
	PNdtot=PNd;
	Mptot=Mp;

#endif

	if(rank==0)
		printf("\nEnd of Pre-Processing\n");

	if (rank == 0){
		// write (At N-1 d) in a file
		write_PNd(PNdtot,npix,dir.tmp_dir);

		//temp
		//		ofstream filee;
		//		string outfile = dir.output_dir + "test_Pnd_Mp.txt";
		//		filee.open(outfile.c_str(), ios::out);
		//		if(!filee.is_open()){
		//			cerr << "File [" << fname << "] Invalid." << endl;
		//			exit(0);
		//		}

		if(rank==0)
			cout << "\nNaive step :" << endl;
		string fnaivname;
		double *map1d;
		long long mi;
		map1d = new double[NAXIS1*NAXIS2];

		//TODO:   addnpix pixel
		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					if(PNdtot[indpix[mi]]==NAN) // replace NAN by 0
						map1d[mi] = 0;
					else
						map1d[mi] = PNdtot[indpix[mi]]/Mptot[indpix[mi]];
					//					filee << PNdtot[indpix[mi]] << " " << Mptot[indpix[mi]] << endl;
					//					getchar();
				} else {
					map1d[mi] = NAN; // replace map[hits[mi] = 0] = NAN
				}
			}
		}

		//		filee.close();
		fnaivname = '!' + dir.output_dir + "naivMap.fits";
		cout << "Output file : " << fnaivname << endl;
		//write_fits(fnaivname, 0, NAXIS1, NAXIS2, tanpix, tancoord, 1, 'd', (void *)map1d);
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Image",0);

		// TODO: Save the map of flag data if needed

		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = PNdtot[indpix[mi]]/hitstot[indpix[mi]];
				} else {
					map1d[mi] = 0;
				}
			}
		}

		fnaivname = dir.output_dir + "naivMap.fits";
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Ultra Naiv",1);

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

		//fnaivname = dir.outdir + "naivMap.fits";
		//		fnaivname = '!' + dir.outdir + "hits.fits";
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Coverage",1);


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


		//		fnaivname = '!' + dir.output_dir + "binMap_noisevar.fits"; // write preconditioner
		//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Error",1);

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

			//			fnaivname = '!' + dir.output_dir + "optimMap_" + "_invnoisevaruncpix.fits";
			//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
			write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d,"Invnoisevaruncpix",1);
		}


		delete [] map1d;
		if(rank==0)
			printf("End of saneNaiv\n\n");

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
	cout << "Total Time : " << t3-t2 << " sec\n";
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
	MPI_Finalize();
#endif

	if(rank==0)
		printf("\nEnd of sanePre \n");

	return 0;
}

//******************************************************************//
//******************************************************************//
//**********************  End of init loop *************************//
//******************************************************************//
//******************************************************************//



