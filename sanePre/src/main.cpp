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
}




#ifdef USE_MPI
#include "mpi.h"
//#include "mpi_architecture_builder.h"
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


	struct user_options u_opt;
	struct samples samples_struct;
	struct input_commons com;
	struct directories dir;
	struct detectors det;


	//DEFAULT PARAMETERS
	com.napod = 0; /*! number of samples to apodize*/
	u_opt.fsamp = 0.0;// 25.0; /*! sampling frequency : BLAST Specific*/



	//Parser parameter (Program options)
	long iframe_min=0, iframe_max=0; /*!  min and max number of frame (used with mpi) */
	u_opt.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	com.NOFILLGAP = 0; /*! fill the gap ? default is YES*/
	com.flgdupl = 0; /*! 1 if flagged data are put in a separate map*/
	int factdupl=1; /*! map duplication factor */
	int flagon = 0; /*! if a data is flagged */
	u_opt.CORRon = 1; /*! correlation included in the analysis (=1), else 0, default 0*/
	u_opt.remove_polynomia = 1; /*! remove a polynomia fitted to the data*/


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
	u_opt.f_lp = 0.0; // low pass filter frequency

	//fftw_complex *fdata_buffer; /*! buffer used to store all the fdata arrays instead of writing on disk */


	double *PNd, *PNdtot; /*!  projected noised data, and global Pnd for mpi utilization */
	double *Mp,*Mptot;
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices*/
	unsigned short *mask;
	long *hits,*hitstot; /*! naivmap parameters : hits count */

	string field; /*! actual boloname in the bolo loop*/
	string fname;
	string prefixe; /*! prefix used for temporary name file creation*/



	/* Parser inputs */
	std::vector<struct box> boxFile;
	std::vector<double> fcut; /* noise cutting frequency */

	// Processing time estimation
	time_t t2, t3;// t4, t5, dt;

	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		int parsed=1;
		parsed=parse_sanePre_ini_file(argv[1],u_opt, dir, samples_struct,com,
				det,boxFile, fcut);

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

	long *frames_index;

	frames_index = new long [samples_struct.ntotscan];


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	if (com.flgdupl) factdupl = 2;

	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs, &NAXIS1, &NAXIS2);

	cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;


	//************************************* Deal with masking the point sources
	mask    = new unsigned short[NAXIS1*NAXIS2];
	indpsrc = new long long[NAXIS1*NAXIS2];

	// Initialize the masks
	addnpix=0;
	npixsrc=0;
	for (long ii=0; ii<NAXIS1*NAXIS2; ii++){
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

	// TODO: Is it really needed here ?
	//projection vector
	indpix=new long long[factdupl*NAXIS1*NAXIS2+2 + addnpix];

	//read projection vector from a file
	read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon);

	// Check indpix readed size = expected size
	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}



#ifdef USE_MPI

	int test=0;
	fname = dir.outdir + parallel_scheme_filename;
	cout << fname << endl;
	//test=define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct.ntotscan,size,samples_struct.nsamples,samples_struct.fitsvect,samples_struct.noisevect,samples_struct.fits_table, samples_struct.noise_table,samples_struct.index_table);
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


	//At N-1 D memory allocation
	PNd = new double[npix];

	// global At N-1 D malloc for mpi

	PNdtot = new double[npix];
	hits=new long[npix];
	hitstot=new long[npix];

	Mp = new double[npix];
	Mptot = new double[npix];



	// initialisation to 0.0
	fill(PNd,PNd+npix,0.0);
	fill(PNdtot,PNdtot+npix,0.0);
	fill(hits,hits+npix,0);
	fill(hitstot,hitstot+npix,0);
	fill(Mp,Mp+npix,0.0);
	fill(Mptot,Mptot+npix,0.0);



	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//

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
		//		ff = fframes[iframe]; //first frame of this scan
		f_lppix = u_opt.f_lp*double(ns)/u_opt.fsamp; // knee freq of the filter in terms of samples in order to compute fft
		f_lppix_Nk = fcut[iframe]*double(ns)/u_opt.fsamp; // noise PS threshold freq, in terms of samples

		printf("[%2.2i] iframe : %ld/%ld\n",rank,iframe+1,iframe_max);


		// if there is correlation between detectors
		if (u_opt.CORRon){
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

			//TODO : WHY fillg is hardcoded to 2 ????
			write_ftrProcesdata(NULL,u_opt,samples_struct,com,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
					npixsrc,addnpix,f_lppix,ns,	iframe);

			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)

			//Processing stops here
			t3=time(NULL);

			//debug : computation time
			cout << "temps : " << t3-t2 << " sec\n";

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
					u_opt.fsamp,ns/*,size_det,rank_det*/,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits/*,fdata_buffer*/);
			// Returns Pnd = (At N-1 d)

			// delete fdata buffer
			//delete [] fdata_buffer;

			/*do_PtNd2(PNd,NULL,extentnoiseSp_all,noiseSppreffile,outdir,prefixe,termin_internal,bolonames,f_lppix_Nk,
								fsamp,ff,ns,ndet,size_det,rank_det,indpix,indpsrc,npixsrc,ntotscan,addnpix,flgdupl,factdupl,
								2,errarcsec,dirfile,scerr_field,flpoint_field,bextension,fextension,shift_data_to_point,
								napod,NORMLIN,NOFILLGAP,remove_polynomia,nn,npix,iframe,NULL,NULL);*/

		} else {


			do_PtNd_nocorr(PNd, dir.tmp_dir,u_opt,samples_struct,com,
					det,f_lppix,f_lppix_Nk,addnpix,
					ns/*,size_det,rank_det*/,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,iframe,NULL);
			// fillgaps + butterworth filter + fourier transform and PNd generation

		}


	} // end of iframe loop




	printf("[%2.2i] End of Pre-Processing\n",rank);

#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

#else
	for(unsigned long ii=0;ii<npix;ii++){
		hitstot[ii]=hits[ii];
		PNdtot[ii]=PNd[ii]; // fill PNdtot with PNd in case mpi is not used
		Mptot[ii]=Mp[ii];
	}
#endif



	if (rank == 0){
		// write (At N-1 d) in a file
		write_PNd(PNdtot,npix,/*termin_internal,*/dir.tmp_dir);


		cout << "naive step" << endl;
		string fnaivname;
		double *map1d;
		long long mi;
		map1d = new double[NAXIS1*NAXIS2];


		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					//if(hits[indpix[mi]]>0)
					map1d[mi] = PNdtot[indpix[mi]]/Mptot[indpix[mi]];
				} else {
					map1d[mi] = NAN;
				}
			}
		}

		fnaivname = '!' + dir.outdir + "naivMap.fits";
		cout << fnaivname << endl;
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = hitstot[indpix[mi]];
				} else {
					map1d[mi] = NAN;
				}
			}
		}

		fnaivname = '!' + dir.outdir + "naivMaphits.fits";
		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


		printf("End of saneNaiv\n");

	}
	/* ---------------------------------------------------------------------------------------------*/

	//Processing stops here
	t3=time(NULL);
	//debug : computation time
	cout << "temps : " << t3-t2 << " sec\n";


	// clean up
	delete [] PNd;
	delete [] PNdtot;

	delete [] Mp;
	delete [] Mptot;

	delete [] hits;
	delete [] hitstot;


	delete [] samples_struct.nsamples;
	delete [] indpix;
	delete [] indpsrc;

	delete [] samples_struct.noise_table;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] frames_index;


	cout << "[" << rank << "] End of Init Loop" << endl;

#ifdef USE_MPI
	MPI_Finalize();
#endif


	return 0;
}

//******************************************************************//
//******************************************************************//
//**********************  End of init loop *************************//
//******************************************************************//
//******************************************************************//



