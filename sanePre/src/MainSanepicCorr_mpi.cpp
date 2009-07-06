#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <fcntl.h>
#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>

#include "todprocess.h"
#include "map_making.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
//#include "mpi_architecture_builder.h"
#include "parsePre.h"
#include "sanePre_preprocess.h"

#include "boloIO.h"
#include "inline_IO.h"
#include "nrutil.h"

extern "C" {
#include <fftw3.h>
}


#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

template<class T> void vector2array(std::vector<T> l, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}

//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//




int main(int argc, char *argv[])
{



	int size, size_det;
	int rank, rank_det;

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
	cout << "Mpi is not used for this step" << endl;
#endif

	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);


	//bool projgaps = 0; //1: data flagged are put in a single pixel
	//   (assume no signal in this pixel),
	//0: data flagged are not reprojected

	//default value of the data to pointing shift
	int shift_data_to_point = 0;


	//DEFAULT PARAMETERS
	long napod = 0; // number of samples to apodize
	double fsamp = 0.0;// 25.0; // sampling frequency : BLAST Specific
	double errarcsec = 15.0; // rejection criteria : scerr[ii] > errarcsec, sample is rejected
	// source error

	//Parser parameter (Program options)
	long iframe_min, iframe_max;
	bool NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	bool NOFILLGAP = 0; // fill the gap ? default is YES
	bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	bool flgdupl = 0; // 1 if flagged data are put in a separate map
	int factdupl=1;
	int flagon = 0;
	bool CORRon = 1; // correlation included in the analysis (=1), else 0, default 0
	bool remove_polynomia = 1; // remove a polynomia fitted to the data

	// data parameters
	long *fframes  ; // first frames table ff_in list -> fframes
	long *nsamples ; // number of samples table nf_in list -> nsamples


	long ntotscan; // total number of scans
	long ndet; // number of channels
	int nnf; // extentnoiseSp_list number of elements
	int samples_per_frames=20; // default = 1, BLAST = 20


	int nn, npix; // nn = side of the map, npix = number of filled pixels
	long npixsrc = 0;
	long addnpix=0;


	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan
	double f_lp, f_lppix, f_lppix_Nk; // frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples
	f_lp = 0.0; // low pass filter frequency




	//unsigned char/* *flag, *flpoint, *rejectsamp,*/ *mask; // samples flags, pointing flags, rejected samples list
	double *PNd, *PNdtot; //
	long *indpix, *indpsrc; // pixels indices, mask pixels indices


	string field; // actual boloname in the bolo loop
	string *extentnoiseSp_all; // ((list -> string*))
	string bolofield; // bolofield = boloname + bextension
	string flagfield; // flagfield = field+fextension;
	string dirfile; // data directory
	string outdir; // output directory
	string poutdir; // current path (pPath) or output dir (outdir)
	string bextension; // bolometer field extension
	string fextension = "NOFLAG"; // flag field extension
	string pextension; // pointing extension
	string termin; // output file suffix
	string noiseSppreffile; // noise file suffix
	string prefixe; // prefix used for temporary name file creation
	string flpoint_field = "FLPOINTING";
	string fname; // parallel scheme filename

	/* DEFAULT PARAMETERS */
	int coordsyst = 1; /// Default is RA/DEC
	int coordsyst2 = 1; // to check coordsyst between sanePre and sanePos



	/* Parser inputs */
	std::vector<string> bolonames; // bolometer list
	std::vector<string> extentnoiseSP; // noise file prefix
	std::vector<long> fframes_vec, nsamples_vec; // first frame list, number of frames per sample
	std::vector<long> xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)
	std::vector<double> fcut;

	// Processing time estimation
	time_t t2, t3;// t4, t5, dt;



	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		//parse_sanePos_ini_file(argv[1]);
		int parsed=1;
		parsed=parse_sanePre_ini_file(argv[1],shift_data_to_point,napod,fsamp,NOFILLGAP,NORMLIN,remove_polynomia,flgdupl,
				CORRon,ntotscan,ndet,nnf,f_lp,dirfile,outdir,poutdir,bextension,fextension,
				pextension,termin,noiseSppreffile,coordsyst,bolonames,fframes_vec,nsamples_vec,fname, xxi, xxf, yyi, yyf, extentnoiseSP, fcut);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Finalize();
#endif
			exit(1);
		}
	}



	// processing begins here
	t2=time(NULL);

	if (CORRon) printf("[%2.2i] CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);
	if (!CORRon) printf("[%2.2i] NO CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);

	//malloc
	fframes       = new long[ntotscan];
	nsamples      = new long[ntotscan];
	extentnoiseSp_all = new string[ntotscan];



	//convert vector to standard C array to speed up memory accesses
	vector2array(nsamples_vec, nsamples);
	vector2array(fframes_vec,  fframes);
	vector2array(extentnoiseSP,  extentnoiseSp_all);
	//cout << fframes[0] << fframes[1] << fframes[2] << endl;
	//cout << nsamples[0] << nsamples[1] << nsamples[2] << endl;


	// some needed parameters calculation
	if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}

	if (flgdupl) factdupl = 2;
	string scerr_field = "ERR"+pextension; //source error filename

	// projection parameters
	double *tanpix, *tancoord;
	tanpix=new double[2];
	tancoord=new double[2];

	// read nn, coordsyst, tanpix, tancoord
	read_info_pointing(nn, outdir, termin, coordsyst2, tanpix, tancoord);

	cout << "Map size :" << nn << "x" << nn << endl;
	//cout << tanpix[0] << " " << tanpix[1] << endl;
	//cout << tancoord[0] << " " << tancoord[1] << endl;

	//ensure that coordsyst for position calculation and preprocess are the same
	if (coordsyst!=coordsyst2){
		cerr << "Error : coordinates systems must be the same for preprocessing and mapmaking" << endl;
		exit(0);
	}

	//******************************** some preprocess again // compute indpsrc and addnpix ****************/
	if(xxi.size()>0)
		addnpix=Compute_indpsrc_addnpix(nn,ntotscan,xxi,xxf,yyi,yyf,indpsrc,npixsrc);
	else{
		addnpix=0;
		npixsrc=0;
		indpsrc = new long[nn*nn];
		for(long ii=0;ii<nn*nn;ii++)
			indpsrc[ii]=-1;
	}

	//projection vector
	indpix=new long[factdupl*nn*nn+2 + addnpix];

	//read projection vector from a file
	read_indpix(factdupl*nn*nn+2 + addnpix, npix, indpix, termin, outdir, flagon);

	cout << "apres indpix read" << endl;

#ifdef USE_MPI
	/********************* Define parallelization scheme   *******/

	if (rank == 0){

		long *frnum ;
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

	}



	/*
	if (rank == 0){
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		for (ii=0;ii<=ntotscan;ii++) fprintf(fp,"frnum[%ld] = %ld \n",ii,frnum[ii]);
		fclose(fp);


	// write parallel schema in a file
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"parallel_for_Sanepic_",termin.c_str(),".bi");
		if ((fp = fopen(testfile,"w"))!=NULL){
		//fprintf(fp,"%d\n",size);
		fwrite(&size,sizeof(int), 1, fp);
		fwrite(ruleorder,sizeof(long),ntotscan,fp);
		fwrite(frnum,sizeof(long),ntotscan,fp);
		fclose(fp);
		}else{
			cerr << "Error : couldn't open file to write parallel options. Exiting" << endl;
			exit(1);
		}
	}*/


	MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);

	iframe_min = frnum[rank];
	iframe_max = frnum[rank+1];
	rank_det = 0;
	size_det = 1;

	delete [] frnum;

#else
	iframe_min = 0;
	iframe_max = ntotscan;
	rank_det = rank;
	size_det = size;

#endif

	cout << npix << endl;
	//At N-1 D memory allocation
	PNd = new double[npix];
	cout << "apres pnd" << endl;
	PNdtot = new double[npix];
	cout << "apres pndtot" << endl;
	fill(PNd,PNd+npix,0.0);
	cout << "apres fill" << endl;
	fill(PNdtot,PNdtot+npix,0.0);
	cout << "apres fill" << endl;
	//init1D_double(PNd,0,npix,0.0);
	//init1D_double(PNdtot,0,npix,0.0);


	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//
	cout <<  ndet << endl;
	cout << nsamples[0] << endl;


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

	prefixe = "fdata"; // Fourier transform of the data file prefixe

	// loop over the scans
	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

		ns = nsamples[iframe]; // number of samples for this scan
		ff = fframes[iframe]; //first frame of this scan
		f_lppix = f_lp*double(ns)/fsamp; // knee freq of the filter in terms of samples in order to compute fft
		f_lppix_Nk = fcut[iframe]*double(ns)/fsamp; // noise PS threshold freq, in terms of samples

		printf("[%2.2i] iframe : %ld/%ld",rank,iframe+1,iframe_max);
		cout << " ( -f " << ff << " -l " << ff+ns/20 << " )" << endl;

		// if there is correlation between detectors
		if (CORRon){
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
			write_ftrProcesdata(NULL,indpix,indpsrc,nn,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
					poutdir,termin,errarcsec,dirfile,scerr_field,flpoint_field,bolonames, bextension,
					fextension,shift_data_to_point,f_lppix,ff,ns,napod,ndet,NORMLIN,NOFILLGAP, remove_polynomia,iframe);
			// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)



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
			do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
					fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
			// return Pnd = At N-1 d



		} else {

			do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
					scerr_field,flpoint_field,bolonames,bextension,fextension,
					shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
					flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
					nn,npix,npixsrc,NORMLIN,NOFILLGAP,remove_polynomia,iframe,NULL);

		}


	} // end of iframe loop



	printf("[%2.2i] End of Pre-Processing\n",rank);

#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
	for(int ii=0;ii<npix;ii++)
		PNdtot[ii]=PNd[ii];
#endif


/*for (int ii=0;ii<20;ii++)
	cout << PNdtot[ii] << " ";
cout << endl;
for (int ii=0;ii<20;ii++)
	cout << PNd[ii] << " ";
cout << endl;*/

	if (rank == 0){
		// write At N-1 d in a file
		write_PNd(PNd,npix,termin,outdir);
	}


	/*
	 *
	 * Noise Power spectra estimation has been moved to another program
	 *
	 *

	// Noise Power spectra	estimation loop
	double *S;
	S = new double[npix];


	for (ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii];

	if (doInitPS == 1){

		if (MixMatfile != "NOFILE"){

			for (iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = nsamples[iframe];
				ff = fframes[iframe];
				extentnoiseSp = extentnoiseSp_all[iframe];

				// estimate noise power spectra from data
				EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,
						S,MixMatfile,bolonames,dirfile,bextension,fextension, //cextension,
						shift_data_to_point,poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,
						extentnoiseSp,outdir);

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
#endif
		cout << "Exit after first EstimPowerSpectra" << endl;
		exit(0);

	}

	 */


	/* ---------------------------------------------------------------------------------------------*/
	//Processing stops here
	t3=time(NULL);

	//debug : computation time
	cout << "temps : " << t3-t2 << " sec\n";

	// Close MPI process


#ifdef USE_MPI
	MPI_Finalize();
#endif

	// clean up
	delete [] PNd;
	delete [] PNdtot;
	delete [] fframes;
	delete [] nsamples;
	delete [] indpix;
	delete [] indpsrc;
	delete [] extentnoiseSp_all;

	delete [] tanpix;
	delete [] tancoord;

	return 0;
}


//  cout << "[" << rank << "] End of Init Loop" << endl;


//******************************************************************//
//******************************************************************//
//**********************  End of init loop *************************//
//******************************************************************//
//******************************************************************//



