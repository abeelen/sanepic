#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <fcntl.h>
#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <algorithm>

#include "todprocess.h"
#include "map_making.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "mpi_architecture_builder.h"
#include "parsePre.h"

#include "positionsIO.h"
#include "imageIO.h"
#include "dataIO.h"

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

	long iframe, ll, ib, ii; // loop indices
	long iframe_min, iframe_max;
	bool NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	bool NOFILLGAP = 0; // fill the gap ? default is YES
	//bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	bool flgdupl = 0; // 1 if flagged data are put in a separate map
	int factdupl=1;
	int flagon = 0;
	bool CORRon = 1; // correlation included in the analysis (=1), else 0, default 0



	// data parameters
	long *fframes  ; // first frames table ff_in list -> fframes
	long *fframesorder ;
	long *nsamples ; // number of samples table nf_in list -> nsamples
	long *nsamplesorder ;
	long *ruleorder ;
	long *frnum ;
	long *hits;


	long ntotscan; // total number of scans
	long ndet; // number of channels
	int nnf; // extentnoiseSp_list number of elements
	int samples_per_frames=20; // default = 1, BLAST = 20



	// map making parameters
	//double pixdeg = 0.00168725828819; // size of pixels (degree)
	double pixdeg;

	int nn, npix; // nn = side of the map, npix = number of filled pixels
	double *tancoord; // tangent point coordinates
	double *tanpix; // tangent pixel
	long npixsrc = 0;
	long addnpix=0;


	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan
	double f_lp, f_lp_Nk, f_lppix, f_lppix_Nk; // frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples


	FILE *fp;

	char testfile[100];
	string fname;

	unsigned char/* *flag, *flpoint, *rejectsamp,*/ *mask; // samples flags, pointing flags, rejected samples list
	double *PNd, *PNdtot; //
	long *indpix, *indpsrc; // pixels indices, mask pixels indices


	string field; // actual boloname in the bolo loop
	// moved array of strings to vector of strings...
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
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation
	string flpoint_field = "FLPOINTING";


	/* DEFAULT PARAMETERS */
	int coordsyst = 1; /// Default is RA/DEC
	int coordsyst2 = 1;



	std::vector<string> bolonames; // bolometer list
	std::vector<long> fframes_vec, nsamples_vec; // first frame list, number of frames per sample
	std::vector<long> xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)





	f_lp = 0.0; // low pass filter frequency
	f_lp_Nk = 0.0; // noise PS frequency threshold



	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		//parse_sanePos_ini_file(argv[1]);
		parse_sanePre_ini_file(argv[1],shift_data_to_point,napod,fsamp,NOFILLGAP,NORMLIN,flgdupl,
				CORRon,ntotscan,ndet,nnf,f_lp,f_lp_Nk,dirfile,outdir,poutdir,bextension,fextension,
				pextension,termin,noiseSppreffile,coordsyst,bolonames,fframes_vec,nsamples_vec, fname, pixdeg);

	}




	fframes       = new long[ntotscan];
	fframesorder  = new long[ntotscan];
	nsamples      = new long[ntotscan];
	nsamplesorder = new long[ntotscan];



	vector2array(nsamples_vec, nsamples);
	vector2array(fframes_vec,  fframes);
	//cout << fframes[0] << fframes[1] << fframes[2] << endl;
	//cout << nsamples[0] << nsamples[1] << nsamples[2] << endl;

	string scerr_field = "ERR"+pextension;

	if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}

	if (flgdupl) factdupl = 2;


	// read nn, coordsyst, tanpix, tancoord
	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"InfoPointing_for_Sanepic_",termin.c_str(),".txt");
	if ((fp = fopen(testfile,"r")) == NULL){
		cerr << "File InfoPointing_for_sanepic :" << testfile << "not found. Exiting" << endl;
		exit(1);
	}

	tanpix=new double[2];
	tancoord=new double[2];
	fscanf(fp,"%d\n",&nn);
	fscanf(fp,"%d\n",&coordsyst2);
	fscanf(fp,"%lf\n",tanpix);
	fscanf(fp,"%lf\n",tanpix+1);
	fscanf(fp,"%lf\n",tancoord);
	fscanf(fp,"%lf\n",tancoord+1);
	fclose(fp);

	cout << "Map size :" << nn << "x" << nn << endl;
	//cout << tanpix[0] << " " << tanpix[1] << endl;
	//cout << tancoord[0] << " " << tancoord[1] << endl;

	if (coordsyst!=coordsyst2){
		cerr << "Error : coordinates systems must be the same for preprocessing and mapmaking" << endl;
		exit(0);
	}

	//******************************** some preprocess again // compute indpsrc and addnpix ****************/

	//************************************* Deal with masking the point sources
	// define the mask
	mask = new unsigned char[nn*nn];
	for (ii=0;ii<nn*nn;ii++)
		mask[ii] = 1;


	if (xxi.size() != 0){
		for (ib = 0;ib < (long)xxi.size(); ib++){ // to avoid warning, mat-27/05
			// for each box crossing constraint removal
			for (ii=xxi[ib];ii<xxf[ib];ii++)
				for (ll=yyi[ib];ll<yyf[ib];ll++)
					mask[ll*nn + ii] = 0;  // mask is initialised to 0
		}
	}



	//long npixsrc = 0;
	indpsrc = new long[nn*nn];
	for (ii=0;ii<nn*nn;ii++){
		if (mask[ii] == 0){
			indpsrc[ii] = npixsrc;
			npixsrc += 1;
		} else {
			indpsrc[ii] = -1;
		}
	}
	addnpix = ntotscan*npixsrc; // addnpix = number of pix to add in pixon = number of scans * number of pix in box crossing constraint removal

	cout  << "addnpix : " << addnpix << endl;

	indpix=new long[factdupl*nn*nn+2 + addnpix];

	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"Indpix_for_conj_grad_",termin.c_str(),".bi");
	cout << testfile << endl;
	if ((fp = fopen(testfile,"r"))!=NULL){
		fread(&flagon,sizeof(int),1,fp);
		fread(&npix,sizeof(int),1,fp);
		fread(indpix,sizeof(long),factdupl*nn*nn+2 + addnpix,fp);
		/*fscanf(fp,"%d",&flagon);
			for(int ii=0;ii<factdupl*nn*nn+2 + addnpix;ii++)
			fscanf(fp,"%ld ",&indpix[ii]);*/
		fclose(fp);
	}else{
		cerr << "Error : cannot find Indpix file " << testfile << endl;
		exit(0);
	}
	cout << "check" << endl;

	int aa=0,npix2;
	for (ii=0;ii<factdupl*nn*nn+2 + addnpix;ii++){
		if (indpix[ii] != -1){
			aa++;
		}
	}
	npix2 = aa;  // npix = number of filled pixels
	if(npix!=npix2)
		cout << "npix " << npix << "npix2 " << npix2 << endl;
	cout << "pndread" << endl;



#ifdef USE_MPI
	/********************* Define parallelization scheme   *******/

	if (rank == 0){

		check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, &frnum);
		// reorder nsamples
		//find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";
		for (ii=0;ii<ntotscan;ii++){
			nsamplesorder[ii] = nsamples[ruleorder[ii]];
			fframesorder[ii] = fframes[ruleorder[ii]];
			//extentnoiseSp_allorder[ii] = extentnoiseSp_all[ii];
		}
		for (ii=0;ii<ntotscan;ii++){
			nsamples[ii] = nsamplesorder[ii];
			fframes[ii] = fframesorder[ii];
			//extentnoiseSp_all[ii] = extentnoiseSp_allorder[ii];
			//printf("frnum[%d] = %d\n",ii,frnum[ii]);
		}

	}




	if (rank == 0){
		/*sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		for (ii=0;ii<=ntotscan;ii++) fprintf(fp,"frnum[%ld] = %ld \n",ii,frnum[ii]);
		fclose(fp);*/


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
	}
#endif


#ifdef USE_MPI

	MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);

	iframe_min = frnum[rank];
	iframe_max = frnum[rank+1];
	rank_det = 0;
	size_det = 1;

#else
	iframe_min = 0;
	iframe_max = ntotscan;
	rank_det = rank;
	size_det = size;

#endif


	PNd = new double[npix];
	PNdtot = new double[npix];
	hits=new long[npix];
	init1D_double(PNd,0,npix,0.0);
	init1D_long(hits,0,npix,0);
	init1D_double(PNdtot,0,npix,0.0);


	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//



	printf("[%2.2i] Pre-processing of the data\n",rank);



	prefixe = "fdata";

	// loop over the scans
	for (iframe=iframe_min;iframe<iframe_max;iframe++){

		ns = nsamples[iframe]; // number of samples for this scan
		ff = fframes[iframe]; //first frame of this scan
		f_lppix = f_lp*double(ns)/fsamp; // knee freq of the filter in terms of samples in order to compute fft
		f_lppix_Nk = f_lp_Nk*double(ns)/fsamp; // noise PS threshold freq, in terms of samples


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
					poutdir,termin,errarcsec,dirfile,scerr_field,flpoint_field,bolonames,
					bextension,fextension,shift_data_to_point,f_lppix,ff,ns,
					napod,ndet,NORMLIN,NOFILLGAP,iframe);// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)


			do_PtNd(PNd,bextension,poutdir,prefixe,termin, dirfile, bolonames,ff,ns,ndet,
					size_det,rank_det,indpix,nn,npix,iframe,NULL,hits);

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

			// return Pnd = At N-1 d

		} else {

			do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
					scerr_field,flpoint_field,bolonames,bextension,fextension,
					/*cextension,*/shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
					flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
					nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,NULL);

		}


	} // end of iframe loop



	cout << "naive step" << endl;
	double *map1d;
	int mi;
	map1d = new double[nn*nn];




	for (int ii=0; ii<nn; ii++) {
		for (int jj=0; jj<nn; jj++) {
			mi = jj*nn + ii;
			if (indpix[mi] >= 0){
				map1d[mi] = hits[indpix[mi]];
			} else {
				map1d[mi] = 0.0;
			}
		}
	}

	fname = '!' + outdir + "naivMap_" + termin + "_hits.fits";
	write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


	for (int ii=0; ii<nn; ii++) {
		for (int jj=0; jj<nn; jj++) {
			mi = jj*nn + ii;
			if (indpix[mi] >= 0){
				if(hits[indpix[mi]]>0)
					map1d[mi] = PNd[indpix[mi]]/(double)hits[indpix[mi]];
			} else {
				map1d[mi] = 0.0;
			}
		}
	}


	fname = '!' + outdir + "naivMap_" + termin;
	fname+= "_naive.fits";
	cout << fname << endl;
	write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


	printf("End of saneNaiv\n");

#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
	PNdtot=PNd;
#endif

	// Close MPI process


#ifdef USE_MPI
	MPI_Finalize();
#endif

	// clean up
	delete [] mask;

	free(testfile);

	delete [] fframesorder; //will be needed
	delete [] nsamplesorder; //will be needed
	delete [] ruleorder; //will be needed
	delete [] frnum; //will be needed

	delete [] PNd; //needed
	delete [] PNdtot; //needed
	delete [] fframes; // needed
	delete [] nsamples; //needed

	delete [] indpix; //needed
	delete [] indpsrc; //needed


	return 0;
}




//******************************************************************//
//******************************************************************//
//**********************  End of Sane Naiv *************************//
//******************************************************************//
//******************************************************************//



