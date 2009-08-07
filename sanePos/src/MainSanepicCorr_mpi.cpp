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


#include "mpi_architecture_builder.h"

#include "boloIO.h"
#include "dataIO.h"
#include "inline_IO2.h"
extern "C" {
#include "nrutil.h"
}


#include "sanePos_map_making.h"
#include "blastSpecific.h"
#include "parsePos.h"
#include "sanePos_preprocess.h"
#include "projection_wcs.h"


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


/*! \mainpage Sanepic for SPIRE
 *
 * \section intro_sec Matthieu HUSSON & Alexandre Beelen
 *
 * Sanepic for SPIRE Manual
 */



int main(int argc, char *argv[])
{



	int size, size_det; /*! size = number of processor used for this step*/
	int rank, rank_det; /*! rank = processor MPI rank*/

#ifdef USE_MPI
	// int tag = 10;
	MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	size = 1;
	rank = 0;
	cout << "Mpi is not used for this step" << endl;
#endif




	//default value of the data to pointing shift
	int shift_data_to_point = 0; /*! default value = 0 */


	//DEFAULT PARAMETERS
	long napod = 0; /*! number of samples to apodize, =0 -> no apodisation */
	double errarcsec = 15.0; /*! source error, rejection criteria : scerr[ii] > errarcsec, sample is rejected */


	long iframe_min, iframe_max; /*! frame number min and max each processor has to deal with */

	int flagon = 0; /*! if rejectsample [ii]==3, flagon=1*/
	bool bfixc = 0; /*! indicates that 4 corners are given for the cross corelation removal box */
	bool pixout = 0; /*! indicates that at least one pixel has been flagged and is out */
	bool NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	bool flgdupl = 0; /*! 1 if flagged data are put in a separate map */


	//set coordinate system
	double *srccoord, *coordscorner; /* srccoord = source coordinates, coordscorner = map corners coordinates*/
	srccoord = new double[2]; // RA/DEC source
	coordscorner = new double[4]; // map min/max RA/DEC coords (-N,-t,-T absents)
	srccoord[0] = -1000; // RA tangent point/source
	srccoord[1] = -1000; // DEC tangent point/source
	double radius = -1.0; /*! map radius (half a side) in degrees */


	// data parameters
	long *fframes  ; /*! first frames number array*/
	long *nsamples ; /*! number of samples table array */


	long ntotscan; /*! total number of scans */
	long ndet; /*! number of channels used*/
	int nnf; /*! number of noise file */
	long addnpix=0; /*!add a number 'n' of pixels to the map */



	// map making parameters
	double pixdeg; /*! size of pixels (degree) */

	int nn, npix; /*! nn = side of the map, npix = number of filled pixels */
	long npixsrc; /*! number of pixels included in CCR */
	double ra_min, ra_max, dec_min, dec_max; /*! ra/dec min/max coordinates of the map*/
	double *offsets, *froffsets, *offmap; /*! bolo offsets / ref bolo, froffsets = frame offsets, offmap = map offsets */
	float *scoffsets; /*! source offsets depending on wavelength */
	scoffsets = new float[6];

	int nfoff; /*! number of offsets */
	foffset *foffsets; /*! tableau d'offsets */

	double *tancoord; /*! tangent point coordinates RA/dec */
	double *tanpix; /*! tangent pixel coordinates in the map */
	double gra_min, gra_max, gdec_min, gdec_max; /*! global ra/dec min and max (to get the min and max of all ra/dec min/max computed by different processors) */

	//internal data params
	long ns; /*! number of samples for this scan */


	string fname; /*! parallel scheme file name */


	char type='d'; /*! returned type of read_data functions, d=64bit double */
	double *ra, *dec, *phi, *scerr; /*! RA/DEC, phi (angle) coordinates of the bolo, source errors */
	unsigned char *flag, *flpoint, *rejectsamp, *mask; /*! samples flags, pointing flags, rejected samples list */
	long *indpix, *indpsrc; /*! pixels indices, CCR mask pixels indices */

	int *xx, *yy; /*! data coordinates in the map */
	long *pixon; /*! this array is used to store the rules for pixels : they are seen or not */
	long *samptopix; /*! sample to pixel conversion array */




	string field; /*! actual boloname in the bolo loop */
	string bolofield; /*! bolofield = boloname + bextension */
	string flagfield; /*! flagfield = field+fextension;*/
	string dirfile; /*! data directory*/
	string outdir; /*! output directory*/
	string poutdir; /*! current path (pPath) or output dir (outdir)*/
	string bextension; /*! bolometer field extension*/
	string fextension = "NOFLAG"; /*! flag field extension*/
	string pextension; /*! pointing extension*/
	string file_offsets; /*! bolometer offsets file*/
	string file_frame_offsets = "NOOFFS"; /*! offset file*/
	//string termin; /*! output file suffix */
	string termin_internal = "internal_data";

	/* DEFAULT PARAMETERS */
	int coordsyst = 1; /*! Default is RA/DEC */


	int samples_per_frames=20;

	/* parser inputs */
	std::vector<string> bolonames/*, extentnoiseSP*/; /*! bolometer list, noise file prefix */
	std::vector<long> fframes_vec, nsamples_vec; /*! first frame list, number of frames per sample */
	std::vector<long> xxi, xxf, yyi, yyf; /*! box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y) */
	//std::vector<double> fcut;
	//std::vector<string> extentnoiseSP; /*! noise file prefix*/


	time_t t2, t3;//, t3, t4, t5, dt;




	pixdeg = -1.0; /*! "Size of pixels (deg)"*/

	// -----------------------------------------------------------------------------//

	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		//parse_sanePos_ini_file(argv[1]);
		int parsed=1;
		parsed=parse_sanePos_ini_file(argv[1],bfixc,shift_data_to_point,napod,NOFILLGAP,flgdupl,
				srccoord,coordscorner,radius,ntotscan,ndet,nnf,
				pixdeg,dirfile,outdir,poutdir,bextension,fextension,
				pextension,file_offsets,file_frame_offsets,coordsyst,bolonames,fframes_vec,nsamples_vec,fname,xxi,xxf,yyi,yyf);

		if (parsed==-1){
#ifdef USE_MPI
			MPI_Finalize();
#endif
			exit(1);
		}

	}


	///////////////: debug ///////////////////////////////
	cout << "ntotscan : " << ntotscan << endl;


	std::vector<long>::iterator it;

	cout << "frames" << endl;
	for(it=fframes_vec.begin();it<fframes_vec.end();it++)
		cout << *it << " ";

	cout << "\nnsamples" << endl;
	for(it=nsamples_vec.begin();it<nsamples_vec.end();it++)
		cout << *it << " ";
	cout << endl;
	///////////////: debug ///////////////////////////////

	// -----------------------------------------------------------------------------//
	t2=time(NULL);

	if (napod){
		printf("[%2.2i] Data are apodized\n", rank);
	} else {
		printf("[%2.2i] Data are not apodized\n", rank);
	}

	printf("[%2.2i] Data written in %s\n",rank, poutdir.c_str());

	// convert lists to regular arrays (MPI_BCas works only on array...
	fframes       = new long[ntotscan];
	nsamples      = new long[ntotscan];


	vector2array(nsamples_vec, nsamples);
	vector2array(fframes_vec,  fframes);
	//cout << fframes[0] << fframes[1] << fframes[2] << endl;
	//cout << nsamples[0] << nsamples[1] << nsamples[2] << endl;



	if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}


	// utilisé lors de la lecture des coord de la map en pixel (dans la f° read_data)
	string ra_field;
	string dec_field;
	string phi_field;
	string scerr_field = "ERR"+pextension;
	string flpoint_field = "FLPOINTING";

	if (coordsyst == 2){
		ra_field = "L"+pextension;
		dec_field = "B"+pextension;
		phi_field = "PHIG"+pextension;
		printf("[%2.2i] Coordinate system: Galactic\n",rank );
	}else{
		ra_field = "RA"+pextension;
		dec_field = "DEC"+pextension;
		phi_field = "PHI"+pextension;
		if (coordsyst == 3){
			printf("[%2.2i] Map in Telescope coordinates. Reference coordinate system is RA/DEC (J2000)\n", rank);
		} else {
			printf("[%2.2i] Coordinate system: RA/DEC (J2000)\n", rank);
		}
	}
	/*

	if (NORMLIN)
		printf("NO BASELINE REMOVED\n");*/


	/*if (projgaps)
		printf("Flaged data are binned. iterative solution to fill gaps with noise only.\n");
	 */

	/*! map offsets*/
	nfoff = map_offsets(file_frame_offsets, ntotscan, scoffsets, foffsets,fframes,rank);





#ifdef USE_MPI
	/********************* Define parallelization scheme   *******/


	if (rank == 0){

		long *frnum ;
		long *ruleorder ;
		long *fframesorder ;
		long *nsamplesorder ;
		//string *extentnoiseSp_allorder;

		check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, &frnum);
		// reorder nsamples
		//find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";


		fframesorder  = new long[ntotscan];
		//extentnoiseSp_allorder = new string[ntotscan];
		nsamplesorder = new long[ntotscan];

		for (long ii=0;ii<ntotscan;ii++){
			nsamplesorder[ii] = nsamples[ruleorder[ii]];
			fframesorder[ii] = fframes[ruleorder[ii]];
			//extentnoiseSp_allorder[ii] = extentnoiseSp_all[ruleorder[ii]];
		}
		for (long ii=0;ii<ntotscan;ii++){
			nsamples[ii] = nsamplesorder[ii];
			fframes[ii] = fframesorder[ii];
			//extentnoiseSp_all[ii] = extentnoiseSp_allorder[ii];
			//printf("frnum[%d] = %d\n",ii,frnum[ii]);
		}

		delete [] fframesorder;
		delete [] nsamplesorder;
		//delete [] extentnoiseSp_allorder;

		delete [] ruleorder;
		delete [] frnum;

	}

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




	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	/************************ Look for distribution failure *******************************/
	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}



	/* END PARAMETER PROCESSING */




	/********** Alocate memory ***********/
	printf("[%2.2i] Allocating Memory\n",rank);

	// seek maximum number of samples
	ns = nsamples[0];
	for(int ii=0;ii<ntotscan;ii++) if (nsamples[ii] > ns) ns = nsamples[ii];

	ra = new double[2*ns]; // RA bolo de ref
	dec = new double[2*ns]; // DEc du bolo de ref
	phi = new double[2*ns]; // (du bolo de ref) angle de la matrice de detecteur par rapport a RA/dec
	scerr = new double[2*ns]; // BLAST SPECIFIC : mesure l'erreur de pointage, si trop grande on flag la donnée
	xx = new int[2*ns]; // sample column coordinates in the map
	yy = new int[2*ns]; // sample row coordinates in the map
	samptopix = new long[2*ns]; // sample to pixel conversion index
	flag = new unsigned char[2*ns]; // flag data => =1
	rejectsamp = new unsigned char[2*ns]; // rejected samples after flag conditions
	flpoint = new unsigned char[2*ns]; // flpoint est un flag du pointage/time. Savoir au temps t, si tu prends ces données là, ou non.
	tancoord = new double[2]; // coordinates in ra/dec of the tangent point
	tanpix = new double[2]; // coordinates in the map of the tangent point

	froffsets = new double[2]; //
	offsets = new double[2];

	offmap = new double[2]; // map offsets


	// default value for map variables
	ra_min  = 1000.0;
	ra_max  = -1000.0;
	dec_min = 1000.0;
	dec_max = -1000.0;

	offmap[0] = 0.0;
	offmap[1] = 0.0;

	/*--------------------------------Default parameters initialization for WCS and fits headers ---------*/
	/*const int NAXIS = 2; // number of dimensions
	const double CRPIX[2] =  {  128.0, 128.0}; // reference pixel coordinates
	const double PC[2][2] = {{    1.0,  0.0}, // linear transformation matrix
			{    0.0,  1.0}};
	const double CDELT[2] =  {-0.00168725828819, 0.00168725828819}; // -pixdeg, pixdeg (coordinate scales)

	char CUNIT[2][18] = {"deg", "deg"}; // en degre // units of CRVAL and CDELT
	char CTYPE[2][18] = {"RA---TAN", "DEC--TAN"}; // X, Y, projection type
	char CNAME[2][18] = {"Right Ascension", "Declination"};

	const double CRVAL[2] = {258.144877,-38.332616}; // ra mean, dec mean
	// => Celestial longitude and latitude of the fiducial point

	const double LONPOLE  = 150.0; // reference pole longitude
	const double LATPOLE  = 150.0; // reference pole latitude

	int NPV = 0;*/
	//struct pvcard PV[2]; // native latitude, longitude of the fiducial point /* Projection parameters are set in main(). */

	/*-------------------------------------------------------------------------------------------------*/




	//********************************************************************************
	//*************  find coordinates of pixels in the map
	//********************************************************************************

	printf("[%2.2i] Finding coordinates of pixels in the map\n",rank);


	//if (coordsyst != 4){ // coordsyst never = 4 => debug mode => delete
	/*!
	 * \fn find_coordinates_in_map : Output : ra_min, ra_max, dec_min, dec_max
	 * -> Compute map coordinates
	 */
	find_coordinates_in_map(ndet,bolonames,bextension,fextension,file_offsets,foffsets,scoffsets,
			/*offsets,*/iframe_min,iframe_max,fframes,nsamples,dirfile,ra_field,dec_field,phi_field,
			scerr_field,flpoint_field,nfoff,pixdeg, xx, yy, nn, coordscorner,
			tancoord, tanpix, bfixc, radius, offmap, srccoord, type,ra,dec,phi,scerr, flpoint,ra_min,ra_max,dec_min,dec_max);


#ifdef USE_MPI
	MPI_Reduce(&ra_min,&gra_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&ra_max,&gra_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(&dec_min,&gdec_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&dec_max,&gdec_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&gra_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gra_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gdec_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gdec_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#else
	gra_min=ra_min;
	gra_max=ra_max;
	gdec_min=dec_min;
	gdec_max=dec_max;
#endif

	//set coordinates
	coordscorner[0] = gra_min; // store ra/dec min/max of the final map
	coordscorner[1] = gra_max;
	coordscorner[2] = gdec_min;
	coordscorner[3] = gdec_max;

	if (rank == 0) {
		printf("[%2.2i] ra  = [ %7.3f, %7.3f ] \n",rank, gra_min, gra_max );
		printf("[%2.2i] dec = [ %7.3f, %7.3f ] \n",rank, gdec_min, gdec_max);
	}

	/// just to set nn in order to compute map-making matrices and vectors
	sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
			tancoord, tanpix, 1, radius, offmap, srccoord,0);


	/*!
	 * \fn write Pointing informations in a file
	 * Write nn : size of the map in pixel
	 * outdir : output directory
	 * termin : generated files prefixe
	 * coordsyst : coordinate system value
	 * tanpix : tangent point coordinates in the map (in pixel)
	 * tancoord : tangent point coordinates in coordsyst coordinate system
	 */
	write_info_pointing(nn, outdir, termin_internal, coordsyst, tanpix, tancoord);


	/*} else {
		// read those parameters from a file : -c = 4 option
		sprintf(testfile,"%s%s%s%s%d%s",outdir.c_str(),"InfoPointing_for_Sanepic_",termin.c_str(),"_", rank,".txt");
		if ((fp = fopen(testfile,"r")) == NULL){
			cerr << "File InfoPointing_for_sanepic... not found. Exiting" << endl;
			exit(1);
		}
		fscanf(fp,"%d\n",&nn);
		fscanf(fp,"%d\n",&coordsyst);
		fscanf(fp,"%lf\n",tanpix);
		fscanf(fp,"%lf\n",tanpix+1);
		fscanf(fp,"%lf\n",tancoord);
		fscanf(fp,"%lf\n",tancoord+1);
		fclose(fp);
	}*/







	//************************************* Deal with masking the point sources
	mask = new unsigned char[nn*nn];
	indpsrc = new long[nn*nn];

	// if a box for crossing constraint removal is given in ini file
	if(xxi.size()>0)
		/*!
		 * \fn compute_indpsrc_addnpix
		 * Input : nn : size of the map in pixels
		 * ntotscan : total number of scans
		 * xxi, xxf, yyi, yyf : box for crossing constraint removal corner coordinates
		 * indpsrc : pixel indice of the pixel included in box for CCR
		 * npixsrc : total number of pix in box for CCR
		 * mask : 0 if the pixel is contained in box for CCR, 1 otherwise
		 */
		addnpix=Compute_indpsrc_addnpix(nn,ntotscan,xxi,xxf,yyi,yyf,indpsrc,npixsrc,mask);
	else{
		for (long ii=0;ii<nn*nn;ii++)
			mask[ii] = 1;
		addnpix=0;
		npixsrc=0;
		for(long ii=0;ii<nn*nn;ii++)
			indpsrc[ii]=-1;
	}


	// map duplication factor
	int factdupl = 1;
	if (flgdupl) factdupl = 2; //  default 0 : if flagged data are put in a duplicated map


	//pixon indicates pixels that are seen
	pixon = new long[factdupl*nn*nn+2 + addnpix];   // last pixel is for flagged samples
	fill(pixon,pixon+(factdupl*nn*nn+2 + addnpix),0);





	//**********************************************************************************
	// get coordinates of pixels that are seen
	//**********************************************************************************
	/*!
	 * \fn Get coordinates of pixels that are seen
	 * Compute the position to pixel projetcion matrices :
	 * One binary file per bolometer and per scan
	 */
	compute_seen_pixels_coordinates(ndet,ntotscan,outdir,bolonames,bextension, fextension, termin_internal,
			file_offsets,foffsets,scoffsets,iframe_min, iframe_max,fframes,
			nsamples,dirfile,ra_field,dec_field,phi_field, scerr_field,
			flpoint_field, nfoff,pixdeg,xx,yy,mask, nn,coordscorner, tancoord,
			tanpix, bfixc, radius, offmap, srccoord, type, ra,dec,
			phi,scerr,flpoint,shift_data_to_point,ra_min,ra_max,dec_min,dec_max, flag,
			napod, errarcsec, NOFILLGAP, flgdupl,factdupl, addnpix, rejectsamp, samptopix, pixon, rank, indpsrc, npixsrc, flagon, pixout);



	bool default_projection = 1;
	string temp = dirfile + "optimMap_sanepic_flux.fits";
	const char *fits_file = temp.c_str();
	fits_header_generation(outdir,fits_file,pixdeg,default_projection,tanpix,tancoord);


	//************** init mapmaking variables *************//

	printf("[%2.2i] Init map making variables\n",rank);


	// pixel indices
	indpix = new long[factdupl*nn*nn+2 + addnpix];

	if(rank==0){

		/*!
		 *  \fn compute indpix : pixels indices used for data projection/deprojection.
		 *  npix is the number of filled pixels
		 *  Pixon : indicates if the pixel is projected on the map or not
		 *  addnpix : number of pixels to add : depends on box for CCR and duplication factor
		 */
		npix=compute_indpix(indpix,factdupl,nn,addnpix,pixon);


		// write in a file for conjugate gradient step
		long ind_size = factdupl*nn*nn+2 + addnpix;

		/*!
		 * Write indpix to a binary file : ind_size = factdupl*nn*nn+2 + addnpix;
		 * npix : total number of filled pixels,
		 * flagon : if some pixels are apodized or outside the map
		 */
		write_indpix(ind_size, npix, indpix, termin_internal, outdir, flagon);
	}
	/*
	testfile2 = outdir + "Indpix_for_conj_grad_" + termin + ".txt";
	fp = fopen(testfile2.c_str(),"w");
	//	long indpix_size = factdupl*nn*nn+2 + addnpix;
	//	fwrite(&indpix_size,sizeof(long),1,fp);
	//fwrite(&flagon,sizeof(int),1,fp); // mat 04/06
	fprintf(fp,"%d ",flagon);
	fprintf(fp,"%d ",npix);
	for(int gg=0;gg<nn*nn+2;gg++)
		fprintf(fp,"%ld ",indpix[gg]);
	fclose(fp);*/

	//  printf("[%2.2i] indpix[nn*nn] = %d\n",rank, indpix[nn*nn]);


	if (pixout)
		printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
	printf("[%2.2i] Total number of detectors : %d\t Total number of Scans : %d \n",rank,(int)ndet, (int) ntotscan);
	printf("[%2.2i] Size of the map : %d x %d\t Total Number of filled pixels : %d\n",rank, nn, nn, npix);


	t3=time(NULL);
	cout << "temps de traitement : " << t3-t2 << " sec" << endl;




#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	/* ---------------------------------------------------------------------------------------------*/


	// Close MPI process


#ifdef USE_MPI
	MPI_Finalize();
#endif

	printf("[%2.2i] Cleaning up\n",rank);


	// clean up
	delete [] ra;
	delete [] dec;
	delete [] phi;
	delete [] scerr; //scerr_field needed
	delete [] xx;
	delete [] yy;
	delete [] flag;
	delete [] rejectsamp;
	delete [] samptopix;
	delete [] flpoint; // flpoint_field needed but not flpoint
	delete [] mask;

	delete [] pixon;

	delete [] scoffsets;
	delete [] offsets;
	delete [] froffsets;
	delete [] offmap;

	delete [] foffsets;
	delete [] srccoord;
	delete [] coordscorner;

	//free(testfile);


	//delete [] frnum;

	delete [] fframes;
	delete [] nsamples;
	delete [] tancoord;
	delete [] tanpix;
	delete [] indpix;
	delete [] indpsrc;

	printf("[%2.2i] End of sanePos\n",rank);

	return 0;
}
