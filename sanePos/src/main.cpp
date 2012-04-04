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
#include <sysexits.h>

#include "mpi_architecture_builder.h"
#include "struct_definition.h"

#include "dataIO.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "error_code.h"

extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
#include "getdata.h"
}


#include "sanePos_map_making.h"
#include "sanePos_preprocess.h"
#include "parser_functions.h"
#include "coord.h"

#ifdef PARA_FRAME
#include "mpi.h"
#endif

using namespace std;

/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *	- Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *  - For each file :
 *      - Generate or clear the dirfile parts that will be filled : indices
 *
 *  - Get input fits META DATA
 *
 *  - Compute map extrema or read it from binary mask file
 *  - Compute map Header or read it from binary mask file
 *  - Compute masked pixels indice table and save it to disk
 *  - Save mapHeader to disk for the other programs
 *  - Compute pixels indice table and save it to disk
 *
 *  - Compute Naïve map including Image map, Coverage map, history table and METADATA header
 *
 */

int main(int argc, char *argv[])
/* Main sanePos function */
{

	int size; /* size = number of processor used for this step*/
	int rank; /* rank = processor MPI rank*/

#ifdef PARA_FRAME
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size); // get mpi number of processors
	MPI_Comm_rank(MPI_COMM_WORLD,&rank); // each processor is given a number
#else
	size = 1; // non-MPI usage : 1 processor with number 0
	rank = 0;
#endif

	if(rank==0)
	  cout << endl << "sanePos: computation of pixel indices" << endl;

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_saneProc proc_param;
	struct param_sanePos pos_param;
	struct param_common dir; /* structure that contains output input temp directories */

	int flagon = 0; /* if rejectsample [ii]==3, flagon=1*/
	bool pixout = 0; /* indicates that at least one pixel has been flagged and is out */

	long long npix; /* npix = number of filled pixels */
	long long npixsrc; /* number of pixels included in Crossing Constraint Removal */
	long long addnpix=0; /* add a number 'n' of pixels to the map */

	int nwcs=1;                      // We will only deal with one wcs....
	struct wcsprm * wcs;             // wcs structure of the image
	long NAXIS = 2, NAXIS1, NAXIS2;  // size of the image

	char * subheader;               // Additionnal header keywords
	int nsubkeys;                   //



	// System should be IEEE 754 complient (TODO : add in the doc)
	double lon_min=NAN, lon_max=NAN, lat_min=NAN, lat_max=NAN; /* min/max coordinates of the map*/
	double glon_min, glon_max, glat_min, glat_max; /* global  min and max (to get the min and max of all min/max computed by different processors) */

	string fname; /* parallel scheme file name */

	// positions variables
	short *mask;
	long long *indpix, *indpsrc; /* pixels indices, CCR mask pixels indices */
	long long *pixon; /* this array is used to store the rules for pixels : they are seen or not */
	long long *pixon_tot=NULL;

	//naiv map data params
	long ns; /* number of samples for this scan, first frame number of this scan*/
	double fhp_pix, fcut_pix; /* frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples*/
	double *PNdNaiv, *PNdtotNaiv=NULL; /*  projected noised data, and global Pnd for mpi utilization */
	long *hitsNaiv, *hitstotNaiv=NULL; /* naivmap parameters : hits count */

	// struct used in the parser
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;
	string parser_output="";

	std::vector<string> key;
	std::vector<int> datatype;
	std::vector<string> val;
	std::vector<string> com;

	std::vector<std::vector<std::string> > bolo_list; // this vector contains all bolonames for all the scans

	uint16_t mask_sanePos = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | PIXDEG_WRONG_VALUE | FILEFORMAT_NOT_FOUND | NAPOD_WRONG_VALUE |
			FHP_PROBLEM | FITS_FILELIST_NOT_FOUND | FCUT_PROBLEM; // 0xc2ff

	string field; /* field = actual boloname in the bolo loop */

	// -----------------------------------------------------------------------------//

	//	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	// Parse ini file
	if (argc<2) {
		compare_to_mask=0x001;
	} else {

		parsed=parser_function(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, struct_sanePic, size, rank);

		compare_to_mask = parsed & mask_sanePos;

		// print parser warning and/or errors
		if (rank == 0)
			cout << endl << parser_output << endl;
	}

	if(compare_to_mask>0x0000){

		switch (compare_to_mask){/* error during parsing phase */

		case 0x0001: printf("Please run %s using a correct *.ini file\n",argv[0]);
		break;

		default : printf("Wrong program options or argument. Exiting !\n");
		break;


		}

#ifdef PARA_FRAME
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return EX_CONFIG;
	}
	//	}


	// -----------------------------------------------------------------------------//
#ifdef DEBUG
	time_t t2=time(NULL);
#endif


#ifdef PARA_FRAME

	MPI_Barrier(MPI_COMM_WORLD);

	if(configure_PARA_FRAME_samples_struct(dir.tmp_dir, samples_struct, rank, size)){
		MPI_Abort(MPI_COMM_WORLD, 1);
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif

	if(channel_list_to_vect_list(samples_struct, bolo_list, rank)){
		cout << "error in channel_list_to_vect_list" << endl;
		return EX_CONFIG;
	}

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic, saneInv_struct);

		cleanup_dirfile_sanePos(dir.tmp_dir, samples_struct, bolo_list);
	}

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
#endif


	// Open the dirfile to read temporary files
	string filedir = dir.tmp_dir + "dirfile";
	samples_struct.dirfile_pointer = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED);

	if (gd_error(samples_struct.dirfile_pointer) != 0) {
		cout << "error opening dirfile : " << filedir << endl;
#ifdef PARA_FRAME
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return 1;
	}

	//	read pointing header
	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, &subheader, &nsubkeys, rank)){ // read keyrec file
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return (EX_IOERR);
	}

	//	if (rank == 1){
	//		print_MapHeader(wcsFake);
	//		cout << "nkeys" << nsubkeys << endl;
	//		char * hptr;
	//		// .... print it
	//		hptr = subheader;
	//		printf("\n\n Sub Header :\n");
	//		for (int ii = 0; ii < nsubkeys; ii++, hptr += 80) {
	//			printf("%.80s\n", hptr);
	//		}
	//	}

	// TODO: Convert Position if needed here....
	if (pos_param.eq2gal || pos_param.gal2eq) {
		if (rank == 0)
			cout << endl<< "Converting coordinates..." << endl;

		int status = convert_Dirfile_LON_LAT(samples_struct, pos_param, bolo_list);

		if (pos_param.eq2gal)
			pos_param.axistype = "GAL";
		if (pos_param.gal2eq)
			pos_param.axistype = "EQ";

		if (status) {
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return 1;
		}
	}


	if (pos_param.maskfile == ""){

		//TODO Find an easy way to provide for a map header (without going through mask...)

		if(rank==0)
			printf("\nDetermining Map Parameters...\n");

		// Populate the wcs structure ...
		// ... add pixel size in deg ...
		for (int ii = 0; ii < NAXIS; ii++) wcs->cdelt[ii] = (ii) ? pos_param.pixdeg: -1*pos_param.pixdeg ;
		for (int ii = 0; ii < NAXIS; ii++) strcpy(wcs->cunit[ii], "deg");

		// ... add axis label ...
		if (pos_param.axistype.compare("EQ") == 0){
			char TYPE[2][5]  = { "RA--", "DEC-"};
			char NAME[2][16] = {"Right Ascension","Declination"};

			for (int ii = 0; ii < NAXIS; ii++) {
				strcpy(wcs->ctype[ii], &TYPE[ii][0]);
				strncat(wcs->ctype[ii],"-",1);
				strncat(wcs->ctype[ii],pos_param.projcode.c_str(), 3);
				strcpy(wcs->cname[ii], &NAME[ii][0]);
			}
		}

		if (pos_param.axistype.compare("GAL") == 0){
			char TYPE[2][5]  = { "GLON", "GLAT"};
			char NAME[2][19] = {"Galactic Longitude", "Galactic Latitude"};

			for (int ii = 0; ii < NAXIS; ii++) {
				strcpy(wcs->ctype[ii], &TYPE[ii][0]);
				strncat(wcs->ctype[ii],"-",1);
				strncat(wcs->ctype[ii],pos_param.projcode.c_str(), 3);
				strcpy(wcs->cname[ii], &NAME[ii][0]);
			}
		}


		// Define a projection center ...
		// ... From the ini file ...
		if ( ! isnan(pos_param.lon) && ! isnan(pos_param.lat) ) {
			wcs->crval[0] = pos_param.lon;
			wcs->crval[1] = pos_param.lat;
		} else {
			// ... or the first valid data point
			double lon_center, lat_center;

			if(rank==0) {

				double *lon, *lat;
				int *flag;
				long ns = samples_struct.nsamples[0];

				lon  = new double[ns];
				lat  = new double[ns];
				flag = new int[ns];

				if(read_LON_from_dirfile(samples_struct.dirfile_pointer, samples_struct.basevect[0], bolo_list[0][0], lon, ns))
					return 1;
				if(read_LAT_from_dirfile(samples_struct.dirfile_pointer, samples_struct.basevect[0], bolo_list[0][0], lat, ns))
					return 1;
				if(read_flag_from_dirfile(samples_struct.dirfile_pointer, samples_struct.basevect[0], bolo_list[0][0], flag, ns))
					return 1;

				int ii = 0;
				while (flag[ii] != 0){ ii++; }
				lon_center = lon[ii];
				lat_center = lat[ii];

				delete [] lon;
				delete [] lat;
				delete [] flag;

			}

#ifdef PARA_FRAME
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&lon_center,   1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&lat_center,  1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			wcs->crval[0] = lon_center;
			wcs->crval[1] = lat_center;

			// In these case, the projection center has not been defined by the user, so we need to find a proper one
			// so...
			if (rank == 0)
				cout << endl << "WW - Projections will have to be made twice, define a projection center..." << endl;
			// set the structure to have the celestial projection routines
			if ( int status = wcsset(wcs) ) {
				printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
			}

			// Project all data once using the temporary projection center...
			// at this stage all shall be in dirfile, so HIPE format...
			if(samples_struct.iframe_min!=samples_struct.iframe_max){
//				switch (pos_param.fileFormat) {
//				case 0:
//					if(computeMapMinima(samples_struct, dir.data_dir,
//							iframe_min,iframe_max, wcs,
//							lon_min,lon_max,lat_min,lat_max, bolo_list)){
//#ifdef PARA_FRAME
//						MPI_Abort(MPI_COMM_WORLD, 1);
//#endif
//						return(EX_OSERR);
//					}
//					break;
//				case 1:
					if(computeMapMinima_HIPE(samples_struct,
							wcs, lon_min,lon_max,lat_min,lat_max, bolo_list)){
#ifdef PARA_FRAME
						MPI_Abort(MPI_COMM_WORLD, 1);
#endif
						return(EX_OSERR);
					}
//					break;
//				}
			}

#ifdef PARA_FRAME

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&lon_min,&glon_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
			MPI_Reduce(&lon_max,&glon_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Reduce(&lat_min,&glat_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
			MPI_Reduce(&lat_max,&glat_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&glon_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&glon_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&glat_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&glat_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#else
			glon_min=lon_min;
			glon_max=lon_max;
			glat_min=lat_min;
			glat_max=lat_max;
#endif


			double x_mean, y_mean, phi, theta, lon_mean, lat_mean;
			int status;

			x_mean = (glon_max+glon_min)/2;
			y_mean = (glat_max+glat_min)/2;

			// We can now define a proper projection center
			if (celx2s(&(wcs->cel), 1, 0, 0, 0, &x_mean, &y_mean, &phi, &theta, &lon_mean, &lat_mean, &status) == 1) {
				printf("ERROR 1: %s\n", prj_errmsg[1]);
			}

			wcs->crval[0] = lon_mean;
			wcs->crval[1] = lat_mean;

			if (rank == 0) {
				cout << "WW - Nominal Projection Center : " << endl;
				cout << "     lon : " << lon_mean << endl;
				cout << "     lat : " << lat_mean << endl << endl;
			}
		}


		// We now have a good projection center, so
		//  ... set the structure to have the celestial projection routines
		if ( int status = wcsset(wcs) ) {
			printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
		}

		// At this stage all should be in DIRFILE so ... HIPE format
		if(samples_struct.iframe_min!=samples_struct.iframe_max){
//			switch (pos_param.fileFormat) {
//			case 0:
//				if(computeMapMinima(samples_struct, dir.data_dir,
//						iframe_min,iframe_max, wcs,
//						lon_min,lon_max,lat_min,lat_max, bolo_list)){
//#ifdef PARA_FRAME
//					MPI_Abort(MPI_COMM_WORLD, 1);
//#endif
//					return(EX_OSERR);
//				}
//				break;
//			case 1:
				if(computeMapMinima_HIPE(samples_struct,
						wcs, lon_min,lon_max,lat_min,lat_max, bolo_list)){
#ifdef PARA_FRAME
					MPI_Abort(MPI_COMM_WORLD, 1);
#endif
					return(EX_OSERR);
				}
//				break;
//			}
		}

#ifdef PARA_FRAME

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&lon_min,&glon_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Reduce(&lon_max,&glon_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Reduce(&lat_min,&glat_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Reduce(&lat_max,&glat_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&glon_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&glon_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&glat_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&glat_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#else
		glon_min=lon_min;
		glon_max=lon_max;
		glat_min=lat_min;
		glat_max=lat_max;
#endif


#ifdef DEBUG
		if (rank == 0) {
			printf("lon  = [ %7.3f, %7.3f ] \n", glon_min/15, glon_max/15);
			printf("lat = [ %7.3f, %7.3f ] \n", glat_min, glat_max);
		}
#endif

		int margingPixel = 1;  // for cosmetic...
		NAXIS1 = ceil(glon_max/pos_param.pixdeg)-floor(glon_min/pos_param.pixdeg)+2*margingPixel;
		NAXIS2 = ceil(glat_max/pos_param.pixdeg)-floor(glat_min/pos_param.pixdeg)+2*margingPixel;

		wcs->crpix[0] =    glon_max/pos_param.pixdeg + margingPixel + 1;
		wcs->crpix[1] = -1*glat_min/pos_param.pixdeg + margingPixel + 1;

		if (int wcsstatus = wcsset(wcs)) {
			printf("wcsset ERROR %d: %s.\n", wcsstatus, wcs_errmsg[wcsstatus]);
		}

		npixsrc = 0;
		addnpix = 0;

		// Initialize empty masks
		mask    = new short[NAXIS1*NAXIS2];
		indpsrc = new long long[NAXIS1*NAXIS2];

		for (long long ii=0; ii<NAXIS1*NAXIS2; ii++){
			mask[ii]    =  0;
			indpsrc[ii] = -1;
		}

	} else {
		// Map header is determined from the mask file
		if(rank==0)
			cout << "Reading Mask map : " << pos_param.maskfile << endl;

		if (read_mask_wcs(dir.input_dir + pos_param.maskfile, "mask", wcs, NAXIS1, NAXIS2, mask )){
			cerr << "Error Reading Mask file" << endl;
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_IOERR);
		}

		npixsrc = 0;
		indpsrc = new long long[NAXIS1*NAXIS2];
		long long ll;

		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				ll = NAXIS1*jj+ii;
				if (mask[ll] != 0)
					indpsrc[ll] = npixsrc++;
				else
					indpsrc[ll] = -1;
			}
		}
//		if(rank==0)
//			cout << "Computing mask intersection... " << endl;
//		count_pixsrc = new long long[npixsrc*samples_struct.ntotscan];
//		fill(count_pixsrc,count_pixsrc+(npixsrc*samples_struct.ntotscan),0);
//
//		if(computeMaskIntersection(samples_struct, pos_param, iframe_min, iframe_max, rank,
//				wcs, NAXIS1, NAXIS2, mask, indpsrc, npixsrc, bolo_list, count_pixsrc)){
//#ifdef PARA_FRAME
//			MPI_Abort(MPI_COMM_WORLD, 1);
//#endif
//			return(EX_OSERR);
//		}
//
//#ifdef PARA_FRAME
//	if(rank==0){
//		count_pixsrc_tot = new long long[npixsrc*samples_struct.ntotscan];
//	}
//	MPI_Reduce(count_pixsrc,count_pixsrc_tot,npixsrc*samples_struct.ntotscan,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
//#else
//	count_pixsrc_tot=count_pixsrc;
//#endif

	}

	// TODO : test and run with real data !!!
	//	if(modify_mask_flag_in_dirfile(dir.tmp_dir, samples_struct, bolo_list, indpsrc,
	//			NAXIS1, NAXIS2, iframe_min, iframe_max)){
	//		cout << "ERROR in  modify_mask_flag_in_dirfile... Exiting...\n";
	//#ifdef PARA_FRAME
	//		MPI_Abort(MPI_COMM_WORLD, 1);
	//#endif
	//		return(EX_CANTCREAT);
	//	}

	if (rank == 0) {
		printf("Map Size         : %ld x %ld pixels\n", NAXIS1, NAXIS2);
		if(save_keyrec(dir.tmp_dir,wcs, NAXIS1, NAXIS2, subheader, nsubkeys)){
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}
	}


	// each frame contains npixsrc pixels with index indsprc[] for which
	// crossing constraint are removed
	// thus
	// addnpix = number of pix to add in pixon
	//         = number of scans * number of pix in box crossing constraint removal
	addnpix = samples_struct.ntotscan*npixsrc;

	// map duplication factor
	int factdupl;
	(pos_param.flgdupl) ? factdupl = 2: factdupl = 1; //  default 1 : if flagged data are put in a duplicated map

	// pixon indicates pixels that are seen
	// factdupl if flagged data are to be projected onto a separate map
	// 1  pixel for flagged data
	// 1  pixel for all data outside the map
	// +1 because this all are indices and C is indexing between 0 and sky_size-1
	long long sky_size = factdupl*NAXIS1*NAXIS2 + 1 + 1 + addnpix + 1;

	pixon = new long long[sky_size];
	fill(pixon,pixon+(sky_size),0);

	//**********************************************************************************
	// Compute pixels indices
	//**********************************************************************************



	if(rank==0)
		cout << endl << "Computing pixel indices..." << endl;

	// At this stage all should be in DIRFILE format, so HIPE...
//	switch (pos_param.fileFormat) {
//	case 0:
//		if(computePixelIndex(dir.tmp_dir, dir.data_dir, samples_struct,
//				proc_param, pos_param, iframe_min, iframe_max,
//				wcs, NAXIS1, NAXIS2,
//				mask,factdupl,
//				addnpix, pixon, rank,
//				indpsrc, npixsrc, flagon, pixout, bolo_list)){
//#ifdef PARA_FRAME
//			MPI_Abort(MPI_COMM_WORLD, 1);
//#endif
//			return(EX_OSERR);
//		}
//		break;
//	case 1:
		if(computePixelIndex_HIPE(samples_struct,
				proc_param, pos_param,
				wcs, NAXIS1, NAXIS2,
				mask,factdupl,
				addnpix, pixon, rank,
				indpsrc, npixsrc, flagon, pixout, bolo_list)){
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_OSERR);
		}
//		break;
//	}



#ifdef PARA_FRAME
	if(rank==0){
		pixon_tot = new long long[sky_size];
		fill(pixon_tot,pixon_tot+sky_size,0);
	}
	MPI_Reduce(pixon,pixon_tot,sky_size,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#else
	pixon_tot=pixon;
#endif

	indpix = new long long[sky_size];
	fill(indpix,indpix+sky_size,0);

	if(rank==0){

		npix = 0;
		for(long long ii=0; ii< sky_size; ii++)
			if (pixon_tot[ii] != 0)
				indpix[ii] = npix++;
			else
				indpix[ii] = -1;


		/* Write indpix to a binary file : ind_size = factdupl*nn*nn+2 + addnpix;
		 * npix : total number of filled pixels,
		 * flagon : if some pixels are apodized or outside the map
		 */
		if(write_indpix(sky_size, npix, indpix, dir.tmp_dir, flagon)){
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}
		if(write_indpsrc((long long) NAXIS1*NAXIS2, npixsrc, indpsrc,  dir.tmp_dir)){
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}

		//		printf("Total number of Scans : %d \n", (int) samples_struct.ntotscan);
		printf("Pixels Indices   : %9lld (%6.1f Mo)\n", sky_size, sky_size*8./1024/1024);
		printf("Filled Pixels    : %9lld (%6.1f Mo)\n", npix, npix*8./1024/1024);

		if (pixout)
			printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
	}


#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast( &npix,       1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(indpix,sky_size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

	//--------------------------------------- NAIVE MAP COMPUTATION ------------------------------------------//

	// (At N-1 D) memory allocation
	PNdNaiv= new double[npix];
	hitsNaiv=new long[npix];

	fill(PNdNaiv,PNdNaiv+npix,0.0);
	fill(hitsNaiv,hitsNaiv+npix,0);


#ifdef PARA_FRAME

	if(rank==0){	// global (At N-1 D) malloc for mpi
		PNdtotNaiv = new double[npix];
		hitstotNaiv=new long[npix];
	}

#endif


	if(rank==0)
		printf("\nComputing Naïve map...\n");

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	// loop over the scans
	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){

		ns       = samples_struct.nsamples[iframe]; // number of samples for this scan
		fhp_pix  = samples_struct.fhp[iframe]  * double(ns)/samples_struct.fsamp[iframe]; // knee freq of the filter in terms of samples in order to compute fft
		fcut_pix = samples_struct.fcut[iframe] * double(ns)/samples_struct.fsamp[iframe]; // noise PS threshold freq, in terms of samples

		std::vector<string> det_vect = bolo_list[iframe];
		long ndet = (long)det_vect.size();

		int pb=0;
		pb+=do_PtNd_Naiv(samples_struct, PNdNaiv, dir.tmp_dir, det_vect, ndet, proc_param.poly_order, proc_param.napod, fhp_pix, ns, indpix, iframe, hitsNaiv);
		if(pb>0){
			cout << "Problem after do_PtNd_Naiv. Exiting...\n";
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return -1;
		}

	} // end of iframe loop


#ifdef PARA_FRAME
	MPI_Reduce(PNdNaiv,PNdtotNaiv,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(hitsNaiv,hitstotNaiv,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);

#else
	PNdtotNaiv=PNdNaiv;
	hitstotNaiv=hitsNaiv;

#endif

	if (rank == 0){

		string fnaivname = dir.output_dir + "naivMap.fits";

		cout << "Output file      : naivMap.fits" << endl;

		long long mi;
		double *map1d_d;
		long *map1d_l;

		map1d_d = new double[NAXIS1*NAXIS2];
		map1d_l = new long[NAXIS1*NAXIS2];

		// TODO: Save the map of flag data if needed

		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d_d[mi] = PNdtotNaiv[indpix[mi]]/(double)hitstotNaiv[indpix[mi]];
				} else {
					map1d_d[mi] = NAN;
				}
			}
		}

		if(write_fits_wcs("!" + fnaivname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d_d,"Image",0, subheader, nsubkeys)){ // open naive Map fits file and fill ultra naive map image
			cerr << "Error Writing Ultra Naiv map ... \n";
		}


		for (long jj=0; jj<NAXIS2; jj++) {
			for (long ii=0; ii<NAXIS1; ii++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d_l[mi] = hitstotNaiv[indpix[mi]];
				} else {
					map1d_l[mi] = 0;
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
							map1d_l[mi] += hitstotNaiv[indpix[ll]];
					}
				}
			}
		}

		if(	write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, 'l', (void *)map1d_l,"Coverage",1, subheader, nsubkeys)){ // open naive Map fits file and fill hit (or coverage) image
			cerr << "Error Writing coverage map  ... \n";
		}

		if( write_fits_inifile(fnaivname, dir, proc_param, pos_param, structPS, struct_sanePic, saneInv_struct)) // write saneProc parameters in naive Map fits file header
			cerr << "WARNING ! No ini file will be included in the file : " << fnaivname << endl;
		if( write_fits_inputfile(fnaivname, samples_struct))
			cerr << "WARNING ! No input files will be included in the file : " << fnaivname << endl;


		if (pos_param.maskfile != "")
			if(copy_fits_mask(fnaivname, dir.input_dir + pos_param.maskfile)) // copy mask in naive map file
				cerr << "Warning ! The mask will not be included in naive map fits file ...\n";

		delete [] map1d_d;
		delete [] map1d_l;

	}
	/* ---------------------------------------------------------------------------------------------*/

#ifdef DEBUG
	if(rank==0){
		//Get processing time
		time_t t3=time(NULL);
		printf("\nProcessing time : %d sec\n",(int)(t3-t2));
	}
#endif


	// clean up
	delete [] PNdNaiv;
	delete [] hitsNaiv;
	delete [] indpix;
	delete [] indpsrc;
	delete [] mask;
	delete [] pixon;


	wcsvfree(&nwcs, &wcs);

	if (gd_close(samples_struct.dirfile_pointer))
		cout << "error closing dirfile : " << filedir << endl;


#ifdef DEBUG
	time_t t3=time(NULL);
	cout << "Total Time : " << t3-t2 << " sec\n";
#endif


#ifdef PARA_FRAME
	if(rank==0 && size > 1){
		// clean up
		delete [] PNdtotNaiv;
		delete [] hitstotNaiv;
		delete [] pixon_tot;
	}

	MPI_Finalize();
#endif

	if(rank==0)
		printf("\nEnd of sanePos\n");

	return EXIT_SUCCESS;
}
