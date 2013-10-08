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

#include "MPIConfiguration.h"
#include "StructDefinition.h"
#include "DataIO.h"
#include "ImageIO.h"
#include "TemporaryIO.h"
#include "InputFileIO.h"
#include "ErrorCode.h"
#include "ParserFunctions.h"
#include "Utilities.h"

extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
#include "getdata.h"
}


#include "SanePosMapMaking.h"
#include "SanePosPreProcess.h"
#include "Coord.h"

#ifdef USE_MPI
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
 *  - Compute Naive map including Image map, Coverage map, history table and METADATA header
 *
 */

int main(int argc, char *argv[]) {


	int      rank,       size; /* MPI processor rank and MPI total number of used processors */
	int  bolo_rank, bolo_size; /* As for parallel scheme */
	int  node_rank, node_size; /* On a node basis, same as *sub* but for frame scheme */

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Comm MPI_COMM_NODE, MPI_COMM_MASTER_NODE, MPI_COMM_MASTER_BOLO;
#else
	size = 1;
	rank = 0;
	bolo_size  = 1;
	bolo_rank  = 0;
	node_size = 1;
	node_rank = 0;
#endif

	if(rank==0)
		cout << endl << "sanePos: computation of map parameters" << endl;

	struct param_common dir; /* contains output input temp directories */
	struct samples samples_struct;  /*  everything about frames, noise files and frame processing order */

	struct param_saneProc Proc_param; /* contains user options about preprocessing properties */
	struct param_sanePos Pos_param; /* contains user options about map projection and properties */

	// struct used in the parser, but not by this program
	struct param_sanePS   PS_param;
	struct param_saneInv Inv_param;
	struct param_sanePic Pic_param;
	string parser_output="";


	uint32_t mask_sanePos = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUTPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | PIXDEG_WRONG_VALUE | FILEFORMAT_NOT_FOUND | NAPOD_WRONG_VALUE |
			FHP_PROBLEM | FITS_FILELIST_NOT_FOUND | FCUT_PROBLEM; // 0xc2ff


	// The next two variable are int because openmpi does not know about bool
	int flagon = 0; /* if rejectsample [ii]==3, flagon=1*/
	int pixout = 0; /* indicates that at least one pixel has been flagged and is out */

	long long npix; /* npix = number of filled pixels */
	long long npixsrc; /* number of pixels included in Crossing Constraint Removal */
	long long addnpix=0; /* add a number 'n' of pixels to the map */

	int nwcs=1;                      // We will only deal with one wcs, but needed afterwards....
	struct wcsprm * wcs;             // wcs structure of the image
	long NAXIS = 2, NAXIS1, NAXIS2;  // size of the image

	char * subheader;               // Additionnal header keywords
	int nsubkeys;                   //

	// System should be IEEE 754 complient (TODO : add in the doc)
	double lon_min=NAN, lon_max=NAN, lat_min=NAN, lat_max=NAN; /* min/max coordinates of the map*/
	double glon_min, glon_max, glat_min, glat_max; /* global  min and max (to get the min and max of all min/max computed by different processors) */

	// positions variables
	short *mask;
	long long *indpix, *indpsrc; /* pixels indices, CCR mask pixels indices */

	long long *pixon; /* this array is used to store the rules for pixels : they are seen or not */
	long long *pixon_tot=NULL;

	//naiv map data params
	long ns; /* number of samples for this scan, first frame number of this scan*/
	double fhp_pix = 0; //, fcut_pix = 0; /* frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples*/
	double *PNdNaiv, *PNdtotNaiv=NULL; /*  projected noised data, and global Pnd for mpi utilization */
	long *hitsNaiv, *hitstotNaiv=NULL; /* naivmap parameters : hits count */

	// -----------------------------------------------------------------------------//

	uint32_t parsed=0x0000; // parser error status
	uint32_t compare_to_mask; // parser error status

	if (argc<2) {/* not enough argument */
		if (rank == 0)
			cerr << "EE - Please run  " << StringOf(argv[0]) << " with a .ini file" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);
	}  else {

		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, Pos_param, Proc_param, PS_param, Inv_param, Pic_param, size, rank);

		compare_to_mask = parsed & mask_sanePos;

		// print parser warning and/or errors
		if (rank==0)
			cout << endl << parser_output << endl;

		if(compare_to_mask>0x0000){

			switch (compare_to_mask){/* error during parsing phase */

			case 0x0001:
				if (rank==0)
					cerr << " EE - Please run " << StringOf(argv[0]) << " using a correct *.ini file" << endl;
				break;

			default :
				if (rank==0)
					cerr << "EE - Wrong program options or argument. Exiting ! " <<  "("<< hex << compare_to_mask << ")" << endl;
				break;

			}

#ifdef USE_MPI
			MPI_Finalize();
#endif
			return EX_CONFIG;
		}
	}

	// -----------------------------------------------------------------------------//
#ifdef DEBUG
	time_t t2=time(NULL);
#endif

	/* ------------------------------------------------------------------------------------*/
	// Start...

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_param);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if(configureMPI(dir.output_dir, samples_struct, rank, size,
			bolo_rank,  bolo_size, node_rank, node_size,
			MPI_COMM_NODE, MPI_COMM_MASTER_NODE)){
		if (rank==0)
			cerr << endl << endl << "Exiting..." << endl;

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_CONFIG;
	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	if ( cleanup_dirfile_sanePos(dir.tmp_dir, samples_struct, bolo_rank) ) {
		cerr << "EE - Error in initializing dirfile" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill bolo_rank 0 has created dirfile architecture.
#endif

	// Read file size once for all
	if (bolo_rank == 0) {
		if ( readFramesFromDirfile(dir.tmp_dir, samples_struct)) {
			cerr << "EE - Error in frame size - Did you run sanePre ?" << endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, EX_CONFIG);
#endif
			return EX_CONFIG;
		}
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
	std::vector<long> nsamples_tot;
	nsamples_tot.assign(samples_struct.ntotscan, 0);

	MPI_Allreduce(&(samples_struct.nsamples[0]), &nsamples_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);

	samples_struct.nsamples = nsamples_tot;
	nsamples_tot.clear();

	//	MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_SUB);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Comm_split(MPI_COMM_WORLD, bolo_rank, 0, &MPI_COMM_MASTER_BOLO);

#endif

	//	read pointing header
	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, &subheader, &nsubkeys, rank)){ // read keyrec file
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
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

	if (Pos_param.eq2gal || Pos_param.gal2eq) {
		if (rank == 0)
			cout << endl<< "Converting coordinates..." << endl;

		if (bolo_rank == 0) {
			// TODO: This routine could also be parallelized by frame or bolo...
			if (convert_Dirfile_LON_LAT(samples_struct, Pos_param) ) {
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return 1;
			}
		}

		if (Pos_param.eq2gal)
			Pos_param.axistype = "GAL";
		if (Pos_param.gal2eq)
			Pos_param.axistype = "EQ";

	}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (Pos_param.maskfile == ""){

		//TODO Find an easy way to provide for a map header (without going through mask...)

		if(rank==0)
			cout << endl << "Determining Map Parameters..." << endl;

		// Populate the wcs structure ...
		// ... add pixel size in deg ...
		for (int ii = 0; ii < NAXIS; ii++) wcs->cdelt[ii] = (ii) ? Pos_param.pixdeg: -1*Pos_param.pixdeg ;
		for (int ii = 0; ii < NAXIS; ii++) strcpy(wcs->cunit[ii], "deg");

		// ... add axis label ...
		if (Pos_param.axistype.compare("EQ") == 0){
			char TYPE[2][5]  = { "RA--", "DEC-"};
			char NAME[2][16] = {"Right Ascension","Declination"};

			for (int ii = 0; ii < NAXIS; ii++) {
				strcpy(wcs->ctype[ii], &TYPE[ii][0]);
				strncat(wcs->ctype[ii],"-",1);
				strncat(wcs->ctype[ii],Pos_param.projcode.c_str(), 3);
				strcpy(wcs->cname[ii], &NAME[ii][0]);
			}
		}

		if (Pos_param.axistype.compare("GAL") == 0){
			char TYPE[2][5]  = { "GLON", "GLAT"};
			char NAME[2][19] = {"Galactic Longitude", "Galactic Latitude"};

			for (int ii = 0; ii < NAXIS; ii++) {
				strcpy(wcs->ctype[ii], &TYPE[ii][0]);
				strncat(wcs->ctype[ii],"-",1);
				strncat(wcs->ctype[ii],Pos_param.projcode.c_str(), 3);
				strcpy(wcs->cname[ii], &NAME[ii][0]);
			}
		}

		// Define a projection center ...
		// ... From the ini file ...
		if ( ! isnan(Pos_param.lon) && ! isnan(Pos_param.lat) ) {
			wcs->crval[0] = Pos_param.lon;
			wcs->crval[1] = Pos_param.lat;
		} else {
			// ... or the first valid data point
			double lon_center, lat_center;

			if(rank==0) {

				double *lon, *lat;
				int *flag;
				long index = samples_struct.iframe_min;
				long ns = samples_struct.nsamples[index];

				lon  = new double[ns];
				lat  = new double[ns];
				flag = new int[ns];
				if(readLonFromDirfile(samples_struct.dirfile_pointers[index], samples_struct.basevect[index], samples_struct.bolo_list[index][0], lon, ns))
					return 1;
				if(readLatFromDirfile(samples_struct.dirfile_pointers[index], samples_struct.basevect[index], samples_struct.bolo_list[index][0], lat, ns))
					return 1;
				if(readFlagFromDirfile(samples_struct.dirfile_pointers[index], samples_struct.basevect[index], samples_struct.bolo_list[index][0], flag, ns))
					return 1;

				int ii = 0;
				while (flag[ii] != 0){ ii++; }
				lon_center = lon[ii];
				lat_center = lat[ii];

				delete [] lon;
				delete [] lat;
				delete [] flag;

			}

			// cout << rank << " there2 " << lon_center << " " << lat_center << endl;

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&lon_center,  1,MPI_DOUBLE,0,MPI_COMM_WORLD);
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
			//	TODO: This routine could also be parallelize by frame or bolo...
			//        but the end MPI operation is a MPI_MIN or MPI_MAX
			//        so you need to insure that all bolo_rank receive work, which is not necessarily the case
			//        if you have very few bolo or several cpus... OR
			//        one could create a temporary MPI_Comm of the rank who actually did something on this round,
			//        but then you need to find a way to MPI_Send/MPI_Reicv the result from the temporary node 0 to all....
			if (bolo_rank == 0){
				if(computeMapMinima_HIPE(samples_struct, wcs, lon_min,lon_max,lat_min,lat_max)){
#ifdef USE_MPI
					MPI_Abort(MPI_COMM_WORLD, 1);
#endif
					return(EX_OSERR);
				}

			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);

			// Find the minimum for the values found by rank0
			if (bolo_rank == 0){
				MPI_Allreduce(&lon_min, &glon_min, 1, MPI_DOUBLE,MPI_MIN,MPI_COMM_MASTER_BOLO);
				MPI_Allreduce(&lon_max, &glon_max, 1, MPI_DOUBLE,MPI_MAX,MPI_COMM_MASTER_BOLO);
				MPI_Allreduce(&lat_min, &glat_min, 1, MPI_DOUBLE,MPI_MIN,MPI_COMM_MASTER_BOLO);
				MPI_Allreduce(&lat_max, &glat_max, 1, MPI_DOUBLE,MPI_MAX,MPI_COMM_MASTER_BOLO);
			}

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&glon_min,1,MPI_DOUBLE,0,MPI_COMM_NODE);
			MPI_Bcast(&glon_max,1,MPI_DOUBLE,0,MPI_COMM_NODE);
			MPI_Bcast(&glat_min,1,MPI_DOUBLE,0,MPI_COMM_NODE);
			MPI_Bcast(&glat_max,1,MPI_DOUBLE,0,MPI_COMM_NODE);
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
				cout << "     lon = " << lon_mean << endl;
				cout << "     lat = " << lat_mean << endl << endl;
			}
		}


		// We now have a good projection center, so
		//  ... set the structure to have the celestial projection routines
		if ( int status = wcsset(wcs) ) {
			printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
		}

		// At this stage all should be in DIRFILE so ... HIPE format
		if(computeMapMinima_HIPE(samples_struct, wcs, lon_min,lon_max,lat_min,lat_max)){
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_OSERR);
		}


#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		// Find the minimum for the values found by rank0
		if (bolo_rank == 0){
			MPI_Allreduce(&lon_min, &glon_min, 1, MPI_DOUBLE,MPI_MIN,MPI_COMM_MASTER_BOLO);
			MPI_Allreduce(&lon_max, &glon_max, 1, MPI_DOUBLE,MPI_MAX,MPI_COMM_MASTER_BOLO);
			MPI_Allreduce(&lat_min, &glat_min, 1, MPI_DOUBLE,MPI_MIN,MPI_COMM_MASTER_BOLO);
			MPI_Allreduce(&lat_max, &glat_max, 1, MPI_DOUBLE,MPI_MAX,MPI_COMM_MASTER_BOLO);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&glon_min,1,MPI_DOUBLE,0,MPI_COMM_NODE);
		MPI_Bcast(&glon_max,1,MPI_DOUBLE,0,MPI_COMM_NODE);
		MPI_Bcast(&glat_min,1,MPI_DOUBLE,0,MPI_COMM_NODE);
		MPI_Bcast(&glat_max,1,MPI_DOUBLE,0,MPI_COMM_NODE);

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

		//		int margingPixel = 0;  // for cosmetic...
		NAXIS1 = ceil(glon_max/Pos_param.pixdeg)-floor(glon_min/Pos_param.pixdeg)+1; // +2*margingPixel;
		NAXIS2 = ceil(glat_max/Pos_param.pixdeg)-floor(glat_min/Pos_param.pixdeg)+1; // +2*margingPixel;

		wcs->crpix[0] = ceil(    glon_max/Pos_param.pixdeg) + 1; // + margingPixel;
		wcs->crpix[1] = floor(-1*glat_min/Pos_param.pixdeg) + 1; // + margingPixel;

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
			cout << endl << "Reading Mask map : " << Pos_param.maskfile << endl;

		// TODO: Mask read by all rank... bad idea....
		if (read_mask_wcs(dir.input_dir + Pos_param.maskfile, "mask", wcs, NAXIS1, NAXIS2, mask )){
			cerr << "Error Reading Mask file" << endl;
#ifdef USE_MPI
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
		//		if(computeMaskIntersection(samples_struct, Pos_param, iframe_min, iframe_max, rank,
		//				wcs, NAXIS1, NAXIS2, mask, indpsrc, npixsrc, bolo_list, count_pixsrc)){
		//#ifdef USE_MPI
		//			MPI_Abort(MPI_COMM_WORLD, 1);
		//#endif
		//			return(EX_OSERR);
		//		}
		//
		//#ifdef USE_MPI
		//	if(rank==0){
		//		count_pixsrc_tot = new long long[npixsrc*samples_struct.ntotscan];
		//	}
		//	MPI_Reduce(count_pixsrc,count_pixsrc_tot,npixsrc*samples_struct.ntotscan,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		//#else
		//	count_pixsrc_tot=count_pixsrc;
		//#endif

	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// TODO : test and run with real data !!!
	//	if(modify_mask_flag_in_dirfile(dir.tmp_dir, samples_struct, bolo_list, indpsrc,
	//			NAXIS1, NAXIS2, iframe_min, iframe_max)){
	//		cout << "ERROR in  modify_mask_flag_in_dirfile... Exiting...\n";
	//#ifdef USE_MPI
	//		MPI_Abort(MPI_COMM_WORLD, 1);
	//#endif
	//		return(EX_CANTCREAT);
	//	}

	if (rank == 0) {
		printf("Map Size         : %ld x %ld pixels\n", NAXIS1, NAXIS2);

		if(save_keyrec(dir.tmp_dir,wcs, NAXIS1, NAXIS2, subheader, nsubkeys)){
#ifdef USE_MPI
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
	(Pos_param.flgdupl) ? factdupl = 2: factdupl = 1; //  default 1 : if flagged data are put in a duplicated map

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

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(computePixelIndex_HIPE(samples_struct, Proc_param, Pos_param,
			wcs, NAXIS1, NAXIS2,
			mask,factdupl, addnpix,  indpsrc, npixsrc,
			pixon, flagon, pixout,
			bolo_rank, bolo_size)){
#ifdef USE_MPI
		MPI_Finalize();
#endif
		return(EX_OSERR);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){
		pixon_tot = new long long[sky_size];
		fill(pixon_tot,pixon_tot+sky_size,0);
	}


	MPI_Allreduce(&flagon, &flagon,    1,        MPI_INT,    MPI_LOR, MPI_COMM_WORLD);
	MPI_Allreduce(&pixout, &pixout,    1,        MPI_INT,    MPI_LOR, MPI_COMM_WORLD);

	MPI_Reduce(pixon,  pixon_tot, sky_size, MPI_LONG_LONG,MPI_SUM, 0, MPI_COMM_WORLD);

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
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}
		if(writeIndexCCR((long long) NAXIS1*NAXIS2, npixsrc, indpsrc,  dir.tmp_dir)){
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}

		//		printf("Total number of Scans : %d \n", (int) samples_struct.ntotscan);
		cout << "Pixels Indices   : " << fixed << sky_size << " (" << prettyPrintSize(sky_size*8.) << ")" << endl;
		cout << "Filled Indices   : " << fixed << npix << " (" << prettyPrintSize(npix*8.) << ")" << endl;

		if (pixout != 0) {
			cout << endl;
			cout << "WW - There are samples outside the map limits : " << endl;
			cout << "     Assuming constant sky emission for those samples" << endl;
			cout << "     They are put in a single pixel" << endl;
		}
	}


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast( &npix,       1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(indpix,sky_size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

	//--------------------------------------- NAIVE MAP COMPUTATION ------------------------------------------//

	if(rank==0)
		cout << endl << "Computing Naive map..." << endl;

	// (At N-1 D) memory allocation
	PNdNaiv= new double[npix];
	hitsNaiv=new long[npix];

	fill(PNdNaiv,PNdNaiv+npix,0.0);
	fill(hitsNaiv,hitsNaiv+npix,0);



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// loop over the scans
	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){

		ns       = samples_struct.nsamples[iframe]; // number of samples for this scan
		fhp_pix  = samples_struct.fhp[iframe]  * double(ns)/samples_struct.fsamp[iframe]; // knee freq of the filter in terms of samples in order to compute fft
		//			fcut_pix = samples_struct.fcut[iframe] * double(ns)/samples_struct.fsamp[iframe]; // noise PS threshold freq, in terms of samples

		std::vector<string> det_vect = samples_struct.bolo_list[iframe];
		long ndet = (long)det_vect.size();

		int pb=0;
		pb+=do_PtNd_Naiv(samples_struct, Proc_param ,PNdNaiv, dir.tmp_dir, det_vect, ndet,
				Proc_param.poly_order, Proc_param.napod, fhp_pix, ns, indpix, iframe, hitsNaiv,
				bolo_rank, bolo_size);
		if(pb>0){
			cerr << "EE - Problem after do_PtNd_Naiv. Exiting..." << endl;
#ifdef USE_MPI
			MPI_Finalize();
#endif
			return EXIT_FAILURE;
		}

	} // end of iframe loop


#ifdef USE_MPI

	if(bolo_rank==0){	// global (At N-1 D) malloc for mpi
		PNdtotNaiv = new double[npix];
		hitstotNaiv=new long[npix];
	}

	MPI_Reduce(PNdNaiv,  PNdtotNaiv,  npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(hitsNaiv, hitstotNaiv, npix, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);
#else
	PNdtotNaiv  = PNdNaiv;
	hitstotNaiv = hitsNaiv;
#endif

	if (rank == 0){

		string fname = dir.output_dir + "naivMap.fits";

		cout << "Output file      : naivMap.fits" << endl;

		double *S;
		S = new double[npix];
		fill(S,S+npix,0.0);

		for (long ii=0; ii<npix; ii++) {
			if (hitstotNaiv[ii] != 0)
				S[ii] = PNdtotNaiv[ii]/(double)hitstotNaiv[ii];
		}


		// Write the raw maps (without flag and without crossing constrains)
		if (writeRawMapToFits(fname, S, NAXIS1, NAXIS2, indpix,	wcs, subheader, nsubkeys, false))
			cerr << "EE - Error writing Naive map... " << endl;

		if (writeHitMapToFits(fname, hitstotNaiv, NULL, addnpix, NAXIS1, NAXIS2,
				 indpix, indpsrc, npixsrc, factdupl,samples_struct.ntotscan,
				 wcs, subheader, nsubkeys, true) )
			cerr << "WW - Error writing Hit Map" << endl;

		if( write_fits_inifile(fname, dir, Proc_param, Pos_param, PS_param, Pic_param, Inv_param)) // write saneProc parameters in naive Map fits file header
			cerr << "WW - No ini file will be included in the file : " << fname << endl;
		if( write_fits_inputfile(fname, samples_struct))
			cerr << "WW - No input files will be included in the file : " << fname << endl;


		if (Pos_param.maskfile != "")
			if(copy_fits_mask(fname, dir.input_dir + Pos_param.maskfile)) // copy mask in naive map file
				cerr << "Warning ! The mask will not be included in naive map fits file ...\n";


	}
	/* ---------------------------------------------------------------------------------------------*/

	// clean up
	delete [] PNdNaiv;
	delete [] hitsNaiv;
	delete [] indpix;
	delete [] indpsrc;
	delete [] mask;
	delete [] pixon;


	wcsvfree(&nwcs, &wcs);

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// Close previously openened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])){
				cerr << "EE - error closing dirfile...";
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}


#ifdef DEBUG
	time_t t3=time(NULL);
	cout << "Total Time : " << t3-t2 << " sec\n";
#endif


#ifdef USE_MPI
	if(bolo_rank==0 && size > 1){
		// clean up
		delete [] PNdtotNaiv;
		delete [] hitstotNaiv;
		delete [] pixon_tot;
	}
	MPI_Comm_free(&MPI_COMM_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_BOLO);

	MPI_Finalize();
#endif

	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;
}
