#include <iostream>
#include <string>

#include "imageIO.h"
#include "struct_definition.h"

#include <sstream>

extern "C" {
#include <fitsio.h>
#include <nrutil.h>
#include "wcslib/wcslib.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}



using namespace std;


void print_fits_error(int status){
	if(status){
		fits_report_error(stderr, status); /* print error report */
		exit(status);    /* terminate the program, returning error status */
	}
	return;
}


//TODO : check and optimize
void write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data, string table_name ,bool fits_already_exist)
{
	// all angles in degrees

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...

	int naxis = 2;                  // number of dimensions
	long naxes[] = {NAXIS1, NAXIS2}; // size of dimensions
	long fpixel[] = {1, 1};          // index for write_pix
	long long ndata = NAXIS1 * NAXIS2;            // number of data points

	char *header, *hptr;
	int nkeyrec;

	if(fits_already_exist){
		if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status))
			fits_report_error(stderr, fits_status);
	}else{
		// create fits file
		if ( fits_create_file(&fp, fname.c_str(), &fits_status) )
			print_fits_error(fits_status);
	}

	// create fits image (switch on data type)
	switch (dtype) {
	case 'd':    // double
		if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
			print_fits_error(fits_status);
		break;
	case 'l':    // long
		if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) )
			print_fits_error(fits_status);
		break;
	default:
		printf("write_fits: data type %c not supported. Exiting.\n",dtype);
		exit(1);
	}

	// Transform wcsprm struture to header
	if ( (fits_status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) ){
		printf("wcshdo ERROR %d: %s.\n", fits_status, wcs_errmsg[fits_status]);
		exit(fits_status);
	}

	hptr = header;
	// write it to the fits file
	for (int keyrec = 0; keyrec < nkeyrec; keyrec++, hptr += 80){
		if ( fits_write_record(fp, (const char*) hptr, &fits_status))
			print_fits_error(fits_status);
	}


	fits_update_key(fp, TSTRING, (char *)"EXTNAME", (void*)(table_name.c_str()),
			(char *) "table name", &fits_status);

	free(header);

	// write date to file
	if ( fits_write_date(fp, &fits_status) )
		print_fits_error(fits_status);

	// write map data
	switch (dtype) {
	case 'd':    // double
		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data, &fits_status) )
			print_fits_error(fits_status);
		break;
	case 'l':    // long
		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status) )
			print_fits_error(fits_status);
		break;
	}


	// close file
	if(fits_close_file(fp, &fits_status))
		print_fits_error(fits_status);
}



void write_fits_hitory(string fname,long NAXIS1, long NAXIS2, struct param_process proc_param, struct param_positions pos_param, std::vector<double> fcut, struct detectors det, struct samples samples_struct, long ncomp){

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	std::ostringstream oss;

	//	int naxis = 2;                  // number of dimensions
	//	long naxes[] = {NAXIS1, NAXIS2};

	if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status))
		fits_report_error(stderr, fits_status);

	fits_create_img(fp, 8, 0, 0, &fits_status);

	fits_update_key(fp, TSTRING, (char *)"EXTNAME", (void*)"History",
			(char *) "table name", &fits_status);

	for(int num=0;num<(int)samples_struct.ntotscan;num++){
		oss << "Source" << num;
		string keyname = oss.str();
		string value = samples_struct.fits_table[num];
		string comm = "Data Source Fits Files";
		if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
			print_fits_error(fits_status);
		oss.str("");
	}


	oss << proc_param.napod;
	string keyname = "NAPOD";
	string value = oss.str();
	string comm = "Number of samples to apodize";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	oss << proc_param.poly_order;
	keyname = "POLYORDER";
	value = oss.str();
	comm = "Fitted polynomia order";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	oss << proc_param.fsamp;
	keyname = "SAMPLINGFREQUENCY";
	value = oss.str();
	comm = "sampling frequency (Hz)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	oss << proc_param.f_lp;
	keyname = "FILTERFREQUENCY";
	value = oss.str();
	comm = "Butterworth filter frequency (Hz)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);


	keyname = "FILLGAPS";
	if(proc_param.NOFILLGAP)
		value = "no";
	else
		value = "yes";
	comm = "Do we fill the gaps in the timeline with White noise + baseline ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "NORMLIN";
	if(proc_param.NORMLIN)
		value = "no";
	else
		value = "yes";
	comm = "Do we remove a baseline from the data ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "CORRELATION";
	if(proc_param.CORRon)
		value = "yes";
	else
		value = "no";
	comm = "Correlations between detectors are not included in the analysis ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "POLYNOMIASUBTRACTION";
	if(proc_param.remove_polynomia)
		value = "yes";
	else
		value = "no";
	comm = "Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	oss << pos_param.pixdeg;
	keyname = "PIXELSIZE";
	value = oss.str();
	comm = "SIZE OF THE PIXEL (deg)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "DUPLICATEDMAP";
	if(pos_param.flgdupl)
		value = "yes";
	else
		value = "no";
	comm = "flagged data are put in a separate map";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "GAPSPROJECTION";
	if(pos_param.projgaps)
		value = "yes";
	else
		value = "no";
	comm = "gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "SOURCESFILEFORMAT";
	if(pos_param.fileFormat)
		value = "HIPE";
	else
		value = "SANEPIC";
	comm = "Sources file format";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	keyname = "MASKFILE";
	oss << pos_param.maskfile;
	value = oss.str();
	comm = "name of the fits file that was used to mask radiant sources";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);


	oss << det.ndet;
	keyname = "NUMBEROFDETECTORS";
	value = oss.str();
	comm = "Number of detectors that were used for the analysis";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
		print_fits_error(fits_status);

	if(ncomp>0){
		oss << ncomp;
		keyname = "COMPONENTNUMBER";
		value = oss.str();
		comm = "number of noise component to estimate in sanePS";
		oss.str("");

		if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status))
			print_fits_error(fits_status);
	}


	// close file
	if(fits_close_file(fp, &fits_status))
		print_fits_error(fits_status);
}


int read_mask_wcs(string fname, string extname, /* char dtype,*/ struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2,  short *& data)
/*
 * Read the extension 'extname' from the 'fname' fits file, extension must be a 2D image
 * Return a wcs structure, the size of the image, the image data type and image itself cast to int, float or double
 */
{
	fitsfile *fptr;
	int status = 0, anynul, wcsstatus[NWCSFIX];
	char *header;
	int nkeyrec, nwcs, nreject;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	// Open the fits file...
	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ... and move to the 'extname'
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) extname.c_str(), NULL, &status))
		fits_report_error(stderr, status);

	// Retrieve the image size
	if (fits_get_img_size(fptr, 2, naxes, &status))
		fits_report_error(stderr, status);

	NAXIS1 = naxes[0];
	NAXIS2 = naxes[1];
	// Allocate the image container and read its depending on the type
	//	switch (dtype) {
	//	case 's':
	data = new short[NAXIS1*NAXIS2];
	if (fits_read_pix(fptr, TSHORT, fpixel, (long long) NAXIS1*NAXIS2, 0, data, &anynul, &status))
		fits_report_error(stderr, status);
	//		break;
	//	case 'i':
	//		data = new int[NAXIS1*NAXIS2];
	//		if (fits_read_pix(fptr, TINT, fpixel, (long long) NAXIS1*NAXIS2, 0, data, &anynul, &status))
	//			fits_report_error(stderr, status);
	//		break;
	//	case 'f':
	//		data = new float[NAXIS1*NAXIS2];
	//		if (fits_read_pix(fptr, TFLOAT, fpixel, (long long) NAXIS1*NAXIS2, 0, data, &anynul, &status))
	//			fits_report_error(stderr, status);
	//		break;
	//	case 'd':
	//		data = new double[NAXIS1*NAXIS2];
	//		if (fits_read_pix(fptr, TDOUBLE, fpixel, (long long) NAXIS1*NAXIS2, 0, data, &anynul, &status))
	//			fits_report_error(stderr, status);
	//		break;
	//	default:
	//		print_fits_error(BAD_DATATYPE);
	//	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status))
		fits_report_error(stderr, status);

	if (( status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status, wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status= wcsset(wcs))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}

	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	free(header);
	return(0);
}


//void write_fits_naivmap(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data, char maptype)
//{
//	// all angles in degrees
//
//	fitsfile *fp;
//	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
//
//	int naxis = 2;                  // number of dimensions
//	long naxes[] = {NAXIS1, NAXIS2}; // size of dimensions
//	long fpixel[] = {1, 1};          // index for write_pix
//	long long ndata = NAXIS1 * NAXIS2;            // number of data points
//
//	char *header, *hptr;
//	int nkeyrec;
//
//	if(maptype=='n')
//		write_fits_wcs(fnaivname, wcs, NAXIS1, NAXIS2, dtype, map1d);
//	else
//		if((maptype=='h')||(maptype=='v')||(maptype==i)){
//
//			if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
//				fits_report_error(stderr, status);
//
//			// create fits image (switch on data type)
//			switch (dtype) {
//			case 'd':    // double
//				if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
//					print_fits_error(fits_status);
//				break;
//			case 'l':    // long
//				if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) )
//					print_fits_error(fits_status);
//				break;
//			default:
//				printf("write_fits: data type %c not supported. Exiting.\n",dtype);
//				exit(1);
//			}
//
//			// Transform wcsprm struture to header
//			if ( (fits_status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) ){
//				printf("wcshdo ERROR %d: %s.\n", fits_status, wcs_errmsg[fits_status]);
//				exit(fits_status);
//			}
//
//			hptr = header;
//			// write it to the fits file
//			for (int keyrec = 0; keyrec < nkeyrec; keyrec++, hptr += 80)
//				if ( fits_write_record(fp, (const char*) hptr, &fits_status))
//					print_fits_error(fits_status);
//
//			free(header);
//
//			// write date to file
//			if ( fits_write_date(fp, &fits_status) )
//				print_fits_error(fits_status);
//		}
//
//	// write map data
//	switch (dtype) {
//	case 'd':    // double
//		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data, &fits_status) )
//			print_fits_error(fits_status);
//		break;
//	case 'l':    // long
//		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status) )
//			print_fits_error(fits_status);
//		break;
//	}
//
//	// close file
//	if(fits_close_file(fp, &fits_status))
//		print_fits_error(fits_status);
//}



//TODO : This function should be more generalized, now make the assumption that indpix has already been computed on the map
void read_fits_signal(string fname, double *S, long long* indpix, long &NAXIS1, long &NAXIS2, long long npix)
/*
 * This function read the sanePic generated map and converts it into S (only seen pixels)
 */
{
	fitsfile *fptr;
	int status = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long mi;
	double *map;


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	//	if (fits_movabs_hdu(fptr, 1, NULL, &status))
	//		fits_report_error(stderr, status);

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*)"Map", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 2, naxes, &status);

	NAXIS1=(long)naxes[0];
	NAXIS2=(long)naxes[1];

	// Initialize the data container
	map = new double [NAXIS1*NAXIS2];

	if (fits_read_pix(fptr, TDOUBLE, fpixel, (long long) NAXIS1*NAXIS2, 0, map, &anynul, &status))
		fits_report_error(stderr, status);


	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				S[indpix[mi]]= map[mi];
			}
		}
	}

	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}


void save_MapHeader(string outdir, struct wcsprm * wcs, long NAXIS1, long NAXIS2){

	FILE *fout;
	int nkeyrec, status;
	char *header, *hptr;

	if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		exit(0);
	}


	outdir=outdir + "mapHeader.keyrec";
	fout = fopen(outdir.c_str(),"w");
	if (fout==NULL) {fputs ("Creation error : File error on mapHeader.keyrec\n",stderr); exit (1);}

	fprintf(fout,"NAXIS1  = %20ld / %-47s\n",NAXIS1,"length of data axis 1");
	fprintf(fout,"NAXIS2  = %20ld / %-47s\n",NAXIS2,"length of data axis 2");

	hptr = header;
	for (int i = 0; i < nkeyrec; i++, hptr += 80) {
		fprintf(fout,"%.80s\n", hptr);
	}
	fclose(fout);
	free(header);
}

void print_MapHeader(struct wcsprm *wcs){

	int nkeyrec;
	char * header, *hptr ;
	if (int status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		exit(0);
	}
	hptr = header;
	printf("\n\n Map Header :\n");
	for (int ii = 0; ii < nkeyrec; ii++, hptr += 80) {
		printf("%.80s\n", hptr);
	}
	free(header);

}

void read_MapHeader(string outdir, struct wcsprm * & wcs, long * NAXIS1, long * NAXIS2){

	outdir = outdir + "mapHeader.keyrec";

	FILE *fin;
	char *memblock;
	int size, nkeyrec, nreject, nwcs, status;
	size_t result;

	fin = fopen(outdir.c_str(),"r");
	if (fin==NULL) {fputs ("Read error : File error on mapHeader.keyrec",stderr); exit (1);}

	fseek(fin, 0L, SEEK_END);     /* Position to end of file */
	size = ftell(fin);            /* Get file length */
	rewind(fin);                  /* Back to start of file */


	nkeyrec = size/81;

	char comment[47];

	// Read the two first lines, NAXIS1/NAXIS2
	result = fscanf(fin,"NAXIS1  = %20ld / %47c\n",NAXIS1,(char *) &comment);
	result = fscanf(fin,"NAXIS2  = %20ld / %47c\n",NAXIS2,(char *) &comment);

	memblock = new char [(nkeyrec-2)*80];
	for (int ii = 0; ii < nkeyrec; ii++) {
		result = fread(&memblock[ii*80], 80, sizeof(char), fin);
		fseek(fin, 1, SEEK_CUR); // skip newline char
	}
	fclose (fin);
	/* Parse the primary header of the FITS file. */
	/* -2 to handle the firts two NAXIS? keyword */
	if ((status = wcspih(memblock, nkeyrec-2, WCSHDR_all, 2, &nreject, &nwcs, &wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
	}
	delete[] memblock;
	//	free(comment);

}

