#include <iostream>
#include <string>

#include "imageIO.h"


using namespace std;

extern "C" {
#include <fitsio.h>
#include <nrutil.h>
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}


void print_fits_error(int status){
	if(status){
		fits_report_error(stderr, status); /* print error report */
		exit(status);    /* terminate the program, returning error status */
	}
	return;
}


/*void write_fits(string fname, double pixsize, long nx, long ny,
		double *tancoord, double *tanpix, int coordsyst, char dtype, void *data)
{
	// all angles in degrees
	// coordcenter is a 2-element array containing RA/DEC (or l/b) of the central pixel

	fitsfile *fp;
	int fits_status = 0;

	long naxis = 2;           // number of dimensions
	long naxes[] = {nx, ny};  // size of dimensions
	long fpixel[] = {1, 1};   // index for write_pix
	long ndata = nx * ny;     // number of data points

	double dtmp;
	char *strx, *stry;


	// create fits file
	if ( fits_create_file(&fp, fname.c_str(), &fits_status) )
		print_fits_error(fits_status);

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

	// write date to file
	if ( fits_write_date(fp, &fits_status) )
		print_fits_error(fits_status);

	// write map parameters (keywords)
	if ( fits_write_key(fp, TLONG, (char*)"NROW", &nx, (char*)"Number of rows", &fits_status) )
		print_fits_error(fits_status);

	if ( fits_write_key(fp, TLONG, (char*)"NCOL", &ny, (char*)"Number of columns", &fits_status) )
		print_fits_error(fits_status);

	if ( fits_write_key(fp, TDOUBLE, (char*)"PIXSIZE", &pixsize, (char*)"Size of pixels (deg)", &fits_status) )
		print_fits_error(fits_status);

	if ( fits_write_comment(fp, "Galactic coordinates",  &fits_status) )
		print_fits_error(fits_status);

	dtmp = (tanpix[0]); // 0-based index to 1
	if ( fits_write_key(fp, TDOUBLE, (char*)"CRPIX1", &dtmp, (char*)"X PIXEL OF TANGENT POINT", &fits_status) )
		print_fits_error(fits_status);

	dtmp = (tanpix[1]); // 0-based index to 1
	if ( fits_write_key(fp, TDOUBLE, (char*)"CRPIX2", &dtmp, (char*)"Y PIXEL OF TANGENT POINT", &fits_status) )
		print_fits_error(fits_status);

	dtmp = -pixsize;
	if ( fits_write_key(fp, TDOUBLE, (char*)"CDELT1", &dtmp, (char*)"COORD VALUE INCR DEG/PIXEL AT ORIGIN ON LINE AXIS",
			&fits_status) )
		print_fits_error(fits_status);

	if ( fits_write_key(fp, TDOUBLE, (char*)"CDELT2", &pixsize, (char*)"COORD VALUE INCR DEG/PIXEL AT ORIGIN ON LINE AXIS",
			&fits_status) )
		print_fits_error(fits_status);

	if (coordsyst == 2){
		if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL1", tancoord, (char*)"GLON AT TANGENT POINT (DEG)", &fits_status) )
			print_fits_error(fits_status);

		if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL2", tancoord+1, (char*)"GLAT AT TANGENT POINT (DEG)", &fits_status) )
			print_fits_error(fits_status);

		strx = (char *)"GLON-TAN";
		stry = (char *)"GLAT-TAN";

	} else {
		if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL1", tancoord, (char*)"RA AT TANGENT POINT (DEG)", &fits_status) )
			print_fits_error(fits_status);

		if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL2", tancoord+1, (char*)"DEC AT TANGENT POINT (DEG)", &fits_status) )
			print_fits_error(fits_status);

		strx = (char *)"RA---TAN";
		stry = (char *)"DEC--TAN";

	}

	if ( fits_write_key(fp, TSTRING, (char*)"CTYPE1", strx, (char*)"TANGENT PLANE PROJECTION", &fits_status) )
		print_fits_error(fits_status);

	if ( fits_write_key(fp, TSTRING, (char*)"CTYPE2", stry, (char*)"TANGENT PLANE PROJECTION", &fits_status) )
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

}*/

//TODO : check and optimize
void write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data)
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
	// create fits file
	if ( fits_create_file(&fp, fname.c_str(), &fits_status) ) // TODO : add exit procedure
		print_fits_error(fits_status);

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
		fits_delete_file(fp, &fits_status);
		exit(fits_status);
	}

	hptr = header;
	// write it to the fits file
	for (int keyrec = 0; keyrec < nkeyrec; keyrec++, hptr += 80)
		if ( fits_write_record(fp, (const char*) hptr, &fits_status))
			print_fits_error(fits_status);

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

//TODO : This function should be more generalized
void read_fits_signal(string fname, double *S, long long* indpix, long &NAXIS1, long &NAXIS2, long long npix)
/*
 * This function read the sanePic generated map and converts it into S (only seen pixels)
 */
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long mi;
	double **map;


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List

	if (fits_movabs_hdu(fptr, 1, NULL, &status))
		fits_report_error(stderr, status);

	//fits_get_num_rows(fptr, &nn2, &status);

	//cout << "nn2 : " << nn2;

	fits_get_img_size(fptr, 2, naxes, &status);
	cout << naxes[0] << " " << naxes[1] << endl;

	NAXIS1=(long)naxes[0];
	NAXIS2=(long)naxes[1];

	// Initialize the data container
	map = dmatrix(0, NAXIS1 - 1, 0, NAXIS2 - 1);

	for (int i = 0; i < NAXIS1; i++) {
		fpixel[1] = i + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, NAXIS2, NULL, map[i], NULL, &status);
	}

//	cout << "map" << endl;
//	cout << map[100][100] << endl;

	int uu=1;

	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				S[indpix[mi]]= - map[jj][ii]; // TODO : suppress - // added minus mat 28_07
				uu++;
			}
		}
	}
	cout << npix << " " << uu << endl;


	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}


void save_MapHeader(string outdir, struct wcsprm wcs, long NAXIS1, long NAXIS2){

	FILE *fout;
	int nkeyrec, status;
	char *header, *hptr;

	if ((status = wcshdo(WCSHDO_all, &wcs, &nkeyrec, &header))) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		exit(0);
	}


	outdir=outdir + "mapHeader.keyrec";
	fout = fopen(outdir.c_str(),"w");
	if (fout==NULL) {fputs ("File error on mapHeader.keyrec",stderr); exit (1);}

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
	if (fin==NULL) {fputs ("File error on mapHeader.keyrec",stderr); exit (1);}

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

