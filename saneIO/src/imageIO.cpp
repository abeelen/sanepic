#include <iostream>
#include <string>

#include "imageIO.h"
#include "struct_definition.h"
#include "inputFileIO.h"

#include <sstream>

extern "C" {
#include <fitsio.h>
#include <nrutil.h>
#include "wcslib/wcslib.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}



using namespace std;


int write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data, string table_name ,bool fits_already_exist)
{
	// all angles in degrees

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	int start=1;
	int naxis = 2;                  // number of dimensions
	long naxes[] = {NAXIS1, NAXIS2}; // size of dimensions
	long fpixel[] = {1, 1};          // index for write_pix
	long long ndata = NAXIS1 * NAXIS2;            // number of data points

	char *header, *hptr;
	int nkeyrec;

	if(fits_already_exist){
		if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}
	}else{
		// create fits file
		if ( fits_create_file(&fp, fname.c_str(), &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		start=0;
	}

	for(int ii=start;ii<2;ii++){
		naxis=2*ii;
		naxes[0]=NAXIS1*ii;
		naxes[1]=NAXIS2*ii;

		// create fits image (switch on data type)
		switch (dtype) {
		case 'd':    // double
			if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) ){
				fits_report_error(stderr, fits_status);
				return 1;
			}
			break;
		case 'l':    // long
			if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) ){
				fits_report_error(stderr, fits_status);
				return 1;
			}
			break;
		default:
			printf("write_fits: data type %c not supported. Exiting.\n",dtype);
			exit(1);
		}

		// Transform wcsprm struture to header
		if ( (fits_status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) ){
			printf("wcshdo ERROR %d: %s.\n", fits_status, wcs_errmsg[fits_status]);
			return 1;
		}

		// write date to file
		if ( fits_write_date(fp, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}

		hptr = header;
		// write it to the fits file
		for (int keyrec = 0; keyrec < nkeyrec; keyrec++, hptr += 80){
			if ( fits_write_record(fp, (const char*) hptr, &fits_status)){
				fits_report_error(stderr, fits_status);
				return 1;
			}
		}
		free(header);
	}

	if(fits_update_key(fp, TSTRING, (char *)"EXTNAME", (void*)(table_name.c_str()),
			(char *) "table name", &fits_status))
		return 1;


	// write map data
	switch (dtype) {
	case 'd':    // double
		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	case 'l':    // long
		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	}


	// close file
	if(fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;
}



int write_fits_hitory(string fname,long NAXIS1, long NAXIS2, string path, struct param_process proc_param, struct param_positions pos_param, std::vector<double> fcut, struct detectors det, struct samples samples_struct, long ncomp){

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	std::ostringstream oss;

	if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status))
		fits_report_error(stderr, fits_status);

	fits_create_img(fp, 8, 0, 0, &fits_status);

	fits_update_key(fp, TSTRING, (char *)"EXTNAME", (void*)"History",
			(char *) "table name", &fits_status);

	string keyname = "PATHNAME";
	string value = path;
	string comm = "Source data path";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	for(int num=0;num<(int)samples_struct.ntotscan;num++){
		oss << "Source" << num;
		string keyname = oss.str();
		oss.str("");
		string value = FitsBasename(samples_struct.fits_table[num]) + ".fits";
		string comm = "Data Source Fits File";
		if (fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}

	}


	oss << proc_param.napod;
	keyname = "NAPOD";
	value = oss.str();
	comm = "Number of samples to apodize";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	oss << proc_param.poly_order;
	keyname = "POLYORDR";
	value = oss.str();
	comm = "Fitted"; // polynomia order";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	// use only 8 characters for keyname to avoid HIERARCH keyword addition by cfitsio...
	oss << proc_param.fsamp;
	keyname = "SAMPFREQ";
	value = oss.str();
	comm = "sampling frequency (Hz)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	oss << proc_param.f_lp;
	keyname = "FILTFREQ";
	value = oss.str();
	comm = "Butterworth filter frequency (Hz)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	keyname = "FILLGAPS";
	if(proc_param.NOFILLGAP)
		value = "no";
	else
		value = "yes";
	comm = "Do we fill the gaps in the timeline with White noise + baseline ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "NORMLIN";
	if(proc_param.NORMLIN)
		value = "no";
	else
		value = "yes";
	comm = "Do we remove a baseline from the data ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "CORREL";
	if(proc_param.CORRon)
		value = "yes";
	else
		value = "no";
	comm = "Correlations between detectors are not included in the analysis ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "POLYNSUB";
	if(proc_param.remove_polynomia)
		value = "yes";
	else
		value = "no";
	comm = "Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	oss << pos_param.pixdeg;
	keyname = "PIXSIZE";
	value = oss.str();
	comm = "SIZE OF THE PIXEL (deg)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "DUPLMAP";
	if(pos_param.flgdupl)
		value = "yes";
	else
		value = "no";
	comm = "flagged data are put in a separate map";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "GAPSPROJ";
	if(pos_param.projgaps)
		value = "yes";
	else
		value = "no";
	comm = "gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "FORMAT";
	if(pos_param.fileFormat)
		value = "HIPE";
	else
		value = "SANEPIC";
	comm = "Sources file format";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "MASKFILE";
	oss << pos_param.maskfile;
	value = oss.str();
	comm = "name of the fits file that was used to mask radiant sources";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	oss << det.ndet;
	keyname = "NUMDET";
	value = oss.str();
	comm = "Number of detectors that were used for the analysis";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if(ncomp>0){
		oss << ncomp;
		keyname = "COMPONEN";
		value = oss.str();
		comm = "number of noise component to estimate in sanePS";
		oss.str("");

		if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}
	}


	// close file
	if(fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;
}


int write_fits_mask(std::string fnaivname, std::string maskfile)
{

	fitsfile *fptr, *outfptr;
	int fits_status = 0;
	long naxes[2] = { 1, 1 };

	if (fits_open_file(&fptr, maskfile.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_open_file(&outfptr, fnaivname.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	if(fits_copy_header(fptr, outfptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	// Retrieve the image size
	if (fits_get_img_size(fptr, 2, naxes, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	//	long NAXIS1 = naxes[0];

	if(fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "mask", NULL, &fits_status))


		//	for(long col=1;col<NAXIS1;col++)
		//		if(fits_copy_col(fptr, outfptr,  col, col,	0, &fits_status))
		//			fits_report_error(stderr, fits_status);
		if (fits_copy_data(fptr, outfptr, &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}


	// close files
	if(fits_close_file(fptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if(fits_close_file(outfptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

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
	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ... and move to the 'extname'
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) extname.c_str(), NULL, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// Retrieve the image size
	if (fits_get_img_size(fptr, 2, naxes, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	NAXIS1 = naxes[0];
	NAXIS2 = naxes[1];
	// Allocate the image container and read its depending on the type
	//	switch (dtype) {
	//	case 's':
	data = new short[NAXIS1*NAXIS2];
	if (fits_read_pix(fptr, TSHORT, fpixel, (long long) NAXIS1*NAXIS2, 0, data, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

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

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	free(header);
	return(0);
}


int read_fits_signal(string fname, double *S, long long* indpix, long &NAXIS1, long &NAXIS2, struct wcsprm * wcs)
/*
 * This function read the sanePic generated map and converts it into S (only seen pixels)
 */
{
	fitsfile *fptr;
	int status = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long mi;
	double *map;
	char comment[80];
	double crpix1=0, crpix2=0;
	double crval1=0.0, crval2=0.0;
	double cdelt1=0.0, cdelt2=0.0;
	double lonpole=0.0, latpole=0.0;


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// TODO: Change that...
	// - read the whole header
	// - and make a wcs struct out of it
	// - and write a routine to compare two wcs struct
	//
	// Read fits header keys
	if (fits_read_key(fptr,TDOUBLE, (char *) "CRPIX1", &crpix1, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "CRPIX2", &crpix2, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "CRVAL1", &crval1, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "CRVAL2", &crval2, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "CDELT1", &cdelt1, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "CDELT2", &cdelt2, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "LONPOLE", &lonpole, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_read_key(fptr,TDOUBLE, (char *) "LATPOLE", &latpole, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}


	// compatibility verifications !

	if(crpix1!=wcs->crpix[0] || crpix2!=wcs->crpix[1]){
		cout << "CRPIX are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(crval1!=wcs->crval[0] || crval2!=wcs->crval[1]){
		cout << "CRVAL are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(cdelt1!=wcs->cdelt[0] || cdelt2!=wcs->cdelt[1]){
		cout << "CDELT are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(lonpole!=wcs->lonpole || latpole!=wcs->latpole){
		cout << "LONPOLE and/or LATPOLE are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*)"Image", NULL, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if(fits_get_img_size(fptr, 2, naxes, &status))
		return 1;

	NAXIS1=(long)naxes[0];
	NAXIS2=(long)naxes[1];



	// Initialize the data container
	map = new double [NAXIS1*NAXIS2];

	if (fits_read_pix(fptr, TDOUBLE, fpixel, (long long) NAXIS1*NAXIS2, 0, map, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// seems to work correctly
	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				S[indpix[mi]]= map[mi];
			}
		}
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;
}


int save_MapHeader(string outdir, struct wcsprm * wcs, long NAXIS1, long NAXIS2){

	FILE *fout;
	int nkeyrec, status;
	char *header, *hptr;

	if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		exit(0);
	}


	outdir=outdir + "mapHeader.keyrec";
	fout = fopen(outdir.c_str(),"w");
	if (fout==NULL) {fputs ("Creation error : File error on mapHeader.keyrec\n",stderr); return (1);}

	fprintf(fout,"NAXIS1  = %20ld / %-47s\n",NAXIS1,"length of data axis 1");
	fprintf(fout,"NAXIS2  = %20ld / %-47s\n",NAXIS2,"length of data axis 2");

	hptr = header;
	for (int i = 0; i < nkeyrec; i++, hptr += 80) {
		fprintf(fout,"%.80s\n", hptr);
	}
	fclose(fout);
	free(header);

	return 0;
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

