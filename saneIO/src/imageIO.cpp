#include <iostream>
#include <string>

#include "imageIO.h"
#include "struct_definition.h"
#include "inputFileIO.h"
#include "utilities.h"

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
			return 1;
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

int write_fits_hitory2(std::string fname,long NAXIS1, long NAXIS2, string path, struct param_sanePre proc_param, struct param_sanePos pos_param, std::vector<double> fcut, struct samples samples_struct, long ncomp)
{

	fitsfile *fptr;
	int fits_status = 0;
	//	long naxes[2] = { 1, 1 };
	std::vector<string> key_vect;
	std::vector<string> comm_vect;
	std::vector<string> value_vect;

	//fill the key_vect
	key_vect.push_back("EXTNAME"); key_vect.push_back("PATHNAME");
	for(long iframe=0;iframe < samples_struct.ntotscan; iframe ++)
		key_vect.push_back("SOURCES");
	key_vect.push_back("NAPOD"); key_vect.push_back("POLY ORDER"); key_vect.push_back("SAMPLING FREQ");
	key_vect.push_back("FILTER FREQ"); key_vect.push_back("FILLGAPS"); key_vect.push_back("NORMLIN");
	key_vect.push_back("CORRELATION"); key_vect.push_back("POLY SUBTRACTION"); key_vect.push_back("PIXSIZE");
	key_vect.push_back("MAP DUPLICATION"); key_vect.push_back("GAPS PROJECTED"); key_vect.push_back("FILE FORMAT");
	key_vect.push_back("MASK FILE");
	if(ncomp>0)
		key_vect.push_back("NOISE COMPONENT");


	comm_vect.push_back("file name");
	comm_vect.push_back("Source data path");
	for(long iframe=0;iframe < samples_struct.ntotscan; iframe ++)
		comm_vect.push_back("Data Source Fits File");
	comm_vect.push_back("Number of samples to apodize");
	comm_vect.push_back("polynomia order");
	comm_vect.push_back("sampling frequency (Hz)");
	comm_vect.push_back("Butterworth filter frequency (Hz)");
	comm_vect.push_back("Do we fill the gaps in the timeline with White noise + baseline ?");
	comm_vect.push_back("Do we remove a baseline from the data ?");
	comm_vect.push_back("Correlations between detectors are not included in the analysis ?");
	comm_vect.push_back("Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment ?");
	comm_vect.push_back("SIZE OF THE PIXEL (deg)");
	comm_vect.push_back("flagged data are put in a separate map");
	comm_vect.push_back("gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively");
	comm_vect.push_back("Sources files format");
	comm_vect.push_back("name of the fits file that was used to mask radiant sources");
	if(ncomp>0)
		comm_vect.push_back("number of noise component estimated in sanePS");

	value_vect.push_back(fname);
	value_vect.push_back(path);
	for(int num=0;num<(int)samples_struct.ntotscan;num++){
		value_vect.push_back(FitsBasename(samples_struct.fitsvect[num]) + ".fits");
	}
	value_vect.push_back(StringOf(proc_param.napod));
	value_vect.push_back(StringOf(proc_param.poly_order));
	value_vect.push_back(StringOf(proc_param.fsamp));
	value_vect.push_back(StringOf(proc_param.f_lp));
	value_vect.push_back(StringOf(proc_param.NOFILLGAP ? "False" : "True"));
	value_vect.push_back(StringOf(proc_param.NORMLIN ? "False" : "True"));
	value_vect.push_back(StringOf(proc_param.CORRon ? "True" : "False"));
	value_vect.push_back(StringOf(proc_param.remove_polynomia ? "True" : "False"));
	value_vect.push_back(StringOf(pos_param.pixdeg));
	value_vect.push_back(StringOf(pos_param.flgdupl ? "True" : "False"));
	value_vect.push_back(StringOf(pos_param.projgaps ? "True" : "False"));
	value_vect.push_back(StringOf(pos_param.fileFormat ? "HIPE" : "SANEPIC"));
	value_vect.push_back(pos_param.maskfile);
	if(ncomp>0)
		value_vect.push_back(StringOf(ncomp));



	if (fits_open_file(&fptr, fname.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}
	// ---------------------------------------------
	// write the Channel List

	char *ttype[] = { (char*) "KEY", (char*) "values", (char*) "comments" };
	char *tform[] = { tableFormat(key_vect), tableFormat(value_vect), tableFormat(comm_vect) };
	char *tunit[] = { (char*) "None", (char*) "None", (char*) "None" };
	char **data, **data2, **data_value;

	data = vString2carray(key_vect);
	data2 = vString2carray(comm_vect);
	data_value = vString2carray(value_vect);

	if (fits_create_tbl(fptr, BINARY_TBL, key_vect.size(), 3, ttype, tform, tunit,
			(char*)"History", &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 1, 1, 1, key_vect.size(), data, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 2, 1, 1, value_vect.size(), data_value, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 3, 1, 1, comm_vect.size(), data2, &fits_status))
		return 1;
	if (fits_write_key(fptr, TSTRING, (char *) "TUNIT1", (char *) "NONE",
			(char *) "physical unit of the field", &fits_status))
		return 1;

	// close file
	if(fits_close_file(fptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	for(long ii=0;ii<(long)key_vect.size();ii++){
		delete [] data[ii];
		delete [] data2[ii];
		delete [] data_value[ii];
	}
	delete [] data2;
	delete [] data_value;
	delete [] data;
	delete [] *tform;

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


int read_fits_signal(string fname, double *S, long long* indpix, long NAXIS1, long NAXIS2, struct wcsprm * wcs)
/*
 * This function read the sanePic generated map and converts it into S (only seen pixels)
 */
{
	fitsfile *fptr;
	int status = 0, anynul, wcsstatus[NWCSFIX];
	char *header;
	int nkeyrec, nwcs, nreject;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long imNAXIS1, imNAXIS2;
	long mi;
	double *map;
	struct wcsprm *wcs_fits;

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*)"Image", NULL, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (( status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs_fits))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs_fits, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status, wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status= wcsset(wcs_fits))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}


	if(fits_get_img_size(fptr, 2, naxes, &status))
		return 1;

	imNAXIS1=(long)naxes[0];
	imNAXIS2=(long)naxes[1];

	if(compare_wcs(fname, wcs, wcs_fits, NAXIS1, NAXIS2, imNAXIS1, imNAXIS2))
		return 1;

	// Initialize the data container
	map = new double [imNAXIS1*imNAXIS2];

	if (fits_read_pix(fptr, TDOUBLE, fpixel, (long long) imNAXIS1*imNAXIS2, 0, map, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// seems to work correctly
	for (long ii=0; ii<imNAXIS1; ii++) {
		for (long jj=0; jj<imNAXIS2; jj++) {
			mi = jj*imNAXIS1 + ii;
			if (indpix[mi] >= 0){
				S[indpix[mi]]= map[mi];
			}
		}
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	//	nwcs=1;
	wcsvfree(&nwcs, &wcs_fits);
	free(header);
	return 0;
}


int save_keyrec(string outdir, struct wcsprm * wcs, long NAXIS1, long NAXIS2){

	FILE *fout;
	int nkeyrec, status;
	char *header, *hptr;

	if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		return 1;
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

int print_MapHeader(struct wcsprm *wcs){

	int nkeyrec;
	char * header, *hptr ;
	if (int status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		return 1;
	}
	hptr = header;
	printf("\n\n Map Header :\n");
	for (int ii = 0; ii < nkeyrec; ii++, hptr += 80) {
		printf("%.80s\n", hptr);
	}
	free(header);

	return 0;
}

void read_keyrec(string outdir, struct wcsprm * & wcs, long * NAXIS1, long * NAXIS2){

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


int compare_wcs(std::string fname, struct wcsprm *wcs, struct wcsprm *wcs_fits, long NAXIS1, long NAXIS2, long imNAXIS1, long imNAXIS2){

	// compatibility verifications !

	if(wcs_fits->crpix[0]!=wcs->crpix[0] || wcs_fits->crpix[1]!=wcs->crpix[1]){
		cout << "CRPIX are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(wcs_fits->crval[0]!=wcs->crval[0] || wcs_fits->crval[1]!=wcs->crval[1]){
		cout << "CRVAL are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(wcs_fits->cdelt[0]!=wcs->cdelt[0] || wcs_fits->cdelt[1]!=wcs->cdelt[1]){
		cout << "CDELT are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(wcs_fits->lonpole!=wcs->lonpole || wcs_fits->latpole!=wcs->latpole){
		cout << "LONPOLE and/or LATPOLE are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(NAXIS1!=imNAXIS1 || NAXIS2!=imNAXIS2){
		cout << "NAXIS are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}


	return 0;

}

