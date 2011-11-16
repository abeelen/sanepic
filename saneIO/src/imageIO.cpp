#ifdef HAVE_CONFIG_H
#include "../../config.h"
#else
#define PACKAGE_VERSION "Unknown"
#endif

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "struct_definition.h"
#include "inputFileIO.h"
#include "utilities.h"
#include "parser_functions.h" // for exporting the structures

extern "C" {
#include <fitsio.h>
#include "wcslib/wcslib.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "imageIO.h"

using namespace std;

int get_fits_META(string fname, struct wcsprm * &wcs, char ** subheader, int *nsubkeys, int rank) {

	int NAXIS = 2;
	if(rank==0){
		// Construct the wcsprm structure

		wcs = (struct wcsprm *) malloc(sizeof(struct wcsprm));
		wcs->flag = -1;
		wcsini(1, NAXIS, wcs);

		// A number of field are part of the wcsprm structure
		// http://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/structwcsprm.html


		fitsfile *fp;
		int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...

		double dvalue;
		char value[80], comment[80], card[162];
		char *subheadptr;

		if (fits_open_file(&fp, fname.c_str(), READONLY, &fits_status)) {
			fits_report_error(stderr, fits_status); return 1;
		}
		//
		//	// retrieve the wcs header ...
		//	char * header, *hptr;
		//	int nkeyrec;
		//	if (fits_hdr2str(fp, 1, NULL, 0, &header, &nkeyrec, &fits_status)) {
		//		fits_report_error(stderr, fits_status); return 1;
		//	}
		//	// .... print it
		//	hptr = header;
		//	printf("\n\n Header :\n");
		//	for (int ii = 0; ii < nkeyrec; ii++, hptr += 80) {
		//		printf("%.80s\n", hptr);
		//	}


		// As wcspih only deals with WCS coordinates when need to go around...


		if (fits_read_key_dbl(fp, (char *) "EQUINOX", &dvalue, comment, &fits_status)) {
			if (fits_status == KEY_NO_EXIST) { dvalue = 2000.0; fits_status = 0;
			} else { fits_report_error(stderr, fits_status); return 1; }
		}
		wcs->equinox = dvalue ;

		if (fits_read_key_str(fp, (char *) "DATE-OBS", value, comment, &fits_status)) {
			if (fits_status == KEY_NO_EXIST) { strcpy(value,"Unknown"); fits_status = 0;
			} else { fits_report_error(stderr, fits_status); return 1; }
		}
		strcpy(wcs->dateobs, value);

		if (fits_read_key_str(fp, (char *) "RADESYS", value, comment, &fits_status)) {
			if (fits_status == KEY_NO_EXIST) { strcpy(value,"Unknown"); fits_status = 0;
			} else { fits_report_error(stderr, fits_status); return 1; }
		}
		strcpy(wcs->radesys, value);

		if (fits_read_key_dbl(fp, (char *) "RESTWAV", &dvalue, comment,	&fits_status)) {
			if (fits_status == KEY_NO_EXIST) { dvalue = -1; fits_status = 0;
			} else { fits_report_error(stderr, fits_status); return 1; }
		}
		wcs->restwav = dvalue;

		if (fits_read_key_dbl(fp, (char *) "RESTFRQ", &dvalue, comment,	&fits_status)) {
			if (fits_status == KEY_NO_EXIST) { dvalue = -1; fits_status = 0;
			} else { fits_report_error(stderr, fits_status); return 1; 	}
		}
		wcs->restfrq = dvalue;

		/*
	// Removed feature, so that it does not add a problem in case of change of coordinates....

	if (fits_read_key_dbl(fp, (char *) "RA", &dvalue, comment,	&fits_status)) {
		if (fits_status == KEY_NO_EXIST) { dvalue = 361; fits_status = 0;
		} else { fits_report_error(stderr, fits_status); return 1; 	}
	}
	wcs->crval[0] = dvalue;

	if (fits_read_key_dbl(fp, (char *) "DEC", &dvalue, comment,	&fits_status)) {
		if (fits_status == KEY_NO_EXIST) { dvalue = 361; fits_status = 0;
		} else { fits_report_error(stderr, fits_status); return 1; 	}
	}
	wcs->crval[1] = dvalue;
		 */

		// Get the rest of the keys, not in the wcs struct....
		std::vector<string> keys;

		keys.push_back("TIMESYS");
		keys.push_back("CREATOR");
		keys.push_back("INSTRUME");
		keys.push_back("OBJECT");
		keys.push_back("TELESCOP");
		keys.push_back("OBSERVER");
		keys.push_back("UNIT");

		// Inspired by "ffhdr2str" from cfitsio

		*subheader = (char *) calloc ( (keys.size()+ 1) * 80 + 1, 1);
		*nsubkeys = 0 ;

		subheadptr = *subheader;

		for (unsigned int ii = 0; ii < keys.size(); ii++) {

			if ( fits_read_card(fp,  (char *) keys[ii].c_str() , card, &fits_status) ) {
				if (fits_status == KEY_NO_EXIST) {
					// If the key is not found, just drop it...
					fits_status = 0;
				} else {
					fits_report_error(stderr, fits_status);
					return 1;
				}
			} else {
				/* pad record with blanks so that it is at least 80 chars long */
				strcat(card, "                                                                                ");

				// If the key is found, then add it
				strcpy(subheadptr, card);
				subheadptr += 80;
				(*nsubkeys)++;
			}
		}

		*subheadptr = '\0';   /* terminate the header string */
		/* minimize the allocated memory */
		*subheader = (char *) realloc(*subheader, (*nsubkeys *80) + 1);


	}

	//TODO: distribue to other processes..

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
	char * header ;
	int nkeys, nreject;

	if (rank ==0){
		if (int status = wcshdo(WCSHDO_all, wcs, &nkeys, &header)) {
			printf("%4d: %s.\n", status, wcs_errmsg[status]);
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&nkeys,   1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nsubkeys,1,MPI_INT,0,MPI_COMM_WORLD);

	if (rank != 0 ){
		header    = (char *) calloc ( (nkeys    + 1) * 80 + 1, 1);
		subheader = (char *) calloc ( (nsubkeys + 1) * 80 + 1, 1);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(   header,(   nkeys+1)*80+1,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(subheader,(nsubkeys+1)*80+1,MPI_CHAR,0,MPI_COMM_WORLD);

	// Back into wcs structure for rank != 0
	if (rank != 0) {
		if ( int status = wcspih(header, nkeys, WCSHDR_all, 2, &nreject, &nwcs, &wcs) ) {
			fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
			return 1;
		}
	}

	free(header);

	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.

#endif

	return 0;

}

int write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,
		char dtype, void *data, string table_name, bool fits_already_exist, char * subheader, int nsubkeys) {
	// all angles in degrees

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	int naxis = 2; // number of dimensions
	long naxes[] = { NAXIS1, NAXIS2 }; // size of dimensions
	long fpixel[] = { 1, 1 }; // index for write_pix
	long long ndata = NAXIS1 * NAXIS2; // number of data points

	char *header, *hptr;
	int nkeyrec;

	if (fits_already_exist) {
		if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
	} else {
		// create fits file
		if (fits_create_file(&fp, fname.c_str(), &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}

		// Dummy image
		long dummy_naxes[] = {0, 0};
		if (fits_create_img(fp, DOUBLE_IMG, 0, dummy_naxes, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}

		// Add Comments to the header
		std::vector<std::string> comments;
		comments.push_back(" ");
		comments.push_back("This fits file was generated by SANEPIC version "+ (string) PACKAGE_VERSION);
		comments.push_back("For more informations about SANEPIC and for SANEPIC updates, please");
		comments.push_back("check our dedicated website at http://www.ias.u-psud.fr/sanepic ");
		comments.push_back(" ");

		for (u_int ii = 0; ii < comments.size(); ii++) {
			if (fits_write_comment(fp, (char *) comments[ii].c_str(), &fits_status))
				return 1;
		}

		if (fits_write_chksum(fp, &fits_status)) {
			cout << "error checksum !\n";
			return 1;
		}

	}

	// create fits image (switch on data type)
	switch (dtype) {
	case 'd': // double
		if (fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	case 'l': // long
		if (fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	default:
		printf("write_fits: data type %c not supported. Exiting.\n", dtype);
		return 1;
	}


	// Transform wcsprm struture to header
	if ((fits_status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
		printf("wcshdo ERROR %d: %s.\n", fits_status,
				wcs_errmsg[fits_status]);
		return 1;
	}

	hptr = header;
	// write it to the fits file
	for (int keyrec = 0; keyrec < nkeyrec; keyrec++, hptr += 80) {
		if (fits_write_record(fp, (const char*) hptr, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
	}
	free(header);

	// ... add the subheader
	hptr = subheader;
	for (int keyrec = 0; keyrec < nsubkeys; keyrec++, hptr += 80) {
		if (fits_write_record(fp, (const char*) hptr, &fits_status)) {
			fits_report_error(stderr, fits_status); return 1;
		}
	}

	// write date to file
	if (fits_write_date(fp, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_update_key(fp, TSTRING, (char *) "EXTNAME",
			(void*) (table_name.c_str()), (char *) "table name", &fits_status))
		return 1;

	// write map data
	switch (dtype) {
	case 'd': // double
		if (fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data,
				&fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	case 'l': // long
		if (fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	}

	// Add Comments to the header
	std::vector<std::string> comments;
	comments.push_back(" ");
	comments.push_back("This fits file was generated by SANEPIC version "+ (string) PACKAGE_VERSION);
	comments.push_back("For more informations about SANEPIC and for SANEPIC updates, please");
	comments.push_back("check our dedicated website at http://www.ias.u-psud.fr/sanepic ");
	comments.push_back(" ");

	for (u_int ii = 0; ii < comments.size(); ii++) {
		if (fits_write_comment(fp, (char *) comments[ii].c_str(), &fits_status))
			return 1;
	}

	if (fits_write_chksum(fp, &fits_status)) {
		cout << "error checksum !\n";
		return 1;
	}

	// close file
	if (fits_close_file(fp, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;
}


int write_fits_inifile(std::string fname,
		struct param_common dir, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct param_sanePS PS_param,
		struct param_sanePic Pic_param, struct param_saneInv Inv_param) {

	fitsfile *fptr;
	int fits_status = 0;
	//	long naxes[2] = { 1, 1 };
	std::vector<string> key_vect;
	std::vector<string> val_vect;
	std::vector<string> com_vect;



	key_vect.push_back("[common]");   val_vect.push_back("");  com_vect.push_back("");
	export_param_common(dir,          key_vect, val_vect, com_vect);
	key_vect.push_back("[sanePos]");  val_vect.push_back("");  com_vect.push_back("");
	export_param_sanePos(pos_param,   key_vect, val_vect, com_vect);
	key_vect.push_back("[saneProc]");  val_vect.push_back("");  com_vect.push_back("");
	export_param_saneProc(proc_param, key_vect, val_vect, com_vect);
	key_vect.push_back("[saneInv]");  val_vect.push_back("");  com_vect.push_back("");
	export_param_saneInv(Inv_param,   key_vect, val_vect, com_vect);
	key_vect.push_back("[sanePic]");  val_vect.push_back("");  com_vect.push_back("");
	export_param_sanePic(Pic_param,   key_vect, val_vect, com_vect);
	key_vect.push_back("[sanePS]");   val_vect.push_back("");  com_vect.push_back("");
	export_param_sanePS(PS_param,     key_vect, val_vect, com_vect);


	if (fits_open_file(&fptr, fname.c_str(), READWRITE, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	char *ttype[] = { (char*) "key", (char*) "values", (char*) "comments" };
	char *tform[] = { tableFormat(key_vect), tableFormat(val_vect),	tableFormat(com_vect) };
	char *tunit[] = { (char*) "None", (char*) "None", (char*) "None" };

	char **key_char, **val_char, **com_char;

	key_char = vString2carray(key_vect);
	val_char = vString2carray(val_vect);
	com_char = vString2carray(com_vect);

	if (fits_create_tbl(fptr, BINARY_TBL, key_vect.size(), 3, ttype, tform,	tunit, (char*) "IniFile", &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 1, 1, 1, key_vect.size(), key_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 2, 1, 1, val_vect.size(), val_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 3, 1, 1, com_vect.size(), com_char, &fits_status))
		return 1;

	if (fits_write_chksum(fptr, &fits_status)) {
		cout << "error checksum !\n";
		return 1;
	}

	// close file
	if (fits_close_file(fptr, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	for (long ii = 0; ii < (long) key_vect.size(); ii++) {
		delete[] key_char[ii];
		delete[] val_char[ii];
		delete[] com_char[ii];
	}

	delete[] key_char;
	delete[] val_char;
	delete[] com_char;
	delete[] tform[0];
	delete[] tform[1];
	delete[] tform[2];

	return 0;
}

int write_fits_inputfile(std::string fname, struct samples samples_struct) {

	fitsfile *fptr;
	int fits_status = 0;

	long nsamples = samples_struct.ntotscan;

	std::vector<string> fits_svect;
	std::vector<string> bolo_svect;
	std::vector<string> index_svect;
	std::vector<string> fsamp_svect;
	std::vector<string> fhp_svect;


	std::vector<string> noise_svect;
	std::vector<string> fcut_svect;

	for (long ii=0; ii<nsamples; ii++){

		fits_svect.push_back(Basename(samples_struct.fitsvect[ii]));
		bolo_svect.push_back(Basename(samples_struct.bolovect[ii]));
		index_svect.push_back(StringOf(samples_struct.scans_index[ii]));
		fsamp_svect.push_back(StringOf(samples_struct.fsamp[ii]));
		fhp_svect.push_back(StringOf(samples_struct.fhp[ii]));


		noise_svect.push_back(Basename(samples_struct.noisevect[ii]));
		fcut_svect.push_back(StringOf(samples_struct.fcut[ii]));
	}


	if (fits_open_file(&fptr, fname.c_str(), READWRITE, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	char *ttype[] = { (char*) "DataFile", (char *) "fsamp", (char*) "fhp", (char*) "BoloFile", (char*) "ScanIndex", (char*) "NoiseFile", (char*) "fcut"};
	char *tform[] = { tableFormat(fits_svect), tableFormat(fsamp_svect), tableFormat(fhp_svect), tableFormat(bolo_svect),	tableFormat(index_svect), tableFormat(noise_svect), tableFormat(fcut_svect) };
	char *tunit[] = { (char*) "None", (char*) "Hz", (char *) "Hz", (char*) "None", (char*) "None" , (char*) "None", (char*) "Hz"};

	char **fits_char, **fsamp_char, **fhp_char, **bolo_char, **index_char, **noise_char, **fcut_char;

	fits_char  = vString2carray(fits_svect);
	fsamp_char = vString2carray(fsamp_svect);
	fhp_char   = vString2carray(fhp_svect);
	bolo_char  = vString2carray(bolo_svect);
	index_char = vString2carray(index_svect);
	noise_char = vString2carray(noise_svect);
	fcut_char  = vString2carray(fcut_svect);

	if (fits_create_tbl(fptr, BINARY_TBL, nsamples, 7, ttype, tform,	tunit, (char*) "InputFiles", &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 1, 1, 1, nsamples, fits_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 2, 1, 1, nsamples, fsamp_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 3, 1, 1, nsamples, fhp_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 4, 1, 1, nsamples, bolo_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 5, 1, 1, nsamples, index_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 6, 1, 1, nsamples, noise_char, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 7, 1, 1, nsamples, fcut_char, &fits_status))
		return 1;


	if (fits_write_chksum(fptr, &fits_status)) {
		cout << "error checksum !\n";
		return 1;
	}

	// close file
	if (fits_close_file(fptr, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	for (long ii = 0; ii < nsamples; ii++) {
		delete[] fits_char[ii];
		delete[] fsamp_char[ii];
		delete[] fhp_char[ii];
		delete[] bolo_char[ii];
		delete[] index_char[ii];
		delete[] noise_char[ii];
		delete[] fcut_char[ii];
	}

	delete[] fits_char;
	delete[] fsamp_char;
	delete[] fhp_char;
	delete[] bolo_char;
	delete[] index_char;
	delete[] noise_char;
	delete[] fcut_char;

	delete[] tform[0];
	delete[] tform[1];
	delete[] tform[2];
	delete[] tform[3];
	delete[] tform[4];
	delete[] tform[5];
	delete[] tform[6];

	return 0;
}


int copy_fits_mask(std::string fname, std::string maskfile) {

	fitsfile *fptr, *outfptr;
	int fits_status = 0;

	if (fits_open_file(&fptr, maskfile.c_str(), READONLY, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_open_file(&outfptr, fname.c_str(), READWRITE, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", 0, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_copy_hdu(fptr, outfptr, 0, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	// close files
	if (fits_close_file(fptr, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_close_file(outfptr, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

}
int read_wcs(string fname, string extname, struct wcsprm *& wcs)
/*
 * Read the wcs from the extension 'extname' from the 'fname' fits file,
 * Return a wcs structure
 */
{
	fitsfile *fptr;
	int status = 0,  wcsstatus[NWCSFIX];
	char *header;
	int nkeyrec, nwcs, nreject;

	// Open the fits file...
	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	// ... and move to the 'extname'
	if (extname.compare("") != 0) {
		if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) extname.c_str(), 0, &status)) {
			fits_report_error(stderr, status);
			return 1;
		}
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
						wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status = wcsset(wcs))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}

	if (fits_close_file(fptr, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	free(header);
	return (0);
}

int read_mask_wcs(string fname, string extname, struct wcsprm *& wcs,
		long &NAXIS1, long &NAXIS2, short *& data)
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
	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	// ... and move to the 'extname'
	if (extname.compare("") != 0) {
		if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) extname.c_str(), 0,
				&status)) {
			fits_report_error(stderr, status);
			return 1;
		}
	}

	// Retrieve the image size
	if (fits_get_img_size(fptr, 2, naxes, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	NAXIS1 = naxes[0];
	NAXIS2 = naxes[1];
	// Allocate the image container and read its depending on the type
	//	switch (dtype) {
	//	case 's':
	data = new short[NAXIS1 * NAXIS2];
	if (fits_read_pix(fptr, TSHORT, fpixel, (long long) NAXIS1 * NAXIS2, 0,
			data, &anynul, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
						wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status = wcsset(wcs))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}

	if (fits_close_file(fptr, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	free(header);
	return (0);
}

int read_fits_signal(string fname, double *S, long long* indpix, long NAXIS1,
		long NAXIS2, struct wcsprm * wcs)
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

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Image", 0, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs,
			&wcs_fits))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs_fits, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
						wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status = wcsset(wcs_fits))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}

	if (fits_get_img_size(fptr, 2, naxes, &status))
		return 1;

	imNAXIS1 = (long) naxes[0];
	imNAXIS2 = (long) naxes[1];

	if (compare_wcs(fname, wcs, wcs_fits, NAXIS1, NAXIS2, imNAXIS1, imNAXIS2))
		return 1;

	// Initialize the data container
	map = new double[imNAXIS1 * imNAXIS2];

	if (fits_read_pix(fptr, TDOUBLE, fpixel, (long long) imNAXIS1 * imNAXIS2,
			0, map, &anynul, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	// work correctly
	for (long ii = 0; ii < imNAXIS1; ii++) {
		for (long jj = 0; jj < imNAXIS2; jj++) {
			mi = jj * imNAXIS1 + ii;
			if (indpix[mi] >= 0) {
				S[indpix[mi]] = map[mi];
			}
		}
	}

	if (fits_close_file(fptr, &status)) {
		fits_report_error(stderr, status);
		return 1;
	}

	//	nwcs=1;
	wcsvfree(&nwcs, &wcs_fits);
	free(header);
	return 0;
}

int save_keyrec(string path, struct wcsprm * wcs, long NAXIS1, long NAXIS2, char * subheader, int nsubkeys) {

	FILE *fout;
	int nkeyrec, status;
	char *header, *hptr;
	string filename;

	if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		return 1;
	}

	filename = path + "mapHeader.keyrec";

	fout = fopen(filename.c_str(), "w");
	if (fout == NULL) {
		fputs("Creation error : File error on mapHeader.keyrec\n", stderr);
		return (1);
	}

	fprintf(fout, "NAXIS1  = %20ld / %-47s\n", NAXIS1, "length of data axis 1");
	fprintf(fout, "NAXIS2  = %20ld / %-47s\n", NAXIS2, "length of data axis 2");

	hptr = header;
	for (int ii = 0; ii < nkeyrec; ii++, hptr += 80) {
		fprintf(fout, "%.80s\n", hptr);
	}
	fclose(fout);
	free(header);

	filename = path + "subHeader.keyrec";

	fout = fopen(filename.c_str(), "w");
	if (fout == NULL) {
		fputs("Creation error : File error on subHeader.keyrec\n", stderr);
		return (1);
	}
	hptr = subheader;
	for (int ii = 0; ii < nsubkeys; ii++, hptr += 80) {
		fprintf(fout, "%.80s\n", hptr);
	}
	fclose(fout);


	return 0;
}

int print_MapHeader(struct wcsprm *wcs) {

	int nkeyrec;
	char * header, *hptr;
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

int read_keyrec(string tmpdir, struct wcsprm * & wcs, long * NAXIS1,
		long * NAXIS2, char ** subheader, int *nsubkeys, int rank) {

	char *memblock = NULL;
	int nkeyrec = 0, nreject, nwcs, status;
	string filename;
	char * subheadptr;

	if (rank == 0) {

		FILE *fin;
		size_t result;
		int size;

		filename = tmpdir + "mapHeader.keyrec";

		fin = fopen(filename.c_str(), "r");
		if (fin == NULL) {
			fputs("Read error : File error on mapHeader.keyrec", stderr);
			return 1;
		}

		fseek(fin, 0L, SEEK_END); /* Position to end of file */
		size = ftell(fin); /* Get file length */
		rewind(fin); /* Back to start of file */

		nkeyrec = size / 81;

		char comment[47];

		// Read the two first lines, NAXIS1/NAXIS2
		result = fscanf(fin, "NAXIS1  = %20ld / %47c\n", NAXIS1, (char *) &comment);
		if (result != 2 )
			cerr << "EE - Read error on mapHeader.keyrec "  << endl;
		result = fscanf(fin, "NAXIS2  = %20ld / %47c\n", NAXIS2, (char *) &comment);
		if (result != 2 )
			cerr << "EE - Read error on mapHeader.keyrec " << endl;

		memblock = new char[(nkeyrec - 2) * 80+1];

		for (int ii = 0; ii < nkeyrec; ii++) {
			result = fread(&memblock[ii * 80], 80, sizeof(char), fin);
			fseek(fin, 1, SEEK_CUR); // skip newline char
		}
		fclose(fin);

		// Read the subHeader file...
		filename = tmpdir + "subHeader.keyrec";

		fin = fopen(filename.c_str(), "r");
		if (fin == NULL) {
			fputs("WW - Read error : File error on subHeader.keyrec", stderr);
			*nsubkeys = 0 ;
		} else {

			fseek(fin, 0L, SEEK_END); /* Position to end of file */
			size = ftell(fin); /* Get file length */
			rewind(fin); /* Back to start of file */

			*nsubkeys = size / 81;

			if ((*nsubkeys) > 0 ){
				*subheader = new char[(*nsubkeys) * 80+1];
				subheadptr = *subheader;

				for (int ii = 0; ii < (*nsubkeys); ii++) {

					result = fread(&subheadptr[ii * 80], 80, sizeof(char), fin);
					fseek(fin, 1, SEEK_CUR); // skip newline char
				}
			} else {
				fputs("WW - Read error : no keys in subHeader.keyrec", stderr);
			}
			fclose(fin);
		}

	}

#ifdef USE_MPI
	//	int position=0;
	MPI_Bcast(NAXIS1,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(NAXIS2,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&nkeyrec,1,MPI_INT,0,MPI_COMM_WORLD);

	if(rank!=0)
		memblock = new char [(nkeyrec-2)*80];

	MPI_Bcast(memblock, (nkeyrec-2)*80, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

	/* Parse the primary header of the FITS file. */
	/* -2 to handle the first two NAXIS? keyword */
	if ((status = wcspih(memblock, nkeyrec - 2, WCSHDR_all, 2, &nreject, &nwcs,
			&wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
		return 1;
	}
	delete[] memblock;
	//	free(comment);

	return 0;
}

int compare_wcs(std::string fname, struct wcsprm *wcs, struct wcsprm *wcs_fits,
		long NAXIS1, long NAXIS2, long imNAXIS1, long imNAXIS2) {

	// compatibility verifications !

	if (wcs_fits->crpix[0] != wcs->crpix[0] || wcs_fits->crpix[1]
	                                                           != wcs->crpix[1]) {
		cout
		<< "CRPIX are different between mapheader.keyrec and the map_file : "
		<< fname << endl;
		cout
		<< "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if (wcs_fits->crval[0] != wcs->crval[0] || wcs_fits->crval[1]
	                                                           != wcs->crval[1]) {
		cout
		<< "CRVAL are different between mapheader.keyrec and the map_file : "
		<< fname << endl;
		cout
		<< "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if (wcs_fits->cdelt[0] != wcs->cdelt[0] || wcs_fits->cdelt[1]
	                                                           != wcs->cdelt[1]) {
		cout
		<< "CDELT are different between mapheader.keyrec and the map_file : "
		<< fname << endl;
		cout
		<< "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if (wcs_fits->lonpole != wcs->lonpole || wcs_fits->latpole != wcs->latpole) {
		cout
		<< "LONPOLE and/or LATPOLE are different between mapheader.keyrec and the map_file : "
		<< fname << endl;
		cout
		<< "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if (NAXIS1 != imNAXIS1 || NAXIS2 != imNAXIS2) {
		cout
		<< "NAXIS are different between mapheader.keyrec and the map_file : "
		<< fname << endl;
		cout
		<< "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}

	return 0;

}

