
#include "map_mak_tools.h"



int test_map_mak(struct samples samples_struct, long ns,struct param_sanePre proc_param,
		std::string outdir,	std::vector<std::string> det,long ndet, double f_lppix, long iframe, int para_bolo_indice, int para_bolo_size){



	double *data, *data_lp;
	int *flag=NULL;
	long long *samptopix;

	string field1, fits_filename, fname;

	data_lp = new double[ns];
	samptopix = new long long[ns];

	fill(data_lp,data_lp+ns,0.0);
	fill(samptopix,samptopix+ns,0);

	fits_filename = samples_struct.fitsvect[iframe];
	copy_fits(fits_filename, outdir, ns);

	//	getchar();

	for (long idet1=para_bolo_indice*ndet/para_bolo_size;idet1<(para_bolo_indice+1)*ndet/para_bolo_size;idet1++){

		field1 = det[idet1];

		fname = fits_filename + "_" + field1 + ".fits";

		fill(data_lp,data_lp+ns,0.0);

		long test_ns;
		if(read_signal_from_fits(fits_filename, field1, data, test_ns))
			return 1;
		if (test_ns != ns) {
			cerr << "Read signal does not correspond to frame size : Check !!" << endl;
			return 1;
		}

		if(read_flag_from_fits(fits_filename , field1, flag, test_ns))
			return 1;
		if (test_ns != ns) {
			cerr << "Read flag does not correspond to frame size : Check !!" << endl;
			return 1;
		}


		MapMakePreProcessData(data,  flag, ns, proc_param, f_lppix, data_lp, NULL);

		//		write_to_fits !
		write_to_fits_data_lp(fits_filename, data_lp, outdir, field1, ns);

		delete [] flag;
		delete [] data;
	} // idet1


	delete[] data_lp;
	delete[] samptopix;

	return 0;
}

int copy_fits(std::string fits_name, std::string outdir, long ns_total){

	string fname2 = "!" + outdir + FitsBasename(fits_name) + "_test_map_mak.fits"; // output fits filename

	int status=0; // fits error status
	std::vector<string> det;
	long ndet;

	int *mask;
	// fits files pointer
	fitsfile * fptr;
	fitsfile *outfptr;


	// open original fits file
	if (fits_open_file(&fptr, fits_name.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// create ouput fixed fits file
	if (fits_create_file(&outfptr, fname2.c_str(), &status))
		fits_report_error(stderr, status);

	// read channels list
	read_bolo_list(fits_name, det, ndet);

	// Copy primary Header
	fits_copy_header(fptr, outfptr, &status);

	// channels
	copy_channels(fptr, outfptr);


	// copy mask
	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status); // move input pointer to mask
	fits_copy_header(fptr, outfptr, &status); // copy header to ouput

	for(long jj=0;jj<ndet;jj++){ // for each detector (column)
		read_flag_from_fits(fits_name, det[jj], mask, ns_total); // read input mask row
		insert_array_in_image(fptr, outfptr, det[jj], mask, ns_total); // insert the filled mask row in ouput table
		delete [] mask;
	}

	// initiate signal !
	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status); // move input pointer to signal
	fits_copy_header(fptr, outfptr, &status); // copy header to ouput

	fits_write_null_img(outfptr, 1, ns_total*ndet, &status);


	// close both fits files
	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	if (fits_close_file(outfptr, &status))
		fits_report_error(stderr, status);

	return 0;
}


int write_to_fits_data_lp(std::string fits_name, double *data_lp, string outdir,string field, long ns_total){

	string fname2 = outdir + FitsBasename(fits_name) + "_test_map_mak.fits"; // output fits filename

	int status=0; // fits error status

	// fits files pointer
	fitsfile * fptr;
	fitsfile *outfptr;


	// open original fits file
	if (fits_open_file(&fptr, fits_name.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// create ouput fixed fits file
	if (fits_open_file(&outfptr, fname2.c_str(),READWRITE, &status))
		fits_report_error(stderr, status);

	fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "signal", NULL, &status); // move input pointer to signal
	insert_array_in_image(fptr, outfptr, field, data_lp, ns_total);

	// close both fits files
	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	if (fits_close_file(outfptr, &status))
		fits_report_error(stderr, status);

	return 0;
}
