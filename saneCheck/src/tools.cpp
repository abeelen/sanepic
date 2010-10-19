

#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <unistd.h>   // For access()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}



using namespace std;

void check_detector_is_in_fits(struct detectors det,struct detectors bolo_fits, string filename)
/*! this function determines whether the user list of detectors is correct or not */
{

	int mycount=0; /* used to count the number of wrong detectors name */

	for(int jj=0;jj< det.ndet; jj++){
		mycount = (int) count (bolo_fits.boloname.begin(), bolo_fits.boloname.end(), det.boloname[jj]); /* count the number of time a bolometer appears in the detector list */
		if(mycount==0) // if this detector never appears
			// Warn user
			cout << "Warning ! The detector " << det.boloname[jj] << " is not referenced in the fits " << filename << endl;
	}

}

int check_positionHDU(string fname,long ns,struct detectors det, int format, struct checkHDU &check_it)
/*! this function determines whether reference and offsets HDUs are present and have the correct size */
{

	fitsfile *fptr; // fits file pointer
	int status = 0; // fits error status
	int colnum; // fits binary table column index
	long ns_test=0; // used only to read tables
	long ndet_test=0; // used to read detector list size


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open fits file
		fits_report_error(stderr, status);

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status)){ // go to reference table
		fits_report_error(stderr, status); // reference position table was not found ?
		check_it.checkREFERENCEPOSITION=0;
		cout << "\"reference position\" was not found, or his Type should be Binary table" << endl;
	}else{

		fits_get_num_rows(fptr, &ns_test, &status);
		if(ns!=ns_test){ // check reference position table's size
			cout << "\"reference position\" has a wrong number of rows (must be equal to ns : " << ns << " )" << endl;
			return -1;
		}

		fits_get_num_cols(fptr, &colnum, &status);
		if(colnum!=3){ // check reference position table's size
			cout << "\"reference position\" has a wrong number of cols (must be equal to 3 : RA, DEC, PHI )" << endl;
			return -1;
		}else{


			fits_get_colnum(fptr, CASEINSEN, (char*) "RA", &colnum, &status);
			if(colnum!=1){ // check RA table's size
				cout << "\"RA\" was not found in \"reference position\"" << endl;
				return -1;
			}

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) "DEC", &colnum, &status);
			if(colnum!=2){ // check DEC table's size
				cout << "\"DEC\" table was not found in \"reference position\"" << endl;
				return -1;
			}

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) "PHI", &colnum, &status);
			if(colnum!=3){ // check PHI table's size
				cout << "\"PHI\" table was not found in \"reference position\"" << endl;
				return -1;
			}
		}
	}

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status)){ // go to offsets table
		fits_report_error(stderr, status); // offsets table was not found ?
		check_it.checkOFFSETS=0;
		cout << "\"offsets\" was not found, or his Type should be Binary table" << endl;
	}else{

		fits_get_num_rows(fptr, &ndet_test, &status);
		if(det.ndet!=ndet_test){ // check offsets table's size
			cout << "\"offsets\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
			return -1;
		}

		colnum=0;
		fits_get_num_cols(fptr, &colnum, &status);
		if(colnum!=3){ // check offsets table's size
			cout << "\"offsets\" has a wrong number of cols (must be equal to 3 : NAME, DX, DY )" << endl;
			return -1;
		}else{
			string name_table,dx_table,dy_table;

			if(format==1){ // check presence of dx and dy tables in offsets table : name depends on format
				name_table="names";
				dx_table="dX";
				dy_table="dY";
			}else{
				name_table="NAME";
				dx_table="DX";
				dy_table="DY";
			}


			colnum=0; // check offsets table's size
			fits_get_colnum(fptr, CASEINSEN, (char*) name_table.c_str(), &colnum, &status);
			if(colnum!=1){
				cout << "\"NAME\" table was not found in \"offsets\"" << endl;
				return -1;
			}

			colnum=0; // check dx table's presence
			fits_get_colnum(fptr, CASEINSEN, (char*) dx_table.c_str(), &colnum, &status);
			if(colnum!=2){
				cout << "\"DX\" table was not found in \"offsets\"" << endl;
				return -1;
			}

			colnum=0; // check dy table's presence
			fits_get_colnum(fptr, CASEINSEN, (char*) dy_table.c_str(), &colnum, &status);
			if(colnum!=3){
				cout << "\"DY\" table was not found in \"offsets\"" << endl;
				return -1;
			}
		}
	}

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	return 0;

}

int check_commonHDU(string fname,long ns,struct detectors det, struct checkHDU &check_it)
/*! this function determines whether channels, time, signal and mask HDUs are present and have the correct size */
{

	fitsfile *fptr;
	int status = 0;
	int colnum;
	long ndet_test=0;
	int naxis=0;
	long naxes[2] = { 1, 1 };


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return -1;
	}

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status)){ // go to channels table
		fits_report_error(stderr, status); // is the table present ?
		cout << "\"channels\" was not found, or his Type should be Binary table" << endl;
		return -1;
	}else{

		ndet_test=0;
		fits_get_num_rows(fptr, &ndet_test, &status);
		if(ndet_test!=det.ndet){ // Is the size correct ?
			cout << "\"channels\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
			return -1;
		}

		colnum=0;
		fits_get_num_cols(fptr, &colnum, &status);
		if(colnum!=1){  // Is the size correct ?
			cout << "\"channels\" has a wrong number of cols (must be equal to 1 : NAMES )" << endl;
			return -1;
		}else{

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) "NAMES", &colnum, &status);
			if(colnum!=1){ // Table NAMES is present ?
				fits_get_colnum(fptr, CASEINSEN, (char*) "NAME", &colnum, &status);
				if(colnum!=1){ // Table NAME is present ?
					cout << "\"NAMES\" table was not found in \"channels\"" << endl;
					return -1;
				}
			}
		}
	}


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status)){
		fits_report_error(stderr, status); // time table is present ?
		cout << "\"time\" was not found, or his Type should be image" << endl;
		return -1;
	}else{
		if(fits_get_img_dim(fptr, &naxis, &status)){ // check size
			fits_report_error(stderr, status);
			return -1;
		}else{

			if(naxis!=1){
				cout << "\"time\" has a wrong number of columns, should be equal to 1 " << endl;
				return -1;
			}else{
				if (fits_get_img_size(fptr, 2, naxes, &status)){
					fits_report_error(stderr, status);
					return -1;
				}
				if(naxes[0]!=ns){
					cout << "\"time\" has a wrong number of elements, should be equal to ns : " << ns << endl;
					return -1;
				}
			}
		}
	}

	status = 0;
	naxes[0]=1;
	naxes[1]=1;


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status)){
		fits_report_error(stderr, status); // signal image is present ?dr
		cout << "\"signal\" was not found, or his Type should be image" << endl;
		return -1;
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status)){ // check size
			fits_report_error(stderr, status);
			return -1;
		}
		if(naxis != 2){
			fits_report_error(stderr,BAD_NAXIS);
			cout << "\"signal\" must have 2 dimensions" << endl;
			return -1;
		}else{
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return -1;
			}else{

				if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
					cout << "\"signal\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
					return -1;
				}
			}
		}
	}

	status = 0;
	naxes[0]=1;
	naxes[1]=1;


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status)){
		fits_report_error(stderr, status); // mask image is present ?
		cout << "\"mask\" was not found, or his Type should be image" << endl;
		return -1;
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status)){ // check size
			fits_report_error(stderr, status);
			return -1;
		}else{
			if(naxis != 2){
				fits_report_error(stderr,BAD_NAXIS);
				cout << "\"mask\" must have 2 dimensions" << endl;
				return -1;
			}else{

				if (fits_get_img_size(fptr, 2, naxes, &status)){
					fits_report_error(stderr, status);
					cout << "\"mask\" must have 2 dimensions" << endl;
					return -1;
				}

				if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
					cout << "\"mask\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
					return -1;
				}
			}
		}
	}

	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return -1;
	}


	return 0;
}


int check_altpositionHDU(string fname,long ns,struct detectors det, struct checkHDU &check_it)
/*! check RA/DEc table presence : only for HIPE format */
{

	fitsfile *fptr; // fits pointer
	int status = 0; // fits error status
	int naxis=0;
	long naxes[2] = { 1, 1 };

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open fits
		fits_report_error(stderr, status);

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status)){ // move to ra table
		fits_report_error(stderr, status); // Does it exists ?
		check_it.checkRA=0;
		cout << "\"ra\" was not found, or his Type should be image" << endl;
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status)){ // get size
			fits_report_error(stderr, status);
			return -1;
		}else{
			if(naxis != 2){
				fits_report_error(stderr,BAD_NAXIS);
				cout << "\"ra\" must have 2 dimensions" << endl;
				return -1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return -1;
			}else{

				if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){ // check size
					cout << "\"ra\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
					return -1;
				}
			}
		}
	}

	status = 0;


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status)){
		fits_report_error(stderr, status); // move to DEC table
		check_it.checkDEC=0;
		cout << "\"dec\" was not found, or his Type should be image" << endl;
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status)){ // get table size
			fits_report_error(stderr, status);
			return -1;
		}else{
			if(naxis != 2){
				fits_report_error(stderr,BAD_NAXIS);
				cout << "\"dec\" must have 2 dimensions" << endl;
				return -1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return -1;
			}else{

				if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){ // check size
					cout << "\"dec\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
					return -1;
				}
			}
		}
	}

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	return 0;
}

int check_NAN_positionHDU(string fname,long ns,struct detectors det, struct checkHDU check_it)
/*! Check presence of non-flagged NANs in position tables */
{

	long ns_test=0;
	double *ra; // ra, dec, phi and offsets table are used to read the fits tables
	double *dec,*phi;
	double **offsets;

	int *flag; // to read mask table
	if(check_it.checkREFERENCEPOSITION){
		for(int ii=0;ii<det.ndet;ii++){
			if(read_ReferencePosition_from_fits(fname, ra, dec, phi, ns_test)) // read ra dec and phi
				return 1;

			if(read_flag_from_fits(fname, det.boloname[ii], flag, ns)) // read mask image
				return 1;

			for(long jj=0;jj<ns_test;jj++){ // check NANs
				if(isnan(ra[jj])&&(flag[jj]==0)){
					cout << "Warning <! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				}
				if(isnan(dec[jj])&&(flag[jj]==0)){
					cout << "Warning <! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				}
				if(isnan(phi[jj])&&(flag[jj]==0)){
					cout << "Warning <! a NAN has been found in \"phi\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				}
			}

			delete [] ra;
			delete [] dec;
			delete [] phi;
			delete [] flag;
		}
	}

	if(check_it.checkOFFSETS){
		// flag indepedent
		if(read_all_bolo_offsets_from_fits(fname, det.boloname, offsets)) // read offsets table
			return 1;

		for(int ii=0;ii<det.ndet;ii++)
			if(isnan(offsets[ii][0])||isnan(offsets[ii][1])){ // check NANs
				cout << "Warning ! a NAN has been found in \"offsets\" table for bolometer n° " << ii << endl;
				cout << "You should not take this detector for the computation of Sanepic\n";
			}

		// clean up
		free_dmatrix(offsets,(long)0,det.ndet-1,(long)0,2-1);
	}

	return 0;

}

int check_NAN_commonHDU(string fname,long ns,struct detectors det, struct checkHDU check_it)
/*! check presence of non-flagged NANs in time, signal and mask tables */
{

	long ns_test=0;
	double *signal;
	int *flag;
	double *time;


	// check nans in mask image
	for(long jj=0;jj<det.ndet;jj++){
		if(read_flag_from_fits(fname, det.boloname[jj], flag, ns))
			return 1;
		for(int kk=0;kk<ns;kk++){
			if(isnan(flag[kk])){
				cout << "Warning <! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
			}
		}
		delete [] flag;
	}

	// check nans in time image
	if(read_time_from_fits(fname, time, ns))
		return 1;
	for(long jj=0;jj<ns;jj++){
		if(isnan(time[jj])&&(flag[jj]==0)){
			cout << "Warning ! a NAN has been found in \"time\" table for sample n° " << jj << endl;
		}
	}

	delete [] time;


	// check nans in signal
	for(int ii=0;ii<det.ndet;ii++){
		if(read_signal_from_fits(fname, det.boloname[ii], signal,ns_test))
			return 1;
		for(long jj=0;jj<ns_test;jj++){
			if(isnan(signal[jj])&&(flag[jj]==0)){
				cout << "Warning <! a NAN has been found in \"signal\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
		}
		delete [] signal;
	}

	return 0;
}

int check_NAN_altpositionHDU(string fname,long ns,struct detectors det, struct checkHDU check_it)
/*! check non-flagged NANs in RA/DEC HIPE format */
{


	long ns_test=0;
	double *ra;
	double *dec;
	int *flag;


	for(int ii=0;ii<det.ndet;ii++){
		if(read_flag_from_fits(fname, det.boloname[ii], flag, ns_test)) // read mask image
			return 1;

		if(read_ra_dec_from_fits(fname, det.boloname[ii], ra, dec, ns_test)) // read RA and DEC tables
			return 1;

		for(long jj=0;jj<ns_test;jj++){
			if(isnan(ra[jj])&&(flag[jj]==0)){ // check for NANs in RA
				cout << "Warning <! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
			if(isnan(dec[jj])&&(flag[jj]==0)){ // check for NANs in DEC
				cout << "Warning <! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
		}
		delete [] ra;
		delete [] dec;
		delete [] flag;

	}

	return 0;
}

bool check_bolos(std::vector<string> bolo_fits_vect, std::vector<string> bolo_fits_0_vect)
/*! Check that a given fits file bolometer list is exactly the same as the first input fits file one */
{


	if((long)bolo_fits_vect.size()!=(long)bolo_fits_0_vect.size()){
		cout << "The channel lists are not the same !\n";
		return 1;
	}

	struct sortclass_string sortobject;
	sort(bolo_fits_vect.begin(), bolo_fits_vect.end(), sortobject);
	sort(bolo_fits_0_vect.begin(), bolo_fits_0_vect.end(), sortobject);

	for(long jj=0; jj< (long)bolo_fits_vect.size(); jj++)
		if(bolo_fits_vect[jj]!=bolo_fits_0_vect[jj]){
			cout << "The channel lists are not the same !\n";
			return 1;
		}

	return 0;
}

int check_flag(string fname,struct detectors det,long ns, string outname,long *&bolos_global,long *&bolos_global_80,double *percent_tab, struct checkHDU check_it)
/*!  Lookfor fully or more than 80% flagged detectors, also flag singletons */
{

	int *flag; // mask table
	long sum=0; // used to compute the percentage of flagged data

	//	int marge = 19; // a singleton is defined by : ['marge' samples flagged] 0 ['marge' samples flagged]
	long ii=1;// singleton seeker indexes
	//	long tt=0;// singleton seeker indexes
	//	long rr=0;// singleton seeker indexes
	cout << "ns : " << ns << endl;

	for(int jj=0;jj<det.ndet;jj++){

		//		cout << det.boloname[jj];
		ii=1;
		sum=0;
		if(read_flag_from_fits(fname, det.boloname[jj], flag, ns))
			return 1;

		while(ii<ns-1){

			if(flag[ii]!=0)
				sum++;
			ii++;

		}

		if(sum==ns){ // fully flagged detector found
			cout << "Warning ! " << det.boloname[jj] << " is totally flagged" << endl;
			bolos_global[jj]=1;
		}else{
			if(sum>80*ns/100){ // valid worst detector found
				double percent = sum/(double)ns*100;
				cout << "Warning ! " << det.boloname[jj] << " is more than 80% flagged : " << percent << endl;
				bolos_global_80[jj]=1;
				percent_tab[jj]=percent;
			}else if(sum>50*ns/100){
				double percent = sum/(double)ns*100;
				cout << "Warning ! " << det.boloname[jj] << " is more than 50% flagged : " << percent << endl;
			}
		}

		delete [] flag;
	}

	return 0;

}


int check_time_gaps(string fname,long ns, double fsamp, struct common dir, struct checkHDU check_it)
/*! check for time gaps in time table */
{


	std::vector<long> indice;

	std::ostringstream oss;
	oss << fname;
	string filename = oss.str();
	string fname2 = dir.tmp_dir + Basename(filename) + "_saneFix_indices.bin"; // output saneFix logfile filename

	double *time,*diff;
	double sum=0.0;
	std::vector<double> freq; // The whole frequency that can be extracted from the time vector

	if(read_time_from_fits(fname, time, ns))
		return 1;

	diff = new double [ns-1]; // differential time vector

	for(long jj=0;jj<ns-1;jj++){
		diff[jj]=(time[jj+1]-time[jj]);
		sum+=diff[jj]; // sum up time gaps
		freq.push_back(1/diff[jj]); // store frequency
	}

	struct sortclass_double sortobject;
	sort(freq.begin(), freq.end(), sortobject); // sort frequency vector

	std::vector<double>::iterator it;
	it = unique(freq.begin(),freq.end()); // get only unique values of frequencies




	long size_tmp = it - freq.begin(); // last unique frequency index

	long *counter;
	counter = new long [size_tmp];

	for (long tt=0;tt<size_tmp;tt++){
		counter[tt]=count(freq.begin(), freq.end(),freq[tt]); // Frequencies occurence
	}

	long maxi = *max_element(counter,counter+size_tmp); // max occurence

	long ind=0;

	for(long tt=0;tt<size_tmp;tt++)
		if(counter[tt]==maxi)
			ind=tt;

	double Populated_freq = freq[ind]; // the most populated frequency only is kept

	double zero_cinq_pourcent=0.5*Populated_freq/100; // compare to the user frequency given in ini file
	if(abs(Populated_freq-fsamp)/fsamp>zero_cinq_pourcent){
		cout << "Warning, the sampling frequency you have mentioned in the ini file seems to be wrong : \n";
		cout << "ini file : " << fsamp << " != " << Populated_freq << endl;
	}


	for(long jj=0;jj<ns-1;jj++){ // locate time gaps
		if((abs(diff[jj])>1.9/Populated_freq)){ //||(abs(diff[jj])<1/Populated_freq/1.95)
			cout << "WARNING ! Time gap at " << jj << " (" << time[jj] <<") : " << fixed <<  setprecision(8) << diff[jj] << endl;
			indice.push_back(jj); // store sample indice : where the time gap is
		}
	}

	sum=0.0;
	// recompute real frequency
	for(long jj=0;jj<ns-1;jj++){
		if(((int) count (indice.begin(), indice.end(), jj))==0)
			sum+=diff[jj];
	}
	double real_freq= (double)(ns-1-(long)indice.size())/sum; // dont take the gaps into account

	// print to std
	cout << "ini file fsamp : " << fsamp << " Most Populated freq : " << Populated_freq << " Recomputed real freq : " << real_freq << endl;

	std::ofstream file; // saneFix indice log file
	file.open(fname2.c_str(), ios::out | ios::trunc);
	if(!file.is_open()){
		cerr << "File [" << fname2 << "] Invalid." << endl;
		return 1;
	}else{
		// store real_freq for saneFix
		file << Populated_freq << " ";

		for(long ii = 0; ii<(long)indice.size();ii++)// store indices for saneFix
			file << indice[ii] << " ";
	}

	file.close();

	cout << endl;

	// clean up
	delete [] time;
	delete [] diff;

	return 0;

}

int check_bolo_gain(string fname,long ns, string bolo_gain_filename, struct detectors det, struct checkHDU check_it){


	long ns_test=0;
	double *signal_tot, *signal_samp;
	double *signal;
	signal_tot=new double[ns];
	signal_samp=new double[det.ndet];
	//	cout << det.ndet << endl;

	// test median !
	//	double *test;
	//	test=new double[4];
	//	test[0]=1;
	//	test[1]=2;
	//	test[2]=10;
	//	test[3]=4;
	//
	//	std::vector<double> test_vec(test, test+4);
	//	double med=median(test_vec);
	//	cout << test_vec[0] << " " << test_vec[1] << " " << test_vec[2] << " " << test_vec[3] <<  endl;
	//	cout << med << endl;
	//
	//	delete [] test;
	//	getchar();
	//---------------
	fill(signal_tot,signal_tot+ns,0.0);

	// sum up all detectors signal
	for(int jj=0;jj<ns;jj++){
		fill(signal_samp,signal_samp+det.ndet,0.0);

		if(read_sample_signal_from_fits(fname, jj+1, signal_samp, det.ndet))
			return 1;

		std::vector<double> signal_vec(signal_samp, signal_samp+det.ndet);
		//		cout << signal_vec[0] << endl;
		signal_tot[jj]=median(signal_vec);
		//		cout << signal_tot[jj] << endl;
	}

	fill(signal_samp,signal_samp+det.ndet,0.0);

	for(long jj=0;jj<10;jj++)
		cout << signal_tot[jj] << " " ;
	getchar();
	for(long idet=0;idet<det.ndet;idet++){
		std::vector<double> signal_vec;
		read_signal_from_fits(fname, det.boloname[idet], signal, ns_test);

		for(long ii=0;ii<ns;ii++)
			signal_vec.push_back(signal_tot[ii]/signal[ii]);

		signal_samp[idet]=median(signal_vec); // store gain in signal_samp

		//		signal_vec.clear();
		delete [] signal;
	}

	cout << "gain \n";
	for(long ii=0;ii<det.ndet;ii++)
		cout << signal_samp[ii] << " ";
	cout << endl;

	//clean up
	delete [] signal_samp;
	delete [] signal_tot;

	return 0;
}

double median(std::vector<double> vec){

	long size;
	struct sortclass_double sortobject;
	sort(vec.begin(), vec.end(), sortobject);
	size= (long)vec.size();

	long mid = size/ 2;

	return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2.0 : vec[mid];
}


void log_gen(long  *bolo_, string outname, struct detectors det, double *percent_tab)
/*! generating log files for user information */
{


	FILE *fp;
	long tot=0;

	fp=fopen(outname.c_str(),"w");

	for(long ii=0; ii< det.ndet; ii++)
		if(bolo_[ii]>0){
			string temp = det.boloname[ii]; // copy the name of the wrong or bad detector in this ascii log file
			if(percent_tab!=NULL)
				fprintf(fp, "%s %lf%%\n", (char*)temp.c_str(),percent_tab[ii]);
			else
				fprintf(fp, "%s\n", (char*)temp.c_str());
			tot++;
		}

	if(tot>0)
		cout << "There are " << tot << " bolometers in this file !\n";

	fclose(fp);

}

int read_sample_signal_from_fits(string filename, int sample, double *& signal_samp, long ndet){
	//TODO : Handle unit to transform to a common internal known unit

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int naxis = 0, anynul;
	long ns;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_img_dim(fptr, &naxis, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if(naxis != 2){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_get_img_size(fptr, 2, naxes, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Allocate Memory
	ns = naxes[0];
	double temp;

	for(int idet=1;idet<=(int)ndet;idet++){
		//		cout << idet << endl;
		// ---------------------------------------------
		// Retrieve the corresponding pix
		fpixel[0] = sample;
		fpixel[1] = idet;
		if (fits_read_pix(fptr, TDOUBLE, fpixel, 1, 0, &temp, &anynul, &status)){
			fits_report_error(stderr, status);
			return 1;
		}
		signal_samp[idet]=temp;
	}

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;

}
