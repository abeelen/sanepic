#include "InputFileIO.h"
#include "MPIConfiguration.h"
#include "DataIO.h"
#include "ParserFunctions.h"
#include "SaneCheckTools.h"

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


#include "gsl/gsl_math.h"
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort.h>


extern "C" {
//#include "nrutil.h"
#include <wcslib/prj.h>
#include <wcslib/cel.h>
#include <wcslib/sph.h>
#include <fitsio.h>
}



using namespace std;

void check_detector_is_in_fits(std::vector<std::string> det, long ndet, std::vector<std::string> bolo_fits, string filename)
/* this function determines whether the user list of detectors is correct or not */
{

	int mycount=0; /* used to count the number of wrong detectors name */

	for(int jj=0;jj< ndet; jj++){
		mycount = (int) count (bolo_fits.begin(), bolo_fits.end(), det[jj]); /* count the number of time a bolometer appears in the detector list */
		if(mycount==0) // if this detector never appears
			// Warn user
			cout << "Warning ! The detector " << det[jj] << " is not referenced in the fits " << filename << endl;
	}

}

long check_NAN_positionHDU(string fname,long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it)
/* Check presence of non-flagged NANs in position tables */
{

	long ns_test=0;
	double *lon, *lat, *phi;
	//	double **offsets = NULL;
	long nan_found=0;

	int *flag; // to read mask table
	if(check_it.checkREFERENCEPOSITION){
		for(int ii=0;ii<ndet;ii++){
			if(read_ReferencePosition_from_fits(fname, lon, lat, phi, ns_test))
				return 1;

			if(read_flag_from_fits(fname, det[ii], flag, ns)) // read mask image
				return 1;

			for(long jj=0;jj<ns_test;jj++){ // check NANs
				if(isnan(lon[jj])&&(flag[jj]==0)){
					cout << "Warning <! a NAN has been found in \"lon\" table for bolometer n° " << ii << " sample n° " << jj << endl;
					nan_found++;
				}
				if(isnan(lat[jj])&&(flag[jj]==0)){
					cout << "Warning <! a NAN has been found in \"lat\" table for bolometer n° " << ii << " sample n° " << jj << endl;
					nan_found++;
				}
				if(isnan(phi[jj])&&(flag[jj]==0)){
					cout << "Warning <! a NAN has been found in \"phi\" table for bolometer n° " << ii << " sample n° " << jj << endl;
					nan_found++;
				}
			}

			delete [] lon;
			delete [] lat;
			delete [] phi;
			delete [] flag;
		}
	}

	if(check_it.checkOFFSETS){
		//TODO: This fails on read.... Check...
		// flag indepedent
		//		if(read_all_bolo_offsets_from_fits(fname, det, offsets)) // read offsets table
		//			return 1;
		//
		//		for(int ii=0;ii<ndet;ii++)
		//			if(isnan(offsets[ii][0])||isnan(offsets[ii][1])){ // check NANs
		//				cout << "Warning ! a NAN has been found in \"offsets\" table for bolometer n° " << ii << endl;
		//				cout << "You should not take this detector for the computation of Sanepic\n";
		//				nan_found++;
		//			}

		// clean up
		//		free_dmatrix(offsets,(long)0,ndet-1,(long)0,2-1);
	}
	return nan_found;

}

long check_NAN_commonHDU(string fname,long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it)
/* check presence of non-flagged NANs in time, signal and mask tables */
{

	long ns_test=0;
	double *signal;
	int *flag;
	double *time;
	long nan_found=0;

	// check nans in mask image
	for(long jj=0;jj<ndet;jj++){
		if(read_flag_from_fits(fname, det[jj], flag, ns))
			return 1;
		for(int kk=0;kk<ns;kk++){
			if(isnan(flag[kk])){
				cout << "Warning <! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
				nan_found++;
			}
		}
		delete [] flag;
	}

	// check nans in time image
	if(read_time_from_fits(fname, time, ns_test))
		return 1;

	for(long jj=0;jj<ns_test;jj++){
		if(isnan(time[jj])&&(flag[jj]==0)){
			cout << "Warning ! a NAN has been found in \"time\" table for sample n° " << jj << endl;
			nan_found++;
		}
	}

	delete [] time;


	// check nans in signal
	for(int ii=0;ii<ndet;ii++){
		if(read_signal_from_fits(fname, det[ii], signal,ns_test))
			return 1;
		for(long jj=0;jj<ns_test;jj++){
			if(isnan(signal[jj])&&(flag[jj]==0)){
				cout << "Warning <! a NAN has been found in \"signal\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				nan_found++;
			}
		}
		delete [] signal;
	}

	return nan_found;
}

long check_NAN_altpositionHDU(string fname,long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it)
/* check non-flagged NANs in LON/LAT HIPE format */
{


	long ns_test=0;
	double *lon;
	double *lat;
	int *flag;
	long nan_found=0;


	for(int ii=0;ii<ndet;ii++){
		if(read_flag_from_fits(fname, det[ii], flag, ns_test)) // read mask image
			return 1;

		if(read_LON_LAT_from_fits(fname, det[ii], lon, lat, ns_test))
			return 1;

		for(long jj=0;jj<ns_test;jj++){
			if(isnan(lon[jj])&&(flag[jj]==0)){ // check for NANs in LON
				cout << "Warning <! a NAN has been found in \"lon\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				nan_found++;
			}
			if(isnan(lat[jj])&&(flag[jj]==0)){ // check for NANs in LAT
				cout << "Warning <! a NAN has been found in \"lat\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				nan_found++;
			}
		}
		delete [] lon;
		delete [] lat;
		delete [] flag;

	}

	return nan_found;
}

bool check_bolos(std::vector<string> bolo_fits_vect, std::vector<string> bolo_fits_0_vect)
/* Check that a given fits file bolometer list is exactly the same as the first input fits file one */
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

int check_Flag(string fname, std::vector<std::string> det, long ndet, long ns, double *percent_tab, long &init_flag_num, long &end_flag_num)
/*  Lookfor fully or more than 80% flagged detectors, also flag singletons */
{

	int *flag; // mask table
	long sum, remove_start, revert_start;

	std::vector<long> init_tab;
	std::vector<long> end_tab;

	for(int jj=0;jj<ndet;jj++){

		sum=0;
		remove_start=0;
		revert_start=ns;

		if(read_flag_from_fits(fname, det[jj], flag, ns))
			return 1;

		while((flag[remove_start]!=0) && (remove_start<ns)){
			remove_start++;
			sum++;
		}

		while((flag[revert_start-1]!=0) && (revert_start>0)){
			revert_start--;
		}


		for(long ii = remove_start; ii<ns; ii++)
			if(flag[ii]!=0)
				sum++;

		init_tab.push_back(remove_start);
		end_tab.push_back((ns-revert_start));

		percent_tab[jj]= ((double) sum*100.) / ns;

		delete [] flag;
	}

	std::vector<long>::iterator it;

	it=min_element(init_tab.begin(),init_tab.end());
	init_flag_num = (*it);

	it=min_element(end_tab.begin(),end_tab.end());
	end_flag_num = (*it);


	return 0;

}

int computeSpeed(std::string fname , unsigned int testExt, long ns, string field, double *& speed){

	double *time, *lon, *lat, *phi;
	long ns_dummy;

	double *x, *y, *dummy1, *dummy2;
	int *stat;

	// We know that we should have a time table... but let check again...
	if (( testExt & EXT_TIME) == EXT_TIME){
		read_time_from_fits(fname, time, ns_dummy);
		if (ns_dummy != ns)
			return BAD_HDU_NUM;
	} else {
		return BAD_HDU_NUM;
	}

	// read the position
	switch(testExt & EXT_POS){
	case SANEPIC_FORMAT:
		read_ReferencePosition_from_fits(fname, lon, lat, phi, ns_dummy);
		if (ns_dummy != ns)
			return BAD_HDU_NUM;
		break;
	case BOTH_FORMAT:
		// default to the LON/LAT case if both extensions are present
		/*no_break*/
	case HIPE_FORMAT:
		read_LON_LAT_from_fits(fname, field, lon, lat, ns_dummy);
		if (ns_dummy != ns)
			return BAD_HDU_NUM;
		break;
	default:
		return BAD_HDU_NUM;
	}

	// Compute the projected speed of the lon/lat vectors detector...
	// ... first define a projection...
	struct celprm celestialProj;
	celini(&celestialProj);

	// Default projection ton TAN...
	tanset(&celestialProj.prj);

	// Dummy reference
	celestialProj.ref[0] = lon[ns/2];
	celestialProj.ref[1] = lat[ns/2];

	if (celset(&celestialProj)) {
		cerr << "problem celset\n";
		return 1;
	}

	// ... Project lon/lat onto x/y (flat sky)
	stat = new int[ns];
	x = new double[ns];
	y = new double[ns];
	dummy1 = new double[ns];
	dummy2 = new double[ns];
	cels2x(&celestialProj, ns, 0, 1, 1,lon, lat, dummy1, dummy2, x, y, stat);
	delete [] lon;
	delete [] lat;
	delete [] dummy1;
	delete [] dummy2;
	delete [] stat;
	celfree(&celestialProj);

	speed   = new double[ns-1];
	// ... Compute the speed on flat sky
	for (long ii=0; ii< ns-1; ii++)
		speed[ii] = sqrt( gsl_pow_2(x[ii+1]-x[ii]) + gsl_pow_2(y[ii+1]-y[ii]) ) / (time[ii+1]-time[ii]) * 3600; // in arcsec/s

	return 0;
}

long check_Speed(param_saneCheck saneCheck_struct, std::string fname , int testExt, long ns, string field, double &meanSpeed, double &belowSpeed, double & aboveSpeed, int *& speedFlag, long & initSpeedFlag, long & endSpeedFlag){

	long nFlagSpeed = 0;

	double * speed=NULL;

	long nGoodSpeed=0;
	double * dev_speed;
	double mad_speed;


	if (computeSpeed(fname, testExt, ns, field, speed))
		return 1;

	// ... Determine the mean value of the speed
	// ... Computing the histogram of the speed is more robust than median...

	double histo_min,   histo_max, binsize = 0.01; // arcsec/sec
	double histo_lower, histo_upper;

	gsl_stats_minmax(&histo_min, &histo_max, speed, 1, ns-1);
	size_t nbins = (histo_max-histo_min)/binsize;

	gsl_histogram * histo;
	histo = gsl_histogram_alloc (nbins);
	gsl_histogram_set_ranges_uniform (histo, histo_min, histo_max);
	for (long ii=0; ii<ns-1; ii++){
		if (speed[ii] > saneCheck_struct.thresholdSpeed){
			gsl_histogram_increment(histo, speed[ii]);
			nGoodSpeed++;
		}
	}

	// ... the mean spead being the most probable speed...
	gsl_histogram_get_range (histo, gsl_histogram_max_bin (histo), &histo_lower, &histo_upper);
	meanSpeed = (histo_lower+histo_upper)/2;


	//	string outfilename = fname+".histo";
	//	FILE *fp;
	//	fp = fopen(outfilename.c_str(),"w");
	//	gsl_histogram_fprintf (fp, histo, "%g", "%g");
	//	fclose(fp);

	gsl_histogram_free (histo);


	// Using MAD to get the deviation
	// http://en.wikipedia.org/wiki/Median_absolute_deviation
	dev_speed = new double[nGoodSpeed];
	long iGoodSpeed = 0;
	for (long ii=0; ii < ns-1; ii++){
		if (speed[ii] > saneCheck_struct.thresholdSpeed ){
			dev_speed[iGoodSpeed++] = abs(speed[ii]-meanSpeed);
		}
	}

	gsl_sort(dev_speed, 1, nGoodSpeed);
	mad_speed = gsl_stats_median_from_sorted_data (dev_speed, 1, nGoodSpeed);
	delete [] dev_speed;

	belowSpeed = meanSpeed - saneCheck_struct.kappaSpeed * mad_speed;
	aboveSpeed = meanSpeed + saneCheck_struct.kappaSpeed * mad_speed;

	// Default to the provided values, if any...
	if (saneCheck_struct.belowSpeed != -1)
		belowSpeed = saneCheck_struct.belowSpeed;

	if (saneCheck_struct.aboveSpeed != -1)
		aboveSpeed = saneCheck_struct.aboveSpeed;

	// Build the flag based on the belowSpeed and aboveSpeed


	for (long ii=0; ii < ns-1; ii++){
		speedFlag[ii] = (int) ( speed[ii] < belowSpeed || speed[ii] > aboveSpeed );
		if (speedFlag[ii] != 0)
			nFlagSpeed++;
	}

	if ( nFlagSpeed > 0)
		// By default flag the last element
		speedFlag[ns-1] = 1;


	// Count the number of flagged element at the beginning and end of the data...
	initSpeedFlag = 0;
	endSpeedFlag = 0;
	while(speedFlag[++initSpeedFlag] != 0);
	while(speedFlag[ns-1-(++endSpeedFlag)] != 0);



	delete [] speed;
	return nFlagSpeed;

}


int check_time_gaps(string fname,long ns, double fsamp, std::vector<long> &indice, double &Populated_freq, struct checkHDU check_it)
/* check for time gaps in time table */
{

	long ns_dummy;
	double *time,*diff;
	double sum=0.0;
	std::vector<double> freq; // The whole frequency that can be extracted from the time vector

	if(read_time_from_fits(fname, time, ns_dummy))
		return 1;

	if (ns_dummy != ns)
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

	Populated_freq = freq[ind]; // the most populated frequency only is kept

	double zero_cinq_pourcent=0.5*Populated_freq/100; // compare to the user frequency given in ini file

	if(abs(Populated_freq-fsamp)/fsamp>zero_cinq_pourcent){
		cout << "\nWarning, the sampling frequency you have mentioned in the ini file seems to be wrong : \n";
		cout << "ini file : " << fsamp << " != " << Populated_freq << endl;
	}

	for(long jj=0;jj<ns-1;jj++){ // locate time gaps
		if((abs(diff[jj])>1.9/Populated_freq) || (diff[jj]<0)){ //||(abs(diff[jj])<1/Populated_freq/1.95)
#ifdef DEBUG
			cout << "WARNING ! Time gap at " << jj << " (" << time[jj] <<") : " << fixed <<  setprecision(8) << diff[jj] << endl;
#endif
			indice.push_back(jj); // store sample indice : where the time gap is
		}
	}

	sum=0.0;
	// recompute real frequency
	for(long jj=0;jj<ns-1;jj++){
		if(((int) count (indice.begin(), indice.end(), jj))==0)
			sum+=diff[jj];
	}


#ifdef DEBUG

	double real_freq= (double)(ns-1-(long)indice.size())/sum; // dont take the gaps into account

	// print to std
	cout << "ini file fsamp : " << fsamp << " Most Populated freq : " << Populated_freq << " Recomputed real freq : " << real_freq << endl;
	cout << endl;
#endif

	// clean up
	delete [] time;
	delete [] diff;
	delete [] counter;

	return 0;

}

int check_bolo_gain(string fname,long ns, string bolo_gain_filename, std::vector<std::string> det, long ndet, struct checkHDU check_it){


	long ns_test=0;
	double *signal_tot, *signal_samp;
	double *signal;
	signal_tot=new double[ns];
	signal_samp=new double[ndet];
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
		fill(signal_samp,signal_samp+ndet,0.0);

		cout << "before read\n";
		if(read_sample_signal_from_fits(fname, jj+1, signal_samp, det, ndet))
			return 1;
		cout << "before alloc\n";
		std::vector<double> signal_vec(signal_samp, signal_samp+ndet);
		cout << signal_vec[0] << endl;
		signal_tot[jj]=median(signal_vec);
		cout << signal_tot[jj] << endl;
	}

	fill(signal_samp,signal_samp+ndet,0.0);

	for(long jj=0;jj<10;jj++)
		cout << signal_tot[jj] << " " ;
	getchar();
	for(long idet=0;idet<ndet;idet++){
		std::vector<double> signal_vec;
		read_signal_from_fits(fname, det[idet], signal, ns_test);

		for(long ii=0;ii<ns;ii++)
			signal_vec.push_back(signal_tot[ii]/signal[ii]);

		signal_samp[idet]=median(signal_vec); // store gain in signal_samp

		//		signal_vec.clear();
		delete [] signal;
	}

	cout << "gain \n";
	for(long ii=0;ii<ndet;ii++)
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

//TODO: Put all that into a fits file

int print_to_bin_file(std::string tmp_dir, std::string filename, long init_flag, long end_flag, double Populated_freq, long ns, int * speedFlag, std::vector<long> indice){

	string fname = tmp_dir + FitsBasename(filename) + "_saneFix_indices.bin"; // output saneFix logfile filename

	std::ofstream file; // saneFix indice log file
	file.open(fname.c_str(), ios::out | ios::trunc);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		return 1;
	}else{
		// store init num of flag to remove for saneFix
		file << init_flag << endl;

		// store last num of flag to remove for saneFix
		file << end_flag << endl;

		// store real_freq for saneFix
		file << Populated_freq << endl;

		file << ns << endl;
		for (long ii=0; ii < ns; ii++)
			file << speedFlag[ii] << " ";

		for(long ii = 0; ii<(long)indice.size();ii++)// store indices for saneFix
			file << indice[ii] << " ";
	}

	file.close();

	return 0;
}


int exportCheckToFits(std::string tmp_dir, std::string filename, long init_flag, long end_flag, double Populated_freq, long ns, int * speedFlag, std::vector<long> indice){

	string fname = "!" + tmp_dir + FitsBasename(filename) + "_saneFix_indices.fits"; // output saneFix logfile filename

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...


	long naxes[] = { 0 }; // size of dimensions
	long fpixel[] = { 1 }; // index for write_pix


	// create fits file
	if (fits_create_file(&fp, fname.c_str(), &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	// Dummy image
	naxes[0] = 0;
	if (fits_create_img(fp, DOUBLE_IMG, 0, naxes, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}

	// Add init_flag, end_flag, populated_freq
	if (fits_write_key(fp, TLONG, (const char*) "initFlag", (long *) & init_flag, (char *) "Number of sample to flag at the beginning", &fits_status)) {
		fits_report_error(stderr, fits_status); return 1;
	}
	if (fits_write_key(fp, TLONG, (const char*) "enFlag", (long *) & end_flag, (char *) "Number of sample to flag at the end", &fits_status)) {
		fits_report_error(stderr, fits_status); return 1;
	}
	if (fits_write_key(fp, TDOUBLE, (const char*) "fsamp", (double *) & Populated_freq, (char *) "[Hz] Sampling frequency", &fits_status)) {
		fits_report_error(stderr, fits_status); return 1;
	}
	// ----------------------------------------------------------------

		naxes[0] = ns;
		if (fits_create_img(fp, SHORT_IMG, 1, naxes, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}


		if (fits_update_key(fp, TSTRING, (char *) "EXTNAME", (char *) "speedFlag", (char *) "table name", &fits_status))
			return 1;

		if (fits_write_pix(fp, TINT, fpixel, ns, (int *) speedFlag, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		// write date to file
		if (fits_write_date(fp, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		if (fits_write_chksum(fp, &fits_status)) {
			cout << "error checksum !\n";
			return 1;
		}

		// ----------------------------------------------------------------

		naxes[0] = (long) indice.size();
		if (fits_create_img(fp, LONG_IMG, 1, naxes, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}


		if (fits_update_key(fp, TSTRING, (char *) "EXTNAME", (char *) "gapsIndices", (char *) "table name", &fits_status))
			return 1;

		if (fits_write_pix(fp, TLONG, fpixel, naxes[0], (long *) &(indice)[0], &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
		// write date to file
		if (fits_write_date(fp, &fits_status)) {
			fits_report_error(stderr, fits_status);
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


void log_gen(long  *bolo_, string outname, std::vector<std::string> det, long ndet, double *percent_tab)
/* generating log files for user information */
{


	FILE *fp;
	long tot=0;

	fp=fopen(outname.c_str(),"w");

	for(long ii=0; ii< ndet; ii++)
		if(bolo_[ii]>0){
			string temp = det[ii]; // copy the name of the wrong or bad detector in this ascii log file
			if(percent_tab!=NULL)
				fprintf(fp, "%s %lf%%\n", (char*)temp.c_str(),percent_tab[ii]);
			else
				fprintf(fp, "%s\n", (char*)temp.c_str());
			tot++;
		}

#ifdef DEBUG
	if(tot>0)
		cout << "There is/are " << tot << " bolometer(s) in this file !\n\n";
#endif

	fclose(fp);

}

int read_sample_signal_from_fits(string filename, int sample, double *& signal_samp, std::vector<std::string> det, long ndet){

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int naxis = 0, anynul;
	//	long ns;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", 0, &status)){
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
	//	ns = naxes[0];
	double temp;

	for(int idet=1;idet<=(int)ndet;idet++){
		// ---------------------------------------------
		// Retrieve the corresponding pix
		fpixel[0] = sample;
		string field = det[idet];
		fpixel[1] = find_channel_index(fptr, (char *)field.c_str());
		if (fits_read_pix(fptr, TDOUBLE, fpixel, 1, 0, &temp, &anynul, &status)){
			fits_report_error(stderr, status);
			return 1;
		}
		signal_samp[idet]=temp;
		cout << signal_samp[idet] << " " << temp << endl;
	}

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status)){
		cout << status << endl;
		fits_report_error(stderr, status);
		return 1;
	}

	cout << "after close\n";

	return 0;

}
