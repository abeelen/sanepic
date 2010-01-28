

#include "covMatrixIO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}



using namespace std;


void read_bolo_list(string fname, struct detectors &det){

	fitsfile *fptr;
	int status = 0;
	//	int colnum;
	long ndet_test=0;
	char **temp_bolo;

	det.boloname.clear();

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	read_channels(fptr,temp_bolo, ndet_test);

	//	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
	//		fits_report_error(stderr, status);
	//
	//
	//	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status)){
	//		fits_report_error(stderr, status);
	//		cout << "\"channels\" was not found, or his Type should be Binary table" << endl;
	//		//		exit(0);
	//	}
	//
	//	colnum=0;
	//	fits_get_num_cols(fptr, &colnum, &status);
	//	if(colnum!=1){
	//		cout << "\"channels\" has a wrong number of cols (must be equal to 1 : NAMES )" << endl;
	//		//		exit(0);
	//	}
	//
	//	ndet_test=0;
	//	fits_get_num_rows(fptr, &ndet_test, &status);
	//	if(ndet_test<=0){
	//		cout << "\"channels\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
	//		exit(0);
	//	}
	//	det.ndet=ndet_test;
	//
	//	colnum=0;
	//	fits_get_colnum(fptr, CASEINSEN, (char*) "NAMES", &colnum, &status);
	//	if(colnum!=1){
	//		cout << "\"NAMES\" table was not found in \"channels\"" << endl;
	//		exit(0);
	//	}
	//
	//	temp_bolo = new string[ndet_test];
	//
	//	fits_read_col(fptr, TSTRING, colnum, 1, 1, ndet_test, NULL, temp_bolo, 0, &status);


	for (long ii=0; ii < ndet_test; ii++){
//		cout << temp_bolo[ii] << endl;
		det.boloname.push_back(temp_bolo[ii]);
	}

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	delete [] temp_bolo;
}

void check_hdu(string fname,long ns,struct detectors det){

	fitsfile *fptr;
	int status = 0;
	int colnum;
	long ns_test=0;
	long ndet_test=0;
	int naxis=0;
	long naxes[2] = { 1, 1 };


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"reference position\" was not found, or his Type should be Binary table" << endl;
		//		exit(0);
	}

	fits_get_num_rows(fptr, &ns_test, &status);
	if(ns!=ns_test){
		cout << "\"reference position\" has a wrong number of rows (must be equal to ns : " << ns << " )" << endl;
		//		exit(0);
	}

	fits_get_num_cols(fptr, &colnum, &status);
	if(colnum!=3){
		cout << "\"reference position\" has a wrong number of cols (must be equal to 3 : RA, DEC, PHI )" << endl;
		//		exit(0);
	}


	fits_get_colnum(fptr, CASEINSEN, (char*) "RA", &colnum, &status);
	if(colnum!=1){
		cout << "\"RA\" was not found in \"reference position\"" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_colnum(fptr, CASEINSEN, (char*) "DEC", &colnum, &status);
	if(colnum!=2){
		cout << "\"DEC\" table was not found in \"reference position\"" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_colnum(fptr, CASEINSEN, (char*) "PHI", &colnum, &status);
	if(colnum!=3){
		cout << "\"PHI\" table was not found in \"reference position\"" << endl;
		//		exit(0);
	}

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"offsets\" was not found, or his Type should be Binary table" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_num_cols(fptr, &colnum, &status);
	if(colnum!=3){
		cout << "\"offsets\" has a wrong number of cols (must be equal to 3 : NAME, DX, DY )" << endl;
		//		exit(0);
	}

	fits_get_num_rows(fptr, &ndet_test, &status);
	if(det.ndet!=ndet_test){
		cout << "\"offsets\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_colnum(fptr, CASEINSEN, (char*) "NAME", &colnum, &status);
	if(colnum!=1){
		cout << "\"NAME\" table was not found in \"offsets\"" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_colnum(fptr, CASEINSEN, (char*) "DX", &colnum, &status);
	if(colnum!=2){
		cout << "\"DX\" table was not found in \"offsets\"" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_colnum(fptr, CASEINSEN, (char*) "DY", &colnum, &status);
	if(colnum!=3){
		cout << "\"DY\" table was not found in \"offsets\"" << endl;
		//		exit(0);
	}


	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"channels\" was not found, or his Type should be Binary table" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_num_cols(fptr, &colnum, &status);
	if(colnum!=1){
		cout << "\"channels\" has a wrong number of cols (must be equal to 1 : NAMES )" << endl;
		//		exit(0);
	}

	ndet_test=0;
	fits_get_num_rows(fptr, &ndet_test, &status);
	if(ndet_test!=det.ndet){
		cout << "\"channels\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
		//		exit(0);
	}

	colnum=0;
	fits_get_colnum(fptr, CASEINSEN, (char*) "NAMES", &colnum, &status);
	if(colnum!=1){
		cout << "\"NAMES\" table was not found in \"channels\"" << endl;
		//		exit(0);
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"time\" was not found, or his Type should be image" << endl;
		//		exit(0);
	}


	if(fits_get_img_dim(fptr, &naxis, &status)){
		fits_report_error(stderr, status);
	}
	if(naxis!=1){
		cout << "\"time\" has a wrong number of elements, should be equal to 1 " << endl;
		//		exit(0);
	}


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"signal\" was not found, or his Type should be image" << endl;
		//		exit(0);
	}

	if (fits_get_img_dim(fptr, &naxis, &status))
		fits_report_error(stderr, status);
	if(naxis != 2){
		fits_report_error(stderr,BAD_NAXIS);
		cout << "\"signal\" must have 2 dimensions" << endl;
	}
	if (fits_get_img_size(fptr, 2, naxes, &status))
		fits_report_error(stderr, status);

	//	cout << naxes[0] << " " << naxes[1] << endl;
	if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
		cout << "\"signal\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
		//		exit(0);
	}

	naxes[0]=1;
	naxes[1]=1;


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"mask\" was not found, or his Type should be image" << endl;
		//		exit(0);
	}

	if (fits_get_img_dim(fptr, &naxis, &status))
		fits_report_error(stderr, status);
	if(naxis != 2){
		fits_report_error(stderr,BAD_NAXIS);
		cout << "\"mask\" must have 2 dimensions" << endl;
	}
	if (fits_get_img_size(fptr, 2, naxes, &status))
		fits_report_error(stderr, status);

	//	cout << naxes[0] << " " << naxes[1] << endl;
	if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
		cout << "\"mask\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
		//		exit(0);
	}


	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}


void check_NaN(string fname,long ns,struct detectors det){




	//	fitsfile *fptr;
	//	int status = 0;
	//	int colnum;
	long ns_test=0;
	//	long ndet_test=0;
	//	int naxis=0;
	//	long naxes[2] = { 1, 1 };
	double *signal;
	short *flag;
	double *ra;
	double *dec,*phi;
	double **offsets;
	double *time;


	for(int ii=0;ii<det.ndet;ii++){
		read_ReferencePosition_from_fits(fname, ra, dec, phi, ns_test);
		//		read_ra_from_fits(fname , det.boloname[ii], ra, ns_test);
		for(long jj=0;jj<ns_test;jj++){
			if(ra[jj]==NAN){
				cout << "Warning ! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(isnan(ra[jj])){
				cout << "Warning <! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(dec[jj]==NAN){
				cout << "Warning ! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(isnan(dec[jj])){
				cout << "Warning <! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(phi[jj]==NAN){
				cout << "Warning ! a NAN has been found in \"phi\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(isnan(phi[jj])){
				cout << "Warning <! a NAN has been found in \"phi\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
		}

		delete [] ra;
		delete [] dec;
		delete [] phi;
	}


	//check nans in offsets

	read_all_bolo_offsets_from_fits(fname, det.boloname, offsets);

	for(int ii=0;ii<det.ndet;ii++)
		if(isnan(offsets[ii][0])||isnan(offsets[ii][1])){
			cout << "Warning ! a NAN has been found in \"offsets\" table for bolometer n° " << ii << endl;
			//			exit(0);
		}


	free_dmatrix(offsets,(long)0,det.ndet-1,(long)0,2-1);

	//check nans in channels

	//	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status)){
	//		fits_report_error(stderr, status);
	//		cout << "\"channels\" was not found, or his Type should be Binary table" << endl;
	//		exit(0);
	//	}

	// check nans in time constant

	read_time_from_fits(fname, time, ns);
	for(long jj=0;jj<ns;jj++){
		//		cout << time[jj] << endl;
		//		getchar();
		if(isnan(time[jj])){
			cout << "Warning ! a NAN has been found in \"time\" table for sample n° " << jj << endl;
			//			exit(0);
		}
	}

	delete [] time;




	// check nans in signal

	for(int ii=0;ii<det.ndet;ii++){
		read_signal_from_fits(fname, det.boloname[ii], signal,ns_test);
		for(long jj=0;jj<ns_test;jj++){
			if(signal[jj]==NAN){
				cout << "Warning ! a NAN has been found in \"signal\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(isnan(signal[jj])){
				cout << "Warning <! a NAN has been found in \"signal\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
		}
		delete [] signal;
	}



	// check nans in flag

	for(long jj=0;jj<det.ndet;jj++){
		//		cout << "loop " << jj << endl;
		read_flag_from_fits(fname, det.boloname[jj], flag, ns);
		for(int kk=0;kk<ns;kk++){
			if(flag[kk]==NAN){
				cout << "Warning ! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
				//				exit(0);
			}
			if(isnan(flag[kk])){
				cout << "Warning <! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
				//				exit(0);
			}
		}
		delete [] flag;
	}





}

int check_flag(string fname,struct detectors det,long ns, string outname,std::vector<std::string> &bolos_global,std::vector<std::string> &bolos_global_80){

	short *flag;
	short sum=0;
	FILE * fp;

	int marge = 20;
	long ii=1;
	long tt=0;
	long rr=0;
	//	long bolo_count=0,bolo_count_80=0;
	//	std::vector<std::string> bolo_bad;
	//	std::vector<std::string> bolo_bad_80;


	fp=fopen(outname.c_str(),"w");
	//	flag = new short[ns];

	//cout << "ndet : "<< det.ndet << endl;

	for(int jj=0;jj<det.ndet;jj++){

		ii=1;
		//		cout << "det : " << det.boloname[jj] << endl;
		sum=0;
		//		cout << "jj : " << jj << endl;
		read_flag_from_fits(fname, det.boloname[jj], flag, ns);

		while(ii<ns-1){

			if((flag[ii]==0)&&(flag[ii+1]==1)&&(flag[ii-1]==1)){
				tt=ii-1;
				rr=ii+1;

				while((tt>0)&&(flag[tt]==1)&&((ii-tt)<marge+2))
					tt--;

				while((rr<ns)&&(flag[rr]==1)&&((rr-ii)<marge+2))
					rr++;
			}

			if(((rr-ii)>=marge)&&((ii-tt)>=marge)){
				cout << ii << " " << rr-ii << " " << ii-tt << endl;
				//getchar();
				flag[ii]=1;
				cout << "singleton trouvé en " << det.boloname[jj] << " sample n° " << ii << endl;
			}

			sum+=flag[ii];
			ii++;

		}

		if(sum==ns){
			cout << "Warning ! " << det.boloname[jj] << " is totally flagged" << endl;
			fprintf(fp, "%s\n", (char*)(det.boloname[jj].c_str()));
			bolos_global.push_back(det.boloname[jj]);
		}else{
			if(sum>80*ns/100){
				cout << "Warning ! " << det.boloname[jj] << " is more than 80% flagged" << endl;
				bolos_global_80.push_back(det.boloname[jj]);
			}else if(sum>50*ns/100)
				cout << "Warning ! " << det.boloname[jj] << " is more than 50% flagged" << endl;
		}

		//		cout << "fin bolo : " << det.boloname[jj] << endl;
		//		getchar();

		delete [] flag;
	}

	fclose(fp);

	return 0;
}


void check_time_gaps(string fname,long ns){


	double *time,*diff;
	double sum=0.0, mean=0.0, std=0.0, three_times_sigma=0.0;


	read_time_from_fits(fname, time, ns);
	diff = new double [ns-1];

	for(long jj=0;jj<ns-1;jj++){
		diff[jj]=(time[jj+1]-time[jj]);
		sum+=diff[jj];

	}
	mean = sum/(ns-1);
	cout << "mean :" << fixed << setprecision(15) <<  mean << endl;

	for(long jj=0;jj<ns-1;jj++)
		std+=(diff[jj]-mean)*(diff[jj]-mean);

	std=sqrt(std/(ns-1));
	cout << "sigma :" << fixed << setprecision(15) <<  std << endl;

	three_times_sigma=std*3;
	cout << "3 * sigma :" << fixed << setprecision(15) << three_times_sigma << endl << endl  << endl;


	for(long jj=0;jj<ns-1;jj++){
		//		if((jj>14800)&&(jj<14900))
		//			cout << fixed << setprecision(15) << diff[jj] << endl;
		if((diff[jj])>three_times_sigma)
			cout << "WARNING ! At sample " << jj << " there is a gaps in the time constant : " << fixed <<  setprecision(8) << (time[jj+1]-time[jj]) << endl;
	}

	cout << endl;

	//	getchar();
	delete [] time;


}


void log_gen(std::vector<string> &bolo_, string outname){


	FILE *fp;

	struct sortclass_string sortobject;
	sort(bolo_.begin(), bolo_.end(), sortobject);

	std::vector<string>::iterator it;
	std::vector<string>::iterator it2;

	// using default comparison:
	it2 = unique(bolo_.begin(), bolo_.end());

	if ((it2-bolo_.begin())>0){

		fp=fopen(outname.c_str(),"w");

		for (it=bolo_.begin(); it != it2; it++) {
			fprintf(fp, "%s\n", (char*)((*it).c_str()));
		}

		fclose(fp);
	}

}

