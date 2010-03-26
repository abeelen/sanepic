

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


//int check_path(string strPath, string path_type){
//
//
//
//	if ( access( strPath.c_str(), 0 ) == 0 )
//	{
//		struct stat status;
//		stat( strPath.c_str(), &status );
//
//		if ( status.st_mode & S_IFDIR )
//		{
//			cout << "The directory " << path_type << " : " << strPath << " exists." << endl;
//			return 0;
//		}
//		else
//		{
//			cout << "Warning : The path " << path_type << " : " << strPath << " is a file." << endl;
//			return 1;
//		}
//	}
//	else
//	{
//		cout << "Warning : Path " << path_type << " : " << strPath << " doesn't exist." << endl;
//		return 1;
//	}
//
//
//
//}

int who_do_it(int size, int rank, int ii){

	if(size==1)
		return 0;

	if(size>=ii)
		return ii;

	if(size<ii){
		while(ii>size)
			ii=ii-size;
		return ii;
	}

	return -1;
}

void check_detector_is_in_fits(struct detectors det,struct detectors bolo_fits, string filename){

	int mycount=0;

	for(int jj=0;jj< det.ndet; jj++){
		mycount = (int) count (bolo_fits.boloname.begin(), bolo_fits.boloname.end(), det.boloname[jj]);
		if(mycount==0)
			cout << "Warning ! The detector " << det.boloname[jj] << " is not referenced in the fits " << filename << endl;
	}


}

void check_positionHDU(string fname,long ns,struct detectors det, int format){

	fitsfile *fptr;
	int status = 0;
	int colnum;
	long ns_test=0;
	long ndet_test=0;


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"reference position\" was not found, or his Type should be Binary table" << endl;
	}else{

		fits_get_num_rows(fptr, &ns_test, &status);
		if(ns!=ns_test){
			cout << "\"reference position\" has a wrong number of rows (must be equal to ns : " << ns << " )" << endl;
		}

		fits_get_num_cols(fptr, &colnum, &status);
		if(colnum!=3){
			cout << "\"reference position\" has a wrong number of cols (must be equal to 3 : RA, DEC, PHI )" << endl;
		}else{


			fits_get_colnum(fptr, CASEINSEN, (char*) "RA", &colnum, &status);
			if(colnum!=1){
				cout << "\"RA\" was not found in \"reference position\"" << endl;
			}

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) "DEC", &colnum, &status);
			if(colnum!=2){
				cout << "\"DEC\" table was not found in \"reference position\"" << endl;
			}

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) "PHI", &colnum, &status);
			if(colnum!=3){
				cout << "\"PHI\" table was not found in \"reference position\"" << endl;
			}
		}
	}

	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"offsets\" was not found, or his Type should be Binary table" << endl;
	}else{

		fits_get_num_rows(fptr, &ndet_test, &status);
		if(det.ndet!=ndet_test){
			cout << "\"offsets\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
		}

		colnum=0;
		fits_get_num_cols(fptr, &colnum, &status);
		if(colnum!=3){
			cout << "\"offsets\" has a wrong number of cols (must be equal to 3 : NAME, DX, DY )" << endl;
		}else{
			string name_table,dx_table,dy_table;

			if(format==1){
				name_table="names";
				dx_table="dX";
				dy_table="dY";
			}else{
				name_table="NAME";
				dx_table="DX";
				dy_table="DY";
			}


			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) name_table.c_str(), &colnum, &status);
			if(colnum!=1){
				cout << "\"NAME\" table was not found in \"offsets\"" << endl;
			}

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) dx_table.c_str(), &colnum, &status);
			if(colnum!=2){
				cout << "\"DX\" table was not found in \"offsets\"" << endl;
			}

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) dy_table.c_str(), &colnum, &status);
			if(colnum!=3){
				cout << "\"DY\" table was not found in \"offsets\"" << endl;
			}
		}
	}

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

void check_commonHDU(string fname,long ns,struct detectors det){

	fitsfile *fptr;
	int status = 0;
	int colnum;
	long ndet_test=0;
	int naxis=0;
	long naxes[2] = { 1, 1 };


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);


	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"channels\" was not found, or his Type should be Binary table" << endl;
	}else{

		ndet_test=0;
		fits_get_num_rows(fptr, &ndet_test, &status);
		if(ndet_test!=det.ndet){
			cout << "\"channels\" has a wrong number of rows (must be equal to ndet : " << det.ndet << " )" << endl;
		}

		colnum=0;
		fits_get_num_cols(fptr, &colnum, &status);
		if(colnum!=1){
			cout << "\"channels\" has a wrong number of cols (must be equal to 1 : NAMES )" << endl;
		}else{

			colnum=0;
			fits_get_colnum(fptr, CASEINSEN, (char*) "NAMES", &colnum, &status);
			if(colnum!=1){
				cout << "\"NAMES\" table was not found in \"channels\"" << endl;
			}
		}
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"time\" was not found, or his Type should be image" << endl;
	}


	if(fits_get_img_dim(fptr, &naxis, &status)){
		fits_report_error(stderr, status);
	}else{

		if(naxis!=1){
			cout << "\"time\" has a wrong number of columns, should be equal to 1 " << endl;
		}else{
			if (fits_get_img_size(fptr, 2, naxes, &status))
				fits_report_error(stderr, status);
			if(naxes[0]!=ns)
				cout << "\"time\" has a wrong number of elements, should be equal to ns : " << ns << endl;
		}
	}


	naxes[0]=1;
	naxes[1]=1;

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"signal\" was not found, or his Type should be image" << endl;
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status))
			fits_report_error(stderr, status);
		if(naxis != 2){
			fits_report_error(stderr,BAD_NAXIS);
			cout << "\"signal\" must have 2 dimensions" << endl;
		}
		if (fits_get_img_size(fptr, 2, naxes, &status))
			fits_report_error(stderr, status);
		else{

			if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
				cout << "\"signal\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
			}
		}
	}

	naxes[0]=1;
	naxes[1]=1;


	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"mask\" was not found, or his Type should be image" << endl;
		//		exit(0);
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status))
			fits_report_error(stderr, status);
		if(naxis != 2){
			fits_report_error(stderr,BAD_NAXIS);
			cout << "\"mask\" must have 2 dimensions" << endl;
		}else{

			if (fits_get_img_size(fptr, 2, naxes, &status))
				fits_report_error(stderr, status);

			//	cout << naxes[0] << " " << naxes[1] << endl;
			if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
				cout << "\"mask\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
				//		exit(0);
			}
		}
	}
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}


void check_altpositionHDU(string fname,long ns,struct detectors det){

	fitsfile *fptr;
	int status = 0;
	int naxis=0;
	long naxes[2] = { 1, 1 };

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"ra\" was not found, or his Type should be image" << endl;
		//		exit(0);
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status))
			fits_report_error(stderr, status);
		if(naxis != 2){
			fits_report_error(stderr,BAD_NAXIS);
			cout << "\"ra\" must have 2 dimensions" << endl;
		}
		if (fits_get_img_size(fptr, 2, naxes, &status))
			fits_report_error(stderr, status);
		else{

			//	cout << naxes[0] << " " << naxes[1] << endl;
			if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
				cout << "\"ra\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
				//		exit(0);
			}
		}
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "\"dec\" was not found, or his Type should be image" << endl;
		//		exit(0);
	}else{

		if (fits_get_img_dim(fptr, &naxis, &status))
			fits_report_error(stderr, status);
		if(naxis != 2){
			fits_report_error(stderr,BAD_NAXIS);
			cout << "\"dec\" must have 2 dimensions" << endl;
		}
		if (fits_get_img_size(fptr, 2, naxes, &status))
			fits_report_error(stderr, status);
		else{

			//	cout << naxes[0] << " " << naxes[1] << endl;
			if((naxes[0]!=ns)&&(naxes[1]!=det.ndet)){
				cout << "\"dec\" has a wrong size, it must be ns*ndet : " << ns << " x " << det.ndet << endl;
				//		exit(0);
			}
		}
	}

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);
}

void check_NAN_positionHDU(string fname,long ns,struct detectors det){

	long ns_test=0;
	double *ra;
	double *dec,*phi;
	double **offsets;

	int *flag;



	//	// check nans in flag
	//	for(long jj=0;jj<det.ndet;jj++){
	//		//		cout << "loop " << jj << endl;
	//		read_flag_from_fits(fname, det.boloname[jj], flag, ns);
	//		for(int kk=0;kk<ns;kk++){
	//			if(flag[kk]==NAN){
	//				cout << "Warning ! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
	//				//				exit(0);
	//			}
	//			if(isnan(flag[kk])){
	//				cout << "Warning <! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
	//				//				exit(0);
	//			}
	//		}
	//
	//	}

	for(int ii=0;ii<det.ndet;ii++){
		read_flag_from_fits(fname, det.boloname[ii], flag, ns);

		read_ReferencePosition_from_fits(fname, ra, dec, phi, ns_test);

		//		read_ra_from_fits(fname , det.boloname[ii], ra, ns_test);
		for(long jj=0;jj<ns_test;jj++){
			if(isnan(ra[jj])){
				cout << "Warning <! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
			if(isnan(dec[jj])){
				cout << "Warning <! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
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
		delete [] flag;
	}

	// flag indepedent
	read_all_bolo_offsets_from_fits(fname, det.boloname, offsets);

	for(int ii=0;ii<det.ndet;ii++)
		if(isnan(offsets[ii][0])||isnan(offsets[ii][1])){
			cout << "Warning ! a NAN has been found in \"offsets\" table for bolometer n° " << ii << endl;
			cout << "You should not take this detector for the computation of Sanepic\n";
			//			exit(0);
		}


	free_dmatrix(offsets,(long)0,det.ndet-1,(long)0,2-1);


}

void check_NAN_commonHDU(string fname,long ns,struct detectors det){

	long ns_test=0;
	double *signal;
	int *flag;
	double *time;


	// check nans in flag
	for(long jj=0;jj<det.ndet;jj++){
		//		cout << "loop " << jj << endl;
		read_flag_from_fits(fname, det.boloname[jj], flag, ns);
		for(int kk=0;kk<ns;kk++){
			if(isnan(flag[kk])){
				cout << "Warning <! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
				//				exit(0);
			}
		}
		delete [] flag;
	}

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
			if(isnan(signal[jj])){
				cout << "Warning <! a NAN has been found in \"signal\" table for bolometer n° " << ii << " sample n° " << jj << endl;
				//				exit(0);
			}
		}
		delete [] signal;
	}


}

void check_NAN_altpositionHDU(string fname,long ns,struct detectors det){


	long ns_test=0;
	double *ra;
	double *dec;
	int *flag;


	for(int ii=0;ii<det.ndet;ii++){
		read_flag_from_fits(fname, det.boloname[ii], flag, ns_test);
		//		read_ra_from_fits(fname, det.boloname[ii], ra, ns_test);
		//		read_dec_from_fits(fname, det.boloname[ii], dec, ns_test);
		read_ra_dec_from_fits(fname, det.boloname[ii], ra, dec, ns_test);
		for(long jj=0;jj<ns_test;jj++){
			if(isnan(ra[jj])){
				cout << "Warning <! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
			if(isnan(dec[jj])){
				cout << "Warning <! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
		}
		delete [] ra;
		delete [] dec;
		delete [] flag;

	}

}

void check_NaN(string fname,long ns,struct detectors det){




	long ns_test=0;
	double *signal;
	int *flag;
	double *ra;
	double *dec,*phi;
	double **offsets;
	double *time;


	for(int ii=0;ii<det.ndet;ii++){
		read_ReferencePosition_from_fits(fname, ra, dec, phi, ns_test);
		for(long jj=0;jj<ns_test;jj++){
			if(isnan(ra[jj])){
				cout << "Warning <! a NAN has been found in \"ra\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
			if(isnan(dec[jj])){
				cout << "Warning <! a NAN has been found in \"dec\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
			if(isnan(phi[jj])){
				cout << "Warning <! a NAN has been found in \"phi\" table for bolometer n° " << ii << " sample n° " << jj << endl;
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
		}


	free_dmatrix(offsets,(long)0,det.ndet-1,(long)0,2-1);

	// check nans in time constant
	read_time_from_fits(fname, time, ns);
	for(long jj=0;jj<ns;jj++){
		if(isnan(time[jj])){
			cout << "Warning ! a NAN has been found in \"time\" table for sample n° " << jj << endl;
		}
	}

	delete [] time;




	// check nans in signal
	for(int ii=0;ii<det.ndet;ii++){
		read_signal_from_fits(fname, det.boloname[ii], signal,ns_test);
		for(long jj=0;jj<ns_test;jj++){
			if(isnan(signal[jj])){
				cout << "Warning <! a NAN has been found in \"signal\" table for bolometer n° " << ii << " sample n° " << jj << endl;
			}
		}
		delete [] signal;
	}



	// check nans in flag
	for(long jj=0;jj<det.ndet;jj++){
		read_flag_from_fits(fname, det.boloname[jj], flag, ns);
		for(int kk=0;kk<ns;kk++){
			if(isnan(flag[kk])){
				cout << "Warning <! there is a NaN in the \"flag\" field of bolometer n° " << jj << " sample n° " << kk << endl;
			}
		}
		delete [] flag;
	}





}

bool check_bolos(std::vector<string> bolo_fits_vect, std::vector<string> bolo_fits_0_vect){


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

void check_flag(string fname,struct detectors det,long ns, string outname,long *&bolos_global,long *&bolos_global_80){

	int *flag;
	short sum=0;

	int marge = 20;
	long ii=1;
	long tt=0;
	long rr=0;



	for(int jj=0;jj<det.ndet;jj++){

		ii=1;
		sum=0;
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
				flag[ii]=1;
				cout << "singleton trouvé en " << det.boloname[jj] << " sample n° " << ii << endl;
			}

			sum+=flag[ii];
			ii++;

		}

		if(sum==ns){
			cout << "Warning ! " << det.boloname[jj] << " is totally flagged" << endl;
			bolos_global[jj]=1;
		}else{
			if(sum>80*ns/100){
				cout << "Warning ! " << det.boloname[jj] << " is more than 80% flagged" << endl;
				bolos_global_80[jj]=1;
			}else if(sum>50*ns/100)
				cout << "Warning ! " << det.boloname[jj] << " is more than 50% flagged" << endl;
		}

		delete [] flag;
	}


}


void check_time_gaps(string fname,long ns, double fsamp, struct common dir){


	std::vector<long> indice;

	std::ostringstream oss;
	oss << fname;
	string filename = oss.str();
	string fname2 = dir.tmp_dir + Basename(filename) + "_saneFix_indices.bin";


	double *time,*diff;
	double sum=0.0, mean=0.0;


	read_time_from_fits(fname, time, ns);
	diff = new double [ns-1];

	for(long jj=0;jj<ns-1;jj++){
		diff[jj]=(time[jj+1]-time[jj]);
		sum+=diff[jj];

	}
	mean = sum/(double)(ns-1);
	cout << "time gaps in mean :" << fixed << setprecision(15) <<  mean << endl;
	cout << "ie. sampling frequency ~ " << 1/mean << " and fsamp = " << fsamp << endl;

	double zero_cinq_pourcent=0.5/mean/100;
	if(abs(1/mean-fsamp)/fsamp>zero_cinq_pourcent){
		cout << "Warning, the sampling frequency you have mentioned in the ini file seems to be wrong : \n";
		cout << "ini file : " << fsamp << " != " << 1/mean << endl;
	}


	for(long jj=0;jj<ns-1;jj++){
		//		if((jj>14800)&&(jj<14900))
		//			cout << fixed << setprecision(15) << diff[jj] << endl;
		if((abs(diff[jj])>1.9*mean)||(abs(diff[jj])<mean/1.9)){
			cout << "WARNING ! At sample " << jj << " there is a gaps in the time constant : " << fixed <<  setprecision(8) << diff[jj] << endl;
			indice.push_back(jj);
		}
	}

	sum=0.0;
	// recompute real frequency
	for(long jj=0;jj<ns-1;jj++){
		if(((int) count (indice.begin(), indice.end(), jj))==0)
			sum+=diff[jj];
	}
	double real_freq= (double)(ns-1-(long)indice.size())/sum;

	std::ofstream file;
	file.open(fname2.c_str(), ios::out | ios::trunc);
	if(!file.is_open()){
		cerr << "File [" << fname2 << "] Invalid." << endl;
	}else{
		// store real_freq for saneFix
		file << real_freq << " ";

		for(long ii = 0; ii<(long)indice.size();ii++)
			file << indice[ii] << " ";
	}

	file.close();

	cout << endl;

	delete [] time;
	delete [] diff;


}


void log_gen(long  *bolo_, string outname, struct detectors det){


	FILE *fp;
	long tot=0;

	fp=fopen(outname.c_str(),"w");

	for(long ii=0; ii< det.ndet; ii++)
		if(bolo_[ii]>0){
			string temp = det.boloname[ii];
			fprintf(fp, "%s\n", (char*)temp.c_str());
			tot++;
		}
	if(tot>0)
		cout << "There are " << tot << " bolometers in this file !\n";

	fclose(fp);

}

