

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
#include <sys/types.h>  // For stat()
#include <sys/stat.h>   // For stat()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}


using namespace std;


int read_indices_file(string fname, struct common dir, std::vector<long> &indice, double &fsamp){

	std::ifstream file;
	long readed;

	std::ostringstream oss;
	oss << fname;
	string filename = oss.str();
	string fname2 = dir.tmp_dir + Basename(filename) + "_saneFix_indices.bin";

	file.open(fname2.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << fname << "] was not found." << endl;
		file.close();
		return(EXIT_FAILURE);
	}

	file >> fsamp;
	while(file >> readed)
		indice.push_back(readed);

	file.close();
	return(EXIT_SUCCESS);
}

long how_many(string fname, long ns, std::vector <long> indice, double *&time, double fsamp,  std::vector <long> &add_sample){

	long total=0;
	long gap=0, sum=0;

	read_time_from_fits(fname, time, ns);
	for(long ii=0; ii < (long)indice.size(); ii++){
		//		cout << indice[ii] << " " <<time[indice[ii]+1]-time[indice[ii]] << " " << round((time[indice[ii]+1]-time[indice[ii]])*fsamp)-1 << endl;
		gap=round((time[indice[ii]+1]-time[indice[ii]])*fsamp)-1;
		if(gap>0){
			sum+=gap;
			add_sample.push_back(gap);
		}else{
			if(gap==0){
				sum++;
				add_sample.push_back(-1);
			}else
				cout << "gap error : " << gap << endl;
		}
	}

	total=(long)sum;

	return total;

}


void fix_time(double *time, double *&time_fixed, std::vector <long> indice, std::vector <long> &add_sample, double fsamp, long nsamples_total){


	long pointer_time=0;
	long jj=0, kk=0, uu=0;

	for(long ii=0; ii < (long)indice.size(); ii++){
		for(jj=pointer_time; jj<=indice[ii]; jj++){
			time_fixed[jj]=time[uu];
			uu++;
		}
		pointer_time=jj;
		if(add_sample[ii]==-1){
			time_fixed[pointer_time]=(time[uu+1]+time[uu])/2;
			add_sample[ii]=1;
		}else{
			for(kk=pointer_time; kk < pointer_time + add_sample[ii]; kk++)
				time_fixed[kk]=time_fixed[kk-1]+1/fsamp;
			pointer_time=kk;
		}
		pointer_time++;
	}
	//	cout << "pt : " << pointer_time << endl;

	if((nsamples_total-pointer_time)>0)
		for(kk=pointer_time; kk <nsamples_total; kk++){
			time_fixed[kk]=time[uu];
			uu++;
		}
}

void fix_row(double *row, double *&row_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total){

	long pointer=0;
	long jj=0, kk=0, uu=0;

	for(long ii=0; ii < (long)indice.size(); ii++){
		for(jj=pointer; jj<=indice[ii]; jj++){
			row_fixed[jj]=row[uu];
			uu++;
		}
		pointer=jj;
		for(kk=pointer; kk < pointer + add_sample[ii]; kk++)
			row_fixed[kk]=NAN;
		pointer=kk;

	}

	if((nsamples_total-pointer)>0)
		for(kk=pointer; kk <nsamples_total; kk++){
			row_fixed[kk]=row[uu];
			uu++;
		}

}

void fix_mask(int *mask, int *&mask_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total){

	long pointer=0;
	long jj=0, kk=0, uu=0;

	for(long ii=0; ii < (long)indice.size(); ii++){
		for(jj=pointer; jj<=indice[ii]; jj++){
			mask_fixed[jj]=mask[uu];
			uu++;
		}
		pointer=jj;
		for(kk=pointer; kk < pointer + add_sample[ii]; kk++)
			mask_fixed[kk]=1;
		pointer=kk;

	}

	if((nsamples_total-pointer)>0)
		for(kk=pointer; kk <nsamples_total; kk++){
			mask_fixed[kk]=mask[uu];
			uu++;
		}

}
void insert_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final){

	int status=0;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time, &status);


}


void insert_row_in_image(fitsfile *fptr, fitsfile *outfptr, string field, double *row_fixed, long ns_total){


	int status=0;


	long fpixel[2]={1,1};
	long rowIndex = find_channel_index(fptr, field.c_str());
	fpixel[1] = rowIndex;
	fits_write_pix(outfptr, TDOUBLE, fpixel, ns_total, row_fixed, &status);

}

void insert_mask_in_image(fitsfile *fptr, fitsfile *outfptr, string field, int *mask_fixed, long ns_total){

	int status=0;


	long fpixel[2]={1,1};
	long rowIndex = find_channel_index(fptr, field.c_str());
	fpixel[1] = rowIndex;
	fits_write_pix(outfptr, TINT, fpixel, ns_total, mask_fixed, &status);
}

void insert_ref_pos_in_fits(fitsfile *fptr, fitsfile *outfptr, double *RA, double *DEC,double *PHI, long ns_total){

	int status = 0;

	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_total, RA, &status);
	fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_total, DEC, &status);
	fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_total, PHI, &status);
}

// TODO : same as in saneSplit, should teleport to a common library ?
void copy_offsets(fitsfile * fptr, fitsfile *outfptr){

	int status=0;

	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);

	for(int col=1;col<4;col++)
		fits_copy_col(fptr, outfptr,  col, col,	0, &status);

}

// TODO : same as in saneSplit, should teleport to a common library ?
void copy_channels(fitsfile * fptr, fitsfile *outfptr){

	int status=0;


	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);


	fits_copy_col(fptr, outfptr,  1, 1,	0, &status);

}

void fix_signal(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample ){

	long ns_temp = 0;
	int status =0;
	double *signal, *signal_fixed;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status);

	signal_fixed = new double [ns_total];

	for(long jj=0;jj<det.ndet;jj++){
		read_signal_from_fits(name, det.boloname[jj], signal, ns_temp);
		fix_row(signal, signal_fixed, indice, add_sample, ns_total);
		insert_row_in_image(fptr, outfptr, det.boloname[jj], signal_fixed, ns_total);
		delete [] signal;
	}
	delete [] signal_fixed;

}


//void fix_RA(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample ){
//
//	long ns_temp = 0;
//	int status =0;
//	double *RA, *RA_fixed;
//
//
//	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status);
//	fits_copy_header(fptr, outfptr, &status);
//	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status);
//
//	RA_fixed = new double [ns_total];
//
//	for(long jj=0;jj<det.ndet;jj++){
//		read_ra_from_fits(name, det.boloname[jj], RA, ns_temp);
//		fix_row(RA, RA_fixed, indice, add_sample, ns_total);
//		insert_row_in_image(fptr, outfptr, det.boloname[jj], RA_fixed, ns_total);
//		//				cout << RA_fixed [167294 + 2 ] << " " << RA_fixed [167294 + 10 ] << " " << RA_fixed [167294 + 1000 ];
//		//				cout << endl;
//		//				getchar();
//		delete [] RA;
//	}
//	delete [] RA_fixed;
//}

void fix_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample){

	long ns_temp = 0;
	int status =0;
	double *RA, *RA_fixed;
	double *DEC, *DEC_fixed;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status);
	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status);

	RA_fixed = new double [ns_total];
	DEC_fixed = new double [ns_total];

	for(long jj=0;jj<det.ndet;jj++){
		read_ra_dec_from_fits(name, det.boloname[jj], RA, DEC, ns_temp);
		for(long nn=0; nn<ns_temp;nn++)
			RA[nn]=RA[nn]*15.0;
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "ra", NULL, &status);
		fix_row(RA, RA_fixed, indice, add_sample, ns_total);
		insert_row_in_image(fptr, outfptr, det.boloname[jj], RA_fixed, ns_total);

		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "dec", NULL, &status);
		fix_row(DEC, DEC_fixed, indice, add_sample, ns_total);
		insert_row_in_image(fptr, outfptr, det.boloname[jj], DEC_fixed, ns_total);
		delete [] RA;
		delete [] DEC;
	}
	delete [] RA_fixed;
	delete [] DEC_fixed;

}

//void fix_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample ){
//
//	long ns_temp = 0;
//	int status =0;
//	double *DEC, *DEC_fixed;
//
//	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status);
//	fits_copy_header(fptr, outfptr, &status);
//	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status);
//
//	DEC_fixed = new double [ns_total];
//
//	for(long jj=0;jj<det.ndet;jj++){
//		read_dec_from_fits(name, det.boloname[jj], DEC, ns_temp);
//		fix_row(DEC, DEC_fixed, indice, add_sample, ns_total);
//		insert_row_in_image(fptr, outfptr, det.boloname[jj], DEC_fixed, ns_total);
//		delete [] DEC;
//	}
//	delete [] DEC_fixed;
//
//
//}

void fix_mask(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample ){

	long ns_temp = 0;
	int status =0;
	int *mask, *mask_fixed;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status);

	mask_fixed = new int [ns_total];

	for(long jj=0;jj<det.ndet;jj++){
		read_flag_from_fits(name, det.boloname[jj], mask, ns_temp);
		fix_mask(mask, mask_fixed, indice, add_sample, ns_total);
		insert_mask_in_image(fptr, outfptr, det.boloname[jj], mask_fixed, ns_total);
		delete [] mask;
	}
	delete [] mask_fixed;

}

void fix_time_table(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, long ns_origin, double fsamp){

	double *time;
	double *time_fixed;
	time_fixed= new double[ns_total];
	read_time_from_fits(name, time, ns_origin);
	fix_time(time,time_fixed,indice, add_sample, fsamp, ns_total);
	insert_time(fptr, outfptr, time_fixed, ns_total);
	delete [] time;
	delete [] time_fixed;

}


void fix_ref_pos(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample ){

	long ns_temp = 0;
	int status =0;
	double *RA, *RA_fixed;
	double *DEC, *DEC_fixed;
	double *PHI, *PHI_fixed;

	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_total, (char*)"Number of rows", &status);

	read_ReferencePosition_from_pointer(fptr, RA, DEC, PHI, ns_temp);
	for(long nn=0; nn<ns_temp;nn++)
		RA[nn]=RA[nn]*15.0;

	RA_fixed = new double [ns_total];
	DEC_fixed = new double [ns_total];
	PHI_fixed = new double [ns_total];

	fix_row(RA, RA_fixed, indice, add_sample, ns_total);
	fix_row(DEC, DEC_fixed, indice, add_sample, ns_total);
	fix_row(PHI, PHI_fixed, indice, add_sample, ns_total);

	insert_ref_pos_in_fits(fptr, outfptr, RA_fixed, DEC_fixed, PHI_fixed, ns_total);

	delete [] RA;
	delete [] DEC;
	delete [] PHI;
	delete [] RA_fixed ;
	delete [] DEC_fixed ;
	delete [] PHI_fixed ;

}
