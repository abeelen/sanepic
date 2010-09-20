

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



int who_do_it(int size, int rank, int ii)
/*! this function determines which processor has to treat the given fits file referenced by his number in the input list */
{

	if(size==1) // if there is only 1 proc, he has to do the job
		return 0;

	if(size>=ii) // if the fits file number is smaller than the number of MPI processors
		return ii;

	if(size<ii){ // if the fits file number is larger than the number of MPI processors
		while(ii>size)
			ii=ii-size; // find the processor that will do the job by substracting iteratively the number of MPI procs
		return ii;
	}

	return -1; // error in the program
}


int read_indices_file(string fname, struct common dir, std::vector<long> &indice, double &fsamp)
/*! read saneCheck log files : get sample indices => where the gaps are */
{

	std::ifstream file;
	long readed;

	std::ostringstream oss;
	oss << fname;
	string filename = oss.str();
	string fname2 = dir.tmp_dir + Basename(filename) + "_saneFix_indices.bin"; // sanecheck log filename

	file.open(fname2.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << fname << "] was not found." << endl;
		file.close();
		return(EXIT_FAILURE);
	}

	file >> fsamp; // read sampling frequency computed by sanecheck
	while(file >> readed)
		indice.push_back(readed); // read the whole list of gaps indexes

	file.close();
	return(EXIT_SUCCESS);
}

long how_many(string fname, long ns, std::vector <long> &indice, double fsamp,  std::vector <long> &add_sample, std::vector <long> & suppress_time_sample)
/*!  compute the number of sample that must be added to fill the gaps and have a continous timeline  */
{

	long total=0; // total number of sample that will be added to input file
	long gap=0, sum=0, gap2=0; // if a gap has been found : gap > 0, if gap = 0 => one sample must be added
	int continue_neg=0;
	double *time;
	std::vector <long> indice_valid;

	read_time_from_fits(fname, time, ns); // read time table from input fits

	for(long ii=0; ii < (long)indice.size(); ii++){ // for each gap
		if(continue_neg>0){
			continue_neg--;
			continue;
		}
		cout << indice[ii] << " " << time[indice[ii]+1]-time[indice[ii]] << " " << round((time[indice[ii]+1]-time[indice[ii]])*fsamp)-1 << endl;
		//  in case gap negatif : sauter à l'indice suivant et remplacer le time faux !
		// calculate the gap size
		gap=round((time[indice[ii]+1]-time[indice[ii]])*fsamp)-1;
		if(gap>0){
			sum+=gap;
			add_sample.push_back(gap); // store the number of samples that must be added to fill this gap
			indice_valid.push_back(indice[ii]);
			suppress_time_sample.push_back(0);
		}else{
			if(gap==0){
				sum++;
				add_sample.push_back(-1); // store -1 means 1 sample must be added => it's a convention
				indice_valid.push_back(indice[ii]);
				suppress_time_sample.push_back(0);
			}else{
				cout << "gap negatif : " << gap << endl;
				continue_neg++;
				long jj=1;
				while((time[indice[ii]+jj]-time[indice[ii]])<0.0)
					jj++;
				gap2=round((time[indice[ii]+jj]-time[indice[ii]])*fsamp)-1;
				cout << "gap2 : " << gap2 << endl;
				if(gap==0){
					sum++;
					add_sample.push_back(-1); // store -1 means 1 sample must be added => it's a convention
					indice_valid.push_back(indice[ii]);
					suppress_time_sample.push_back(jj-1);
				}else{
					sum+=gap2-jj+1;
					add_sample.push_back(gap2);
					indice_valid.push_back(indice[ii]);
					suppress_time_sample.push_back(jj-1);
				}
			}
		}
	}

	cout << "result ! \n";
	for(long ii=0; ii < (long)add_sample.size(); ii++){
		cout << indice_valid[ii] << " " << add_sample[ii] << endl;
	}

	indice=indice_valid;

	total=(long)sum; // total number of samples that must be added

	delete [] time;

	return total;

}


void fix_time(double *time, double *&time_fixed, std::vector <long> indice, std::vector <long> &add_sample, double fsamp, long nsamples_total, std::vector <long> suppress_time_sample)
/*! fix time table : draw a continuous time table by filling gaps using the calculated sampling frequency */
{


	long pointer_time=0;
	long jj=0, kk=0, uu=0;

	for(long ii=0; ii < (long)indice.size(); ii++){ // for each gap
		jj=pointer_time;
		while(uu<=indice[ii]){ // copy time vector to output table untill a gap is found
			time_fixed[jj]=time[uu];
			uu++;
			jj++;
		}
		pointer_time=jj;
		if(add_sample[ii]==-1){ // compute the mean only if there is 1 sample to add
			time_fixed[pointer_time]=(time[uu+1]+time[uu])/2;
			add_sample[ii]=1;
			//			cout << "detected -1\n";
		}else{
			uu+=suppress_time_sample[ii];
			for(kk=pointer_time; kk < pointer_time + add_sample[ii]; kk++){ // fil the gap using sampling frequency
				time_fixed[kk]=time_fixed[pointer_time-1]+(kk-pointer_time+1)/fsamp;
			}
			pointer_time=kk;
		}
	}

	if((nsamples_total-pointer_time)>0) // copy time vector to the end
		for(kk=pointer_time; kk <nsamples_total; kk++){
			time_fixed[kk]=time[uu];
			uu++;
		}
}

void fix_row(double *row, double *&row_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total, std::vector <long> suppress_time_sample)
/*! fix a row vector : given a row, fill this row with average values */
{

	long pointer=0;
	long jj=0, kk=0, uu=0;

	for(long ii=0; ii < (long)indice.size(); ii++){ // for each gap
		jj=pointer;
		while(uu<=indice[ii]){ // copy row untill a gap is found
			row_fixed[jj]=row[uu];
			uu++;
			jj++;
		}
		pointer=jj;
		for(kk=pointer; kk < pointer + add_sample[ii]; kk++) // fill with average values
			row_fixed[kk] = row[uu]+(row[uu+1+suppress_time_sample[ii]] - row[uu])/add_sample[ii]*(kk-pointer);
		pointer=kk;
		uu+=suppress_time_sample[ii];

	}

	if((nsamples_total-pointer)>0) // copy row vector untill its end
		for(kk=pointer; kk <nsamples_total; kk++){
			row_fixed[kk]=row[uu];
			uu++;
		}

}

void fix_mask(int *mask, int *&mask_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total, std::vector <long> suppress_time_sample)
/*! Fill a mask row vector with ones for each gap generated values */
{

	long pointer=0;
	long jj=0, kk=0, uu=0;

	for(long ii=0; ii < (long)indice.size(); ii++){ // for each gap
		for(jj=pointer; jj<=indice[ii]; jj++){ // copy mask vector untill gap
			mask_fixed[jj]=mask[uu];
			uu++;
		}
		pointer=jj;
		for(kk=pointer; kk < pointer + add_sample[ii]; kk++) // fill gaps values with ones
			mask_fixed[kk]=1;
		pointer=kk;
		uu+=suppress_time_sample[ii];

	}

	if((nsamples_total-pointer)>0) // copy vector untill its end
		for(kk=pointer; kk <nsamples_total; kk++){
			mask_fixed[kk]=mask[uu];
			uu++;
		}

}
void insert_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final)
/*! Copy time header from original fits file to output, insert final time table and update header */
{

	int status=0; // fits error status

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status); // move pointer to time HDU
	fits_copy_header(fptr, outfptr, &status); // copy header to ouput
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status); // update sample number

	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time, &status); // write time table in output file


}


void insert_row_in_image(fitsfile *fptr, fitsfile *outfptr, string field, double *row_fixed, long ns_total)
/*! Insert one fixed row in the output fits file */
{
	// NOTE : You must place the output fits pointer to the correct HDU before this !

	int status=0; // fits error status


	long fpixel[2]={1,1};
	long rowIndex = find_channel_index(fptr, field.c_str()); // find row index corresponding to channel name
	fpixel[1] = rowIndex;
	fits_write_pix(outfptr, TDOUBLE, fpixel, ns_total, row_fixed, &status); // write row in the output fits

}

void insert_mask_in_image(fitsfile *fptr, fitsfile *outfptr, string field, int *mask_fixed, long ns_total)
/*! insert a mask row in the output mask image */
{

	int status=0; // fits error status


	long fpixel[2]={1,1};
	long rowIndex = find_channel_index(fptr, field.c_str());
	fpixel[1] = rowIndex;
	fits_write_pix(outfptr, TINT, fpixel, ns_total, mask_fixed, &status);
}

void insert_ref_pos_in_fits(fitsfile *fptr, fitsfile *outfptr, double *RA, double *DEC,double *PHI, long ns_total)
/*! insert "reference position" table in the output fits file */
{

	int status = 0; // fits error status

	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_total, RA, &status); // write columns one by one
	fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_total, DEC, &status);
	fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_total, PHI, &status);
}


void copy_offsets(fitsfile * fptr, fitsfile *outfptr)
/*! copy offsets table from input to output file */
{

	int status=0; // fits error status

	if(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status)){// move input pointer to "offsets" HDU
		cout << "WARNING : offsets table was not found, skipping this table...\n";
		return;
	}
	fits_copy_header(fptr, outfptr, &status); // copy header to output file

	for(int col=1;col<4;col++) // copy the table column by column
		fits_copy_col(fptr, outfptr,  col, col,	0, &status);

}


void copy_channels(fitsfile * fptr, fitsfile *outfptr)
/*! copy offsets table from input to output file */
{

	int status=0; // fits error status


	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status); // move input pointer to "channels" table
	fits_copy_header(fptr, outfptr, &status); // copy header to output
	fits_copy_col(fptr, outfptr,  1, 1,	0, &status); // copy whole channel table to output

}

void fix_signal(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample)
/*! Copy input signal header to output and fill the gaps in signal table */
{

	long ns_temp = 0;
	int status =0; // fits error status
	double *signal, *signal_fixed;


	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status); // move input pointer to "signal" image
	fits_copy_header(fptr, outfptr, &status); // copy header to output
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output header

	signal_fixed = new double [ns_total];

	for(long jj=0;jj<det.ndet;jj++){ // for each detector (column)
		read_signal_from_fits(name, det.boloname[jj], signal, ns_temp); // read signal row
		fix_row(signal, signal_fixed, indice, add_sample, ns_total, suppress_time_sample); // fill the gaps
		insert_row_in_image(fptr, outfptr, det.boloname[jj], signal_fixed, ns_total); // insert in output fits file
		delete [] signal;
	}
	delete [] signal_fixed;

}

void fix_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample)
/*! Copy input RA and DEC header to output and fill the gaps in those tables : HIPE format only */
{

	long ns_temp = 0;
	int status =0; // fits error status
	double *RA, *RA_fixed;
	double *DEC, *DEC_fixed;

	if(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status)){ // move input pointer to RA
		cout << "WARNING : ra table was not found, skipping this table...\n";
		return;
	}
	fits_copy_header(fptr, outfptr, &status); // copy RA header to output file
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output RA header
	if(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status)){ // move input pointer to DEC

		cout << "WARNING : ra table was not found, skipping this table...\n";
		return;
	}

	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output DEC header

	RA_fixed = new double [ns_total];
	DEC_fixed = new double [ns_total];

	for(long jj=0;jj<det.ndet;jj++){ // for each detector (column)
		read_ra_dec_from_fits(name, det.boloname[jj], RA, DEC, ns_temp); // read input RA and DEC row
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "ra", NULL, &status); // move output pointer to RA table
		fix_row(RA, RA_fixed, indice, add_sample, ns_total, suppress_time_sample); // fill gaps in RA row
		insert_row_in_image(fptr, outfptr, det.boloname[jj], RA_fixed, ns_total); // insert the filled RA row in ouput table


		// same process for DEC
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "dec", NULL, &status);
		fix_row(DEC, DEC_fixed, indice, add_sample, ns_total, suppress_time_sample);
		insert_row_in_image(fptr, outfptr, det.boloname[jj], DEC_fixed, ns_total);
		delete [] RA;
		delete [] DEC;
	}
	delete [] RA_fixed;
	delete [] DEC_fixed;

}

void fix_mask(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample)
/*! Copy input mask header to output and fill the gaps with ones */
{

	long ns_temp = 0;
	int status =0; // fits error status
	int *mask, *mask_fixed;
//	long ii=1;// singleton seeker indexes

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status); // move input pointer to mask
	fits_copy_header(fptr, outfptr, &status); // copy header to ouput
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output mask header

	mask_fixed = new int [ns_total];

	for(long jj=0;jj<det.ndet;jj++){ // for each detector (column)
		read_flag_from_fits(name, det.boloname[jj], mask, ns_temp); // read input mask row
		fix_mask(mask, mask_fixed, indice, add_sample, ns_total, suppress_time_sample); // fill gaps in mask row
//		ii=1;
//		while(ii<ns_total-1){
//			if((mask_fixed[ii]==0)&&(mask_fixed[ii+1]!=0)&&(mask_fixed[ii-1]!=0)){
////				mask_fixed[ii]=1; // TODO : do we keep singletons now ??
//				cout << "singleton found : " << det.boloname[jj] << " sample n° " << ii << endl;
//			}
//			ii++;
//		}
		insert_mask_in_image(fptr, outfptr, det.boloname[jj], mask_fixed, ns_total); // insert the filled mask row in ouput table
		delete [] mask;
	}
	delete [] mask_fixed;

}

void fix_time_table(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, long ns_origin, double fsamp, std::vector <long> suppress_time_sample)
/*! Copy input time header to output and fill the gaps with computed values using sampling frequency */
{

	double *time;
	double *time_fixed;
	time_fixed= new double[ns_total];
	read_time_from_fits(name, time, ns_origin); // read input time table
	fix_time(time,time_fixed,indice, add_sample, fsamp, ns_total, suppress_time_sample); // fill gaps in time table
	insert_time(fptr, outfptr, time_fixed, ns_total); // insert table in output fits file
	delete [] time;
	delete [] time_fixed;

}


void fix_ref_pos(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample)
/*! copy "reference position" table to output and fill the gaps with average values */
{

	long ns_temp = 0;
	int status =0; // fits error status
	double *RA, *RA_fixed;
	double *DEC, *DEC_fixed;
	double *PHI, *PHI_fixed;

	if(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status)){ // move input pointer to ref pos
		cout << "WARNING : reference position table was not found, skipping this table...\n";
		return;
	}

	fits_copy_header(fptr, outfptr, &status); // copy header to output
	fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_total, (char*)"Number of rows", &status); // update output header (sample size has changed)

	read_ReferencePosition_from_pointer(fptr, RA, DEC, PHI, ns_temp); // read input RA, DEC and PHI tables
	for(long nn=0; nn<ns_temp;nn++)
		RA[nn]=RA[nn]*15.0;

	RA_fixed = new double [ns_total];
	DEC_fixed = new double [ns_total];
	PHI_fixed = new double [ns_total];

	// fill gaps in each table
	fix_row(RA, RA_fixed, indice, add_sample, ns_total, suppress_time_sample);
	fix_row(DEC, DEC_fixed, indice, add_sample, ns_total, suppress_time_sample);
	fix_row(PHI, PHI_fixed, indice, add_sample, ns_total, suppress_time_sample);

	insert_ref_pos_in_fits(fptr, outfptr, RA_fixed, DEC_fixed, PHI_fixed, ns_total); // insert the 3 tables in output file

	delete [] RA;
	delete [] DEC;
	delete [] PHI;
	delete [] RA_fixed ;
	delete [] DEC_fixed ;
	delete [] PHI_fixed ;

}
