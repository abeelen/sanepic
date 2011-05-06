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
#include <fitsio.h>
}


using namespace std;


int read_indices_file(string fname, struct param_common dir, std::vector<long> &indice, double &fsamp, long &init_num_delete, long &end_num_delete)
/*! read saneCheck log files : get sample indices => where the gaps are */
{

	std::ifstream file;
	long readed;

	string fname2 = dir.tmp_dir + FitsBasename(fname) + "_saneFix_indices.bin"; // sanecheck log filename

	file.open(fname2.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << fname << "] was not found." << endl;
		file.close();
		return(EXIT_FAILURE);
	}

	// getting initial number of flagged data to remove
	file >> init_num_delete;

	// getting final number of flagged data to remove
	file >> end_num_delete;

	file >> fsamp; // read sampling frequency computed by sanecheck
	while(file >> readed)
		indice.push_back(readed); // read the whole list of gaps indexes

	file.close();
	return(EXIT_SUCCESS);
}

void refresh_indice(double &fsamp, long init_num_delete, long end_num_delete, std::vector <long> &indice, long ns){


	std::vector<long>::iterator it;

	long ni = (long)indice.size();
	if(ni>0){
		it=indice.begin();

		for(long ii=0; ii<ni; ii++){
			if(indice[ii]<=init_num_delete){
				//				cout << "erasing : ii = " << ii << " it = " << *(it) << endl;
				indice.erase(it);
			}else
				continue;
		}
	}

	ni = (long)indice.size();
	if(ni>0){
		it=indice.end();
		it--;

		while(*(it)>=ns-end_num_delete){
			indice.pop_back();
			if((int)indice.size()==0)
				return;
			it--;
		}
	}
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

	if((long)indice.size()>0){
		for(long ii=0; ii < (long)indice.size(); ii++){ // for each gap
			if(continue_neg>0){
				continue_neg--;
				continue;
			}
#ifdef DEBUG
			cout << indice[ii] << " " << time[indice[ii]+1]-time[indice[ii]] << " " << round((time[indice[ii]+1]-time[indice[ii]])*fsamp)-1 << endl;
#endif
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
//					cout << "negativ gap : " << gap << endl;
					continue_neg++;
					long jj=1;
					while((time[indice[ii]+jj]-time[indice[ii]])<0.0)
						jj++;
					gap2=round((time[indice[ii]+jj]-time[indice[ii]])*fsamp)-1;
//					cout << "gap2 : " << gap2 << endl;
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

#ifdef DEBUG
		cout << "result ! \n";
		for(long ii=0; ii < (long)add_sample.size(); ii++){
			cout << indice_valid[ii] << " " << add_sample[ii] << " " << suppress_time_sample[ii] << endl;
		}
#endif
		indice.clear();
		indice=indice_valid;

		total=(long)sum; // total number of samples that must be added

	}else
		total=0;

	delete [] time;

	return total;

}


void fix_time(double *time, double *&time_fixed, std::vector <long> indice, std::vector <long> &add_sample, double fsamp, long nsamples_total, std::vector <long> suppress_time_sample, long init_num_delete)
/*! fix time table : draw a continuous time table by filling gaps using the calculated sampling frequency */
{


	long pointer_time=0;
	long jj=0, kk=0, uu=init_num_delete;

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

void fix_row(double *row, double *&row_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total, std::vector <long> suppress_time_sample, long init_num_delete)
/*! fix a row vector : given a row, fill this row with average values */
{

	long pointer=0;
	long jj=0, kk=0, uu=init_num_delete;

	for(long ii=0; ii < (long)indice.size(); ii++){ // for each gap
		jj=pointer;
		while(uu<=indice[ii]){ // copy row untill a gap is found
			row_fixed[jj]=row[uu];
			uu++;
			jj++;
		}
		pointer=jj;
		for(kk=pointer; kk < pointer + add_sample[ii]; kk++) // fill with average values
			row_fixed[kk] = row[uu-1]+(row[uu+1+suppress_time_sample[ii]] - row[uu-1])/add_sample[ii]*(kk+1-pointer);
		pointer=kk;
		uu+=suppress_time_sample[ii];

	}

	if((nsamples_total-pointer)>0) // copy row vector untill its end
		for(kk=pointer; kk <nsamples_total; kk++){
			row_fixed[kk]=row[uu];
			uu++;
		}

}

void fix_mask(int *mask, int *&mask_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total, std::vector <long> suppress_time_sample, long init_num_delete)
/*! Fill a mask row vector with ones for each gap generated values */
{


	for(long s=0; s<(long)add_sample.size(); s++)
		if(add_sample[s]==-1)
			add_sample[s]=1;

	long ind = 0;
	long p_copy = init_num_delete;
	long ii=0;

	while(ii<nsamples_total){
		if(ind>=(long)indice.size()){
			mask_fixed[ii]=mask[p_copy];
			p_copy ++;
			ii++;
		}else{

			mask_fixed[ii]=mask[p_copy]; // p_copy = 0 cannot be a wrong data ! (convention), it has been removed before if time was bad
			p_copy ++;
			ii++;

			if(p_copy==indice[ind]){
				mask_fixed[ii]=mask[p_copy];
				p_copy ++;
				ii++;
				for(long jj=0; jj< add_sample[ind]; jj++)
					mask_fixed[ii+jj]=1;

				ii=ii+add_sample[ind];
				p_copy+=suppress_time_sample[ind];
				ind++;

			}
		}
	}







	//	long ind = 0;
	//	long p_copy = 0;
	//	long ii=0;
	//
	//	while(ii<nsamples_total){
	//		if(ind>=(long)indice.size()){
	//			mask_fixed[ii]=mask[p_copy];
	//			p_copy ++;
	//		}else{
	//			mask_fixed[ii]=mask[p_copy];
	//			p_copy ++;
	//
	//			if(ii==indice[ind]){
	////				cout << "ii : " << ii << " addsample : " << add_sample[ind] << " ind : " << ind << endl;
	////				getchar();
	//				for(long jj=1; jj<= add_sample[ind]; jj++)
	//					mask_fixed[ii+jj]=1;
	//
	//				ii=ii+add_sample[ind];
	//				ind++;
	//				//			}else{
	//				//				mask_fixed[ii]=mask[p_copy];
	//				//								p_copy ++;
	//				//			}
	//			}
	//		}
	////		cout << "ii : " << ii << endl;
	////		cout << "p_copy : " << p_copy << endl;
	//		ii++;
	//	}

}
void insert_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final)
/*! Copy time header from original fits file to output, insert final time table and update header */
{

	int status=0; // fits error status

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", 0, &status); // move pointer to time HDU
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

void fix_signal(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, std::vector<std::string> det, long ndet, std::vector <long> indice,
		std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete)
/*! Copy input signal header to output and fill the gaps in signal table */
{

	long ns_temp = 0;
	int status =0; // fits error status
	double *signal, *signal_fixed;


	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", 0, &status); // move input pointer to "signal" image
	fits_copy_header(fptr, outfptr, &status); // copy header to output
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output header

	signal_fixed = new double [ns_total];

	for(long jj=0;jj<ndet;jj++){ // for each detector (column)
		read_signal_from_fits(name, det[jj], signal, ns_temp); // read signal row
		fix_row(signal, signal_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete); // fill the gaps
		insert_row_in_image(fptr, outfptr, det[jj], signal_fixed, ns_total); // insert in output fits file
		delete [] signal;
	}
	delete [] signal_fixed;

}

void fix_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, std::vector<std::string> det, long ndet,
		std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete)
/*! Copy input RA and DEC header to output and fill the gaps in those tables : HIPE format only */
{

	long ns_temp = 0;
	int status =0; // fits error status
	double *RA, *RA_fixed;
	double *DEC, *DEC_fixed;

	if(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", 0, &status)){ // move input pointer to RA
		cout << "WARNING : ra table was not found, skipping this table...\n";
		return;
	}
	fits_copy_header(fptr, outfptr, &status); // copy RA header to output file
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output RA header
	if(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", 0, &status)){ // move input pointer to DEC

		cout << "WARNING : ra table was not found, skipping this table...\n";
		return;
	}

	fits_copy_header(fptr, outfptr, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output DEC header

	RA_fixed = new double [ns_total];
	DEC_fixed = new double [ns_total];

	for(long jj=0;jj<ndet;jj++){ // for each detector (column)
		read_ra_dec_from_fits(name, det[jj], RA, DEC, ns_temp); // read input RA and DEC row
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "ra", 0, &status); // move output pointer to RA table
		fix_row(RA, RA_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete); // fill gaps in RA row
		insert_row_in_image(fptr, outfptr, det[jj], RA_fixed, ns_total); // insert the filled RA row in ouput table


		// same process for DEC
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "dec", 0, &status);
		fix_row(DEC, DEC_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete);
		insert_row_in_image(fptr, outfptr, det[jj], DEC_fixed, ns_total);
		delete [] RA;
		delete [] DEC;
	}
	delete [] RA_fixed;
	delete [] DEC_fixed;

}

void fix_mask(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, std::vector<std::string> det, long ndet,
		std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete)
/*! Copy input mask header to output and fill the gaps with ones */
{

	long ns_temp = 0;
	int status =0; // fits error status
	int *mask, *mask_fixed;

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", 0, &status); // move input pointer to mask
	fits_copy_header(fptr, outfptr, &status); // copy header to ouput
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_total, (char*)"Number of rows", &status); // update output mask header

	mask_fixed = new int [ns_total];
	for(long ii=0;ii<ns_total;ii++)
		mask_fixed[ii]=-1.0;

	for(long jj=0;jj<ndet;jj++){ // for each detector (column)
		read_flag_from_fits(name, det[jj], mask, ns_temp); // read input mask row
		fix_mask(mask, mask_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete); // fill gaps in mask row
		insert_mask_in_image(fptr, outfptr, det[jj], mask_fixed, ns_total); // insert the filled mask row in ouput table
		delete [] mask;
	}
	delete [] mask_fixed;

}

void fix_time_table(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, std::vector <long> indice, std::vector<long> add_sample, long ns_origin,
		double fsamp, std::vector <long> suppress_time_sample, long init_num_delete)
/*! Copy input time header to output and fill the gaps with computed values using sampling frequency */
{

	double *time;
	double *time_fixed;
	time_fixed= new double[ns_total];
	read_time_from_fits(name, time, ns_origin); // read input time table
	fix_time(time,time_fixed,indice, add_sample, fsamp, ns_total, suppress_time_sample, init_num_delete); // fill gaps in time table
	insert_time(fptr, outfptr, time_fixed, ns_total); // insert table in output fits file
	delete [] time;
	delete [] time_fixed;

}


void fix_ref_pos(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, std::vector <long> indice, std::vector<long> add_sample,
		std::vector <long> suppress_time_sample, long init_num_delete)
/*! copy "reference position" table to output and fill the gaps with average values */
{

	long ns_temp = 0;
	int status =0; // fits error status
	double *RA, *RA_fixed;
	double *DEC, *DEC_fixed;
	double *PHI, *PHI_fixed;

	if(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", 0, &status)){ // move input pointer to ref pos
		cout << "WARNING : reference position table was not found, skipping this table...\n";
		return;
	}

	fits_copy_header(fptr, outfptr, &status); // copy header to output
	fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_total, (char*)"Number of rows", &status); // update output header (sample size has changed)

	read_ReferencePosition_from_fits(name, RA, DEC, PHI, ns_temp);


	for(long nn=0; nn<ns_temp;nn++)
		RA[nn]=RA[nn]*15.0;

	RA_fixed = new double [ns_total];
	DEC_fixed = new double [ns_total];
	PHI_fixed = new double [ns_total];

	// fill gaps in each table
	fix_row(RA, RA_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete);
	fix_row(DEC, DEC_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete);
	fix_row(PHI, PHI_fixed, indice, add_sample, ns_total, suppress_time_sample, init_num_delete);

	insert_ref_pos_in_fits(fptr, outfptr, RA_fixed, DEC_fixed, PHI_fixed, ns_total); // insert the 3 tables in output file

	delete [] RA;
	delete [] DEC;
	delete [] PHI;
	delete [] RA_fixed ;
	delete [] DEC_fixed ;
	delete [] PHI_fixed ;

}
