
#include "covMatrixIO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneSplit.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}




using namespace std;

//class vector2D_long
//{
//public :
//	std::vector<long> sample_limits;
//
//};


void usage(char *name)
{
	cerr << "USAGE: " << name << " inifile.ini [-f <path/filename>] [-m <min time>] [-M <max time>]" << endl;

}


int main(int argc, char *argv[]) {

	int parsed = -1;

	struct samples samples_struct;
	struct directories dir;
	string fname;
	std::vector< double > min_time, max_time;
	int retval;
	int m_count = 0, f_count = 0;
	int mM_count = 0;
	int format_fits=0;
	//	int *format=NULL;
	double *time=NULL;
	double time_min=0, time_max=0;
	struct detectors det;
	//	long nsamples_global=0;


	//	prog_name=(string)argv[0];



	printf("\nBeginning of saneCheck:\n\n");



	if (argc<2)
		parsed=1;
	else
		// Parse ini file
		parsed=parse_saneSplit_ini_file(argv[1],dir);

	if (parsed>0){
		switch (parsed){

		case 1: printf("Please run %s using a *.ini file\n",argv[0]);
		usage(argv[0]);
		break;

		case 2 : printf("Wrong program options or argument. Exiting !\n");
		break;

		default :;
		}

		exit(EXIT_FAILURE);
	}


	// read -m and -M options
	while ( (retval = getopt(argc, argv, "f:m:M:")) != -1) {
		switch (retval) {
		case 'f':
			samples_struct.fitsvect.push_back(optarg);
			readFrames(samples_struct.fitsvect, samples_struct.nsamples);
			cout << "Scan      : " << samples_struct.fitsvect[0] << endl;
			cout << "Containing      : " << samples_struct.nsamples[0] << " samples. " << endl;
			f_count++;
			break;
		case 'm':
			min_time.push_back(atof(optarg));
			m_count++;
			mM_count++;
			break;

		case 'M':
			max_time.push_back(atof(optarg));
			mM_count--;
			break;

		default:;
		}
	}


	if((f_count>1)||(f_count==0)){
		cout << "\nError ! You must give 1 file AND only one ! Exiting\n";
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	// if no m or no M, return
	if((m_count==0)||(mM_count!=0)){
		cout << "\nError ! You must give at least 1 min and max value AND each -M must follow a -m option ! Exiting\n";
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	//	long size_total=0;
	//	for(int ii=0;ii<samples_struct.ntotscan;ii++)
	//		size_total +=samples_struct.nsamples[ii];

	//	if(min_time[0]<0){
	//		cout << "Warning : you must provide Positives cut limits ! Exiting\n";
	//		exit(EXIT_FAILURE);
	//	}


	// if m<0 or M<m, return
	for(int ii=0;ii<m_count;ii++){
		if((min_time[ii]<0)||(max_time[ii]-min_time[ii])<=0){
			cout << "Warning : you must provide a crescent order of Positives cut limits : M>m !\nExiting\n";
			exit(EXIT_FAILURE);
		}
	}


	cout << setprecision(20) << min_time[0] << " " << max_time[0] << endl;
	//	cout << min_time[1] << " " << max_time[1] << endl;

	//	time_min=new double[samples_struct.ntotscan];
	//	time_max=new double[samples_struct.ntotscan];

	// get total numbers of samples (for all the scans)
	//	for(int ii = 0;ii<samples_struct.ntotscan;ii++)
	//		nsamples_global += samples_struct.nsamples[ii];
	// allocate memory for the whole data
	//	time_global=new double[nsamples_global];


	// get min and max time for each fits file (and store the time sequences)
	//	long indice=0;
	//	for(int jj=0;jj<samples_struct.ntotscan;jj++){
	time=new double[samples_struct.nsamples[0]];
	read_time_from_fits(samples_struct.fitsvect[0], time, samples_struct.nsamples[0]);
	//		for(long kk=0;kk<samples_struct.nsamples[jj];kk++){
	//			time_global[indice]=time[kk];
	//			indice++;
	//		}
	time_min=time[0];
	time_max=time[samples_struct.nsamples[0]-1];
	//	delete [] time;
	//}

	// get the min of all mins and max of all maxs
	//double inf_time=*min_element(time_min,time_min+samples_struct.ntotscan);
	//double sup_time=*max_element(time_max,time_max+samples_struct.ntotscan);

	//	cout << "fits limits  : \n";
	//	cout << setprecision(40) << time_min << endl;
	//	cout << setprecision(40) << time_max << endl;
	//	cout << m_count << endl;
	//	cout << setprecision(40) << min_time[0] << endl;
	//	cout << setprecision(40) << max_time[0] << endl;

	// if given max value (-M) is larger than the max time (fits), set M to this max value
	//	for(int ii=0;ii<m_count;ii++){
	int ii=0;
	while(ii<m_count){
		int indic=0;
		if(max_time[ii]>time_max){
			max_time[ii]=time_max;
			indic++;
			//cout << "Warning : You are trying to reach a time sample that is larger than the maximum time sample !\nExiting\n";
			//exit(EXIT_FAILURE);
		}
		// if given min value (-M) is shorter than the min time (fits), set m to this min value
		if(min_time[ii]<time_min){
			min_time[ii]=time_min;
			indic++;
			//cout << "Warning : You are trying to reach a time sample that is smaller than the minimum time sample !\nExiting\n";
			//exit(EXIT_FAILURE);
		}
		if(indic==2){
			cout << "Warning : You are trying to split a file into the exact same file... \nSwitching to next cut limits...\n";
			max_time.erase (max_time.begin()+ii);
			min_time.erase (min_time.begin()+ii);
			m_count--;
			ii--;
		}
		ii++;
	}

	if(m_count==0){
		cout << "Warning : no more cut limits to apply ! Exiting ...\n";
		exit(EXIT_FAILURE);
	}





	//	read_Split_file(fname, cut_sample, samples_struct);

	//	cout << cut_sample[0] << " " << cut_sample[1] << endl;

	//format=new int[samples_struct.ntotscan];
	// select the fits format (HIPE or SANEPIC) and check all fits have the same format
	//for(int ii =0;ii<samples_struct.ntotscan;ii++){
	format_fits=test_format(samples_struct.fitsvect[0]);
	//	if(format[ii]!=format[0]){
	//		cout << "Error ! All fits files must have the same format : HIPE or sanepic!\nExiting...\n";
	//		exit(EXIT_FAILURE);
	//	}
	//
	//}


	//format_fits=format[0];

	read_bolo_list(samples_struct.fitsvect[0], det);

	//	double *image;
	//	long ns, ndet;
	//	read_image_2D_from_fits(samples_struct.fitsvect[0], image, "signal", ns, ndet);
	//	cout << "ns " << ns << endl;
	//	cout << ndet << endl;
	//
	//	cout << image[0] << " " << image[1] << endl;
	//	cout << image[10][1] << " " << image[0][0] << endl;

	int status;
	fitsfile *fptr;
	fitsfile *outfptr;
	//	char section[]={"*,*"};
	std::ostringstream oss;
	oss << samples_struct.fitsvect[0];
	string filename = oss.str();
	string fname2 = Basename(filename) + "_split_";

	oss.str("");

	fname=samples_struct.fitsvect[0];
	for(int ii=0; ii < m_count ; ii++){
		oss << "!" << dir.outdir << fname2 << setprecision(14) << min_time[ii] << "_" << max_time[ii] << ".fits";
		std::string temp = oss.str();
		cout << temp << endl;

		//		cout << fname << endl;
		//		getchar();
		//exit(0);

		long min_sample=0;
		long max_sample=0;

		// find samples index using time :
		if(max_time[ii]==time[samples_struct.nsamples[0]-1])
			max_sample=samples_struct.nsamples[0]-1;
		else{
			//for(long jj=0;jj<samples_struct.nsamples[0];jj++)
			long jj=0;
			while((jj<samples_struct.nsamples[0])){
				if(max_time[ii]<time[jj]){
					max_sample=jj;
					cout << max_time[ii] << " < " << time[jj] << endl;
					break;
				}
				jj++;
			}
		}

		if(min_time[ii]==time[0])
			min_sample=0;
		else{
			long jj=samples_struct.nsamples[0]-1;
			while(jj>=0){
				if(min_time[ii]>time[jj])
					min_sample=jj; // maybe jj+1
				jj--;
			}
		}

		cout << min_sample << " " << max_sample << endl;

		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);

		if (fits_create_file(&outfptr, temp.c_str(), &status))
			fits_report_error(stderr, status);

		//		char *header;
		//		int nkeys=0;

		fits_copy_header(fptr, outfptr, &status);

		double *RA, *DEC, *PHI;
		double *RA_bis, *DEC_bis, *PHI_bis, *time_bis;
		long ns_temp, ns_final;

		read_ReferencePosition_from_fits(samples_struct.fitsvect[0], RA, DEC, PHI, ns_temp);

		ns_final = max_sample - min_sample;
		ns_final++;
		RA_bis = new double [ns_final];
		DEC_bis = new double [ns_final];
		PHI_bis = new double [ns_final];
		time_bis = new double [ns_final];

		for(long ii = 0; ii< ns_final; ii++){
			RA_bis[ii]=RA[ii]*15.0;
			DEC_bis[ii]=DEC[ii];
			PHI_bis[ii]=PHI[ii];
			time_bis[ii]=time[ii];
		}



		fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status);
		fits_copy_header(fptr, outfptr, &status);


		// insert column
		//fits_insert_col(fptr, 1, TDOUBLE, , &status);
		fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, RA_bis, &status);
		fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_final, DEC_bis, &status);
		fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_final, PHI_bis, &status);
		fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_final, (char*)"Number of rows", &status);



		fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", NULL, &status);
		fits_copy_header(fptr, outfptr, &status);

		for(int col=1;col<4;col++)
			fits_copy_col(fptr, outfptr,  col, col,	0, &status);

		fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status);
		fits_copy_header(fptr, outfptr, &status);


		fits_copy_col(fptr, outfptr,  1, 1,	0, &status);

		fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status);
		fits_copy_header(fptr, outfptr, &status);

		fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status);
		fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);


		fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status);
		fits_copy_header(fptr, outfptr, &status);
		fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

		double *signal, *signal_bis;
		signal_bis = new double [ns_final];
		for(long jj=0;jj<det.ndet;jj++){

			read_signal_from_fits(samples_struct.fitsvect[0], det.boloname[jj], signal, ns_temp);
			for(long ii = 0; ii< ns_final; ii++)
				signal_bis[ii]=signal[ii];

			long fpixel[2]={1,jj+1};
			//fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status);
			fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, signal_bis, &status);
		}



		fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status);
		fits_copy_header(fptr, outfptr, &status);
		fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);


		int *mask, *mask_bis;
		mask_bis = new int [ns_final];
		for(long jj=0;jj<det.ndet;jj++){

			read_flag_from_fits(samples_struct.fitsvect[0], det.boloname[jj], mask, ns_temp);
			for(long ii = 0; ii< ns_final; ii++)
				mask_bis[ii]=mask[ii];

			long fpixel[2]={1,jj+1};
			//fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status);
			fits_write_pix(outfptr, TINT, fpixel, ns_final, mask_bis, &status);
		}


		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		if (fits_close_file(outfptr, &status))
			fits_report_error(stderr, status);


		if(format_fits==1){ // HIPE format
			cout << "HIPE format found\n";

			// separate HIPE format tables








		}else{ // sanepic format
			cout << "SANEPIC format found\n";
		}

		delete [] RA;
		delete [] DEC;
		delete [] PHI;
		delete [] signal;
		delete [] mask;
		delete [] RA_bis;
		delete [] DEC_bis;
		delete [] PHI_bis;
		delete [] time_bis;
		delete [] signal_bis;
		delete [] mask_bis;

	}

	//delete [] time_min;
	delete [] time;
	delete [] samples_struct.nsamples;

	//	free_dmatrix(image,(long)0, ndet-1, (long)0, ns-1);
	//	delete [] image;

	//	delete [] time_global;
	//delete [] format;
	// separate common tables
	cout << "End of saneSplit\n";

}




