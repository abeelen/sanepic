/*
 * inline_IO2.cpp
 *
 *  Created on: 24 juin 2009
 *      Author: matthieu
 */


#include <cstdlib>
//#include <fcntl.h>
//#include <unistd.h>
//#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "struct_definition.h"
//#include <fftw3.h>
//#include <fstream>
#include <iostream>

#include "inline_IO2.h"

extern "C" {
#include "nrutil.h"
//#include "getdata.h"

}

using namespace std;


bool compute_dirfile_format_file(std::string outdir, struct detectors det, long ntotscan, int rank){

	if(rank==0){
		std::ostringstream oss;
		FILE *fp;
		string temp = outdir + "format";
		if((fp = fopen(temp.c_str(),"w"))!=NULL){
			fprintf(fp,"/ENDIAN little\n/ENCODING none\n/PROTECT none\n/VERSION 7\n");
			fprintf(fp,"\n/INCLUDE Indexes/format\n");
			fprintf(fp,"/INCLUDE Fourier_data/format\n");
			fprintf(fp,"/INCLUDE Noise_data/format\n");
			fclose(fp);
		}

		temp = outdir + "Indexes/format";
		if((fp = fopen(temp.c_str(),"w"))!=NULL){
			fprintf(fp,"/FRAMEOFFSET 1\n/ENDIAN little\n/ENCODING none\n/PROTECT none\n/VERSION 7\n");

			for(long iframe = 0; iframe< ntotscan; iframe++)
				for(long ii =0; ii<det.ndet; ii++){
					oss << "samptopix_" << iframe << "_" << det.boloname[ii] << ".bi";
					std::string temp = oss.str();
					fprintf(fp,"%s RAW INT64 1\n",temp.c_str());
					oss.str("");
				}
		}else{
			cerr << "ERROR : Could not create " << temp << endl;
			return (0);
		}
		fclose(fp);
	}

	return 1;

}

bool compute_dirfile_format_noisePS(std::string outdir, std::vector<string> det, string suffix){


	std::ostringstream oss;
	FILE *fp;
	string temp = outdir + "/Noise_data/format";
	long ndet=(long)det.size();

	if((fp = fopen(temp.c_str(),"w"))!=NULL){
		fprintf(fp,"/FRAMEOFFSET 1\n/ENDIAN little\n/ENCODING none\n/PROTECT none\n/VERSION 7\n");
		for(long ii =0; ii<ndet; ii++){
			// open file
			string filename = det[ii] + "_" + suffix;
			fprintf(fp,"%s RAW INT64 1\n",filename.c_str());
		}
	}else{
		cerr << "ERROR : Could not create " << temp << endl;
		return (0);
	}

	fclose(fp);

	return 1;


}

bool compute_dirfile_format_fdata(std::string outdir, struct detectors det, long ntotscan, int rank){

	if(rank==0){
		std::ostringstream oss;
		FILE *fp;
		string temp = outdir + "/Fourier_data/format";


		if((fp = fopen(temp.c_str(),"w"))!=NULL){
			fprintf(fp,"/FRAMEOFFSET 1\n/ENDIAN little\n/ENCODING none\n/PROTECT none\n/VERSION 7\n");
			for(long iframe = 0; iframe< ntotscan; iframe++)
				for(long ii =0; ii<det.ndet; ii++){
					oss << "fdata_" << iframe << "_" << det.boloname[ii] << ".bi";
					std::string temp = oss.str();
					fprintf(fp,"%s RAW INT64 1\n",temp.c_str());
					oss.str("");
				}
		}else{
			cerr << "ERROR : Could not create " << temp << endl;
			return (0);
		}


		fclose(fp);
	}

	return 1;

}


void write_samptopix(long ns, long long *&samptopix, string outdir, long iframe, std::string boloname) {
	FILE *fp;
	// créer un flux de sortie
	std::ostringstream oss;

	oss << outdir + "/Indexes/samptopix_" << iframe << "_" << boloname << ".bi";

	// récupérer une chaîne de caractères
	std::string temp = oss.str();
	long long ns2 = (long long)ns;

	if((fp = fopen(temp.c_str(),"w"))!=NULL){
		fwrite(&ns2,sizeof(long long),1,fp);
		fwrite(samptopix,sizeof(long long), ns, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << temp << endl;
		exit(0);
	}

	//	size_t nwrite=0;
	oss.str("");
	oss << "samptopix_" << iframe << "_" << boloname << ".bi";
	string temp2=oss.str();
	//	DIRFILE  * dirf;
	//	int test_close=0;

	//	char *buffer;
	//	size_t buflen=100;
	//	buffer = new char [buflen];
	//
	//	char ffname[100];
	//	char *tab;
	//	tab=new char[100];
	//	strcpy(ffname,outdir.c_str());
	//	dirf=dirfile_open(ffname, GD_RDWR | GD_CREAT /*| GD_PRETTY_PRINT | GD_VERBOSE | GD_UNENCODED | GD_LITTLE_ENDIAN*/); // GD_RDWR |GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_PRETTY_PRINT | GD_UNENCODED
	//	get_error_string(dirf, buffer, buflen);
	//	cout << buffer << endl;


	//	dirfile_flush(dirf,  NULL);
	//	get_error_string(dirf, buffer, buflen);
	//	cout << buffer << endl;
	//
	//	getchar();
	//	tab=(char*)dirfilename(dirf);
	//	//	string hh = (string)tab;
	//	cout << tab << endl;
	//	cout << temp2 << endl;
	//	 int g=get_fragment_index(dirf, "samptopix_0_C1.bi");
	//	 cout << "g : " << g << endl;
	//	char **tab2;
	//	tab2=new char*[100];
	//	for(long uu=0;uu<100;uu++)
	//		tab2[uu]=new char[100];
	//	//	dirfile_metaflush(dirf);
	//	tab2=(char**)get_field_list(dirf);
	//	cout << "fields\n";
	//	cout << tab2[0] << endl;
	//	cout << tab2[1] << endl;
	//	cout << tab2[2] << endl;
	//	cout << tab2[3] << endl;
	//	cout << tab2[4] << endl;

	//	if(boloname=="C1"){
	//		int z= dirfile_include(dirf,  temp2.c_str(),  0,  GD_CREAT);
	//		cout << "z " << z << endl;
	//	}
	//	int r= dirfile_validate(dirf, temp2.c_str());
	//	cout << "r " << r << endl;
	//	get_error_string(dirf, buffer, buflen);
	//	cout << buffer << endl;


	//	if(boloname=="C1"){
	//		int i=dirfile_add_raw(dirf, "samptopix_0_C1.bi", GD_INT64,  1,  0);
	//		cout << "i " << i << endl;
	//		get_error_string(dirf, buffer, buflen);
	//		cout << "diraddraw error : " << buffer << endl;
	//		getchar();
	//	}


	//	cout << samptopix[0] << " " << samptopix[1000] << " " << samptopix[29339] << endl;
	//	if((fp = fopen(temp.c_str(),"w"))!=NULL){
	//		int z = dirfile_add_raw(dirf,  temp2.c_str(),  GD_INT64,  1,  0);

	//		dirfile_flush(dirf,  NULL);
	//		get_error_string(dirf, buffer, buflen);
	//		cout << buffer << endl;
	//
	//		if(boloname=="C6"){
	//			cout << "ds C6\n";
	//			gd_entry_t *entry;
	//			entry = new gd_entry_t[1];
	//			entry->field=(char*)temp2.c_str();
	//			entry->field_type=GD_RAW_ENTRY;
	//			entry->data_type=GD_INT64;
	//			char *pt;
	//			pt=new char[1];
	//			*pt='1';
	//			entry->scalar[1]=pt;
	//
	//			cout << "avant add\n";
	//
	//			int z =dirfile_madd(dirf, entry, "samptopix_0_C1.bi");
	//			cout << "z " << z << endl;
	//		}
	//
	//		get_error_string(dirf, buffer, buflen);
	//		cout << buffer << endl;

	//		int t= put_string(dirf,  "", "RAW INT64 1");
	//		get_error_string(dirf, buffer, buflen);
	//		cout << "t " << t << endl;
	//		cout << buffer << endl;

	//		dirfile_flush(dirf,  NULL);
	//		get_error_string(dirf, buffer, buflen);
	//		cout << buffer << endl;
	//		dirfile_close(dirf);
	//		getchar();


	//		nwrite=putdata(dirf, temp2.c_str() , 0, 0, 0, ns,  GD_INT64, samptopix);
	//
	//		get_error_string(dirf, buffer, buflen);
	//		if(strcmp(buffer,"Success")||(nwrite==0)){
	//			cout << "putdata error : " << buffer << endl;
	//			exit(0); // TODO : mettre ca en MPI
	//		}
	//
	//
	//		//	cout << "i : " << i << endl;
	//		//		cout << "written : " << nwrite << endl;
	//
	//		fclose(fp);
	//	}else{
	//		cerr << "ERROR : Could not open " << temp << endl;
	//		exit(0);
	//	}




	//	long long *data_out;
	//	data_out = new long long[ns];
	//	//
	//	//	get_error_string(dirf, buffer, buflen);
	//	//	cout << buffer << endl;
	//	//
	//	size_t read_s = getdata(dirf,  temp2.c_str(),  0,  0,  0,   ns,  GD_INT64,  data_out);
	//		cout << "read : " << read_s << endl;
	//		cout << data_out[0] << " " << data_out[1000] << " " << data_out[29339] << endl;
	//	//
	//	get_error_string(dirf, buffer, buflen);
	//	if(strcmp(buffer,"Success")||(read_s==0)){
	//		cout << "getdata error " << buffer << endl;
	//		exit(0); // TODO : mettre ca en MPI
	//	}
	//	delete [] data_out;

	//	dirfile_close(dirf);
	//	cout << test_close << endl;

#ifdef DEBUG_PRINT
	oss.str("");

	// debug
	oss << outdir + "samptopix_" << iframe << "_" << boloname  << ".txt";
	temp = oss.str();

	if((fp = fopen(temp.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fprintf(fp,"%ld ",ns);
		for(long ii = 0; ii< ns; ii++)
			fprintf(fp,"%lld ",samptopix[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << temp << endl;
		exit(0);
	}
#endif
}


void read_samptopix(long ns, long long *&samptopix, string outdir, long iframe, std::string boloname) {
	FILE *fp;
	size_t result;
	long long ns2;

	//	DIRFILE  * dirf;
	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "/Indexes/samptopix_" << iframe << "_" << boloname  << ".bi";
	std::string testfile = oss.str();
	//	oss.str("");
	//	oss << "samptopix_" << iframe << "_" << boloname << ".bi";
	//	string temp2=oss.str();
	// récupérer une chaîne de caractères


	// Declare buffer in case of errors
	//	char *buffer;
	//	size_t buflen=100;
	//	buffer = new char [buflen];
	//
	//	char ffname[100];
	//	strcpy(ffname,outdir.c_str());
	//	dirf=dirfile_open(ffname, GD_RDWR);
	//
	//	get_error_string(dirf, buffer, buflen);
	//	if(strcmp(buffer,"Success")){
	//		cout << "open dir " << buffer << endl;
	//		exit(0); // TODO : mettre ca en MPI
	//	}

	if((fp = fopen(testfile.c_str(),"r"))){

		//		size_t read_s = getdata(dirf,  temp2.c_str(),  0,  0,  0,   ns,  GD_INT64,  samptopix);
		//
		//		get_error_string(dirf, buffer, buflen);
		//		if(strcmp(buffer,"Success")||(read_s==0)){
		//			cout << "getdata error " << buffer << endl;
		//			exit(0); // TODO : mettre ca en MPI
		//		}
		//		cout << "read : " << read_s << endl;
		//		cout << samptopix[0] << " " << samptopix[1000] << " " << samptopix[29339] << endl;
		result = fread(&ns2,sizeof(long long),1,fp);
		if(ns==ns2)
			result = fread(samptopix,sizeof(long long),ns,fp);
		else{
			fclose(fp);
			cerr << "ERROR : samptopix size is not correct " << ns << " != " << ns2 << endl;
			exit(0);
		}
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile << endl;
		exit(0);
	}

	//	dirfile_close(dirf);
}


void write_indpsrc(long long map_size, long long  npixsrc, long long * indpsrc, std::string outdir){
	FILE *fp;
	string testfile = outdir+"indpsrc.bin";

	if((fp = fopen(testfile.c_str(),"w"))!=NULL){
		fwrite(&map_size,sizeof(long long),1,fp);
		fwrite(&npixsrc, sizeof(long long),1,fp);
		fwrite(indpsrc,sizeof(long long), map_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile << endl;
		exit(0);
	}
}


void  read_indpsrc(long long &map_size, long long &npixsrc, long long *&indpsrc, std::string outdir){
	FILE *fp;
	string testfile;
	size_t result;

	testfile = outdir + "indpsrc.bin";
	if ((fp = fopen(testfile.c_str(),"r"))!=NULL){
		result = fread(&map_size,sizeof(long long),1,fp);
		result = fread(&npixsrc,sizeof(long long),1,fp);
		indpsrc = new long long[map_size];
		result = fread(indpsrc,sizeof(long long), map_size, fp);
		fclose(fp);
	}else{
		cerr << "Error : cannot find indpsrc.bin file at " << testfile << endl;
		exit(0);
	}
}


void write_indpix(long long ind_size, long long npix, long long *indpix, string outdir, int flagon) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad.bi";

	if((fp = fopen(testfile2.c_str(),"w"))!=NULL){
		fwrite(&flagon,sizeof(int),1,fp); // mat 04/06
		fwrite(&npix,sizeof(long long),1,fp);
		fwrite(&ind_size,sizeof(long long),1,fp);
		fwrite(indpix,sizeof(long long), ind_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile2 << endl;
		exit(0);
	}

#ifdef DEBUG_PRINT
	//Debug
	testfile2 = outdir + "Indpix_for_conj_grad.txt";
	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fprintf(fp,"%d\n",flagon); // mat 04/06
		fprintf(fp,"%lld\n",npix);
		fprintf(fp,"%lld\n",ind_size);
		for(long ii =0;ii<ind_size;ii++)
			fprintf(fp,"%lld ",indpix[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}

#endif
}

void read_indpix(long long &ind_size, long long &npix, long long *&indpix, string outdir, int &flagon) {
	FILE *fp;
	string testfile2;
	size_t result;

	testfile2 = outdir + "Indpix_for_conj_grad.bi";
	if ((fp = fopen(testfile2.c_str(),"r"))!=NULL){
		result = fread(&flagon,sizeof(int),1,fp); // mat 04/06
		result = fread(&npix,sizeof(long long),1,fp);
		result = fread(&ind_size,sizeof(long long),1,fp);
		indpix=new long long[ind_size];
		result = fread(indpix,sizeof(long long), ind_size, fp);
		fclose(fp);
	}else{
		cerr << "Error : cannot find Indpix file " << testfile2 << endl;
		exit(0);
	}
}

void write_PNd(double *PNd, long long npix,  string outdir) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "PNdCorr.bi";

	if((fp = fopen(testfile2.c_str(),"w"))!=NULL){
		//fprintf(fp,"%d\n",npix);
		//for(long ii=0;ii<npix;ii++)
		//fprintf(fp,"%lf ",PNd[ii]);
		fwrite(&npix,sizeof(long long),1,fp);
		fwrite(PNd,sizeof(double), npix, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}

#ifdef DEBUG_PRINT
	//Debug
	testfile2 = outdir + "PNdCorr.txt";

	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		//fprintf(fp,"%d\n",npix);
		//for(long ii=0;ii<npix;ii++)
		//fprintf(fp,"%lf ",PNd[ii]);
		fprintf(fp,"%lld\n",npix);
		for(int ii= 0; ii<npix;ii++)
			fprintf(fp,"%lf ",PNd[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}

#endif
}

// correlation between npix here and npix in Indpix file is done in sanePic (main.cpp)
void read_PNd(double *&PNdtot, long long &npix,  string outdir) {
	FILE *fp;
	string testfile2;
	size_t result;

	testfile2 = outdir + "PNdCorr.bi";

	if ((fp = fopen(testfile2.c_str(),"r"))!=NULL){
		result = fread(&npix,sizeof(long long),1,fp);
		PNdtot= new double[npix]; // mat 04/06
		result = fread(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
	}else{
		cerr << "Error. Unable to read file : " << testfile2 << endl;
		exit(0);
	}

}

void write_fdata(long ns, fftw_complex *fdata, string prefixe, string outdir, long idet, long iframe, std::vector<std::string> bolonames) {
	FILE *fp;
	double data_size;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "/Fourier_data/" + prefixe << iframe << "_" << bolonames[idet] << ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		data_size = (double)((ns/2+1)*2);
		fwrite(&data_size,sizeof(double),1,fp);
		fwrite(fdata,sizeof(double), (long) data_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not write " << testfile << endl;
		exit(0);
	}

#ifdef DEBUG_PRINT
	// Debug
	oss.str("");
	oss << outdir + "fdata_" << iframe << "_" << bolonames[idet] << ".txt";

	// récupérer une chaîne de caractères
	testfile = oss.str();
	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		data_size = (ns/2+1)*2;

		fprintf(fp,"%ld ",data_size);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		for(int ii=0;ii<(ns/2+1);ii++){
			fprintf(fp,"%lf ",fdata[ii][0]);
			fprintf(fp,"%lf \n",fdata[ii][1]);}

		//cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile << endl;
		exit(0);
	}
#endif

}


void read_fdata(long ns, fftw_complex *&fdata, string prefixe,  string outdir, long idet, long iframe, std::vector<std::string> bolonames) {
	FILE *fp;
	size_t result;
	double data_size;
	long data_read;

	// créer un flux de sortie
	std::ostringstream oss;

	oss << outdir + "/Fourier_data/" + prefixe << iframe << "_" << bolonames[idet] << ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"r"))!=NULL){
		result = fread(&data_size,sizeof(double),1,fp);
		data_read=(long) data_size;
		if (data_read!=(ns/2+1)*2){
			cerr << "Error. fdata size does not correspond to expected size : " << data_read << " != " << (ns/2+1)*2 << endl;
			fclose(fp);
			exit(1);
		}
		result = fread(fdata,sizeof(double), data_read, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}
}



//void write_fPs(long ns, fftw_complex *fdata, string outdir, long idet, long iframe, std::vector<std::string> bolonames) {
//	FILE *fp;
//	long data_size;
//
//	// créer un flux de sortie
//	std::ostringstream oss;
//	oss << outdir + "fPs_" << iframe << "_" << bolonames[idet] << ".bi";
//
//	// récupérer une chaîne de caractères
//	std::string testfile = oss.str();
//
//	if((fp = fopen(testfile.c_str(),"w"))!=NULL){
//		data_size = (ns/2+1)*2;
//		fwrite(&data_size,sizeof(long),1,fp);
//		fwrite(fdata,sizeof(double), data_size, fp);
//		fclose(fp);
//	}else{
//		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
//		exit(1);
//	}
//
//#ifdef DEBUG_PRINT
//	oss.str("");
//	oss << outdir + "fPs_" << iframe << "_" << bolonames[idet] << ".txt";
//
//	// récupérer une chaîne de caractères
//	testfile = oss.str();
//	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
//		data_size = (ns/2+1)*2;
//
//		fprintf(fp,"%ld ",data_size);
//		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
//		for(int ii=0;ii<(ns/2+1);ii++){
//			fprintf(fp,"%lf ",fdata[ii][0]);
//			fprintf(fp,"%lf ",fdata[ii][1]);}
//
//		//cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
//		fclose(fp);
//	}else{
//		cerr << "ERROR : Could not open " << testfile << endl;
//		exit(0);
//	}
//#endif
//
//}

void read_mixmat_txt(string MixMatfile, long ndet, long ncomp, double **&mixmat)
{
	//TODO : CAN NOT WORK !
	FILE *fp;
	int result;
	long ncomp2;
	double dummy1; // used to read mixing matrix

	if ((fp = fopen(MixMatfile.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
		exit(1);
	}
	result = fscanf(fp,"%ld",&ncomp2); // modified d => ld to avoid warning, mat-27/05

	mixmat = dmatrix(0,ndet-1,0,ncomp-1); // ajout mat 24/07 pour eviter de declarer 20 comp useless

	for (long ii=0;ii<ndet;ii++){
		for (long jj=0;jj<ncomp2;jj++){
			result = fscanf(fp,"%lf",&dummy1);
			if (jj< ncomp)
				mixmat[ii][jj] = dummy1;
		}
	}
	fclose(fp);

}

