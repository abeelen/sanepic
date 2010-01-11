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

//#include <fftw3.h>
//#include <fstream>
#include <iostream>

#include "inline_IO2.h"

extern "C" {
#include "nrutil.h"

}

using namespace std;



void write_samptopix(long ns, long long *&samptopix, string outdir, long idet, long iframe, std::vector<string> bolonames) {
	FILE *fp;
	// créer un flux de sortie
	std::ostringstream oss;

	oss << outdir + "samptopix_" << iframe << "_" << bolonames[idet] << ".bi";

	// récupérer une chaîne de caractères
	std::string temp = oss.str();

	if((fp = fopen(temp.c_str(),"w"))!=NULL){
		fwrite(&ns,sizeof(long),1,fp);
		fwrite(samptopix,sizeof(long long), ns, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << temp << endl;
		exit(0);
	}


#ifdef DEBUG_PRINT
	oss.str("");

	// debug
	oss << outdir + "samptopix_" << iframe << "_" << bolonames[idet]  << ".txt";
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


void read_samptopix(long ns, long long *&samptopix, string outdir, long idet, long iframe, std::vector<std::string> bolonames) {
	FILE *fp;
	size_t result;
	long ns2;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "samptopix_" << iframe << "_" << bolonames[idet]  << ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"r"))){
		result = fread(&ns2,sizeof(long),1,fp);
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

void write_fdata(long ns, fftw_complex *fdata, string outdir, long idet, long iframe, std::vector<std::string> bolonames) {
	FILE *fp;
	long data_size;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "fdata_" << iframe << "_" << bolonames[idet] << ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		data_size = (ns/2+1)*2;
		fwrite(&data_size,sizeof(long),1,fp);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		//cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
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
	long data_size;

	// créer un flux de sortie
	std::ostringstream oss;

	oss << outdir + prefixe << iframe << "_" << bolonames[idet] << ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"r"))!=NULL){
		result = fread(&data_size,sizeof(long),1,fp);
		if (data_size!=(ns/2+1)*2){
			cerr << "Error. fdata size does not correspond to expected size : " << data_size << " != " << (ns/2+1)*2 << endl;
			fclose(fp);
			exit(1);
		}
		result = fread(fdata,sizeof(double), data_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}
}


void write_fPs(long ns, fftw_complex *fdata, string outdir, long idet, long iframe, std::vector<std::string> bolonames) {
	FILE *fp;
	long data_size;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "fPs_" << iframe << "_" << bolonames[idet] << ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"w"))!=NULL){
		data_size = (ns/2+1)*2;
		fwrite(&data_size,sizeof(long),1,fp);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		fwrite(fdata,sizeof(double), data_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}

#ifdef DEBUG_PRINT
	oss.str("");
	oss << outdir + "fPs_" << iframe << "_" << bolonames[idet] << ".txt";

	// récupérer une chaîne de caractères
	testfile = oss.str();
	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		data_size = (ns/2+1)*2;

		fprintf(fp,"%ld ",data_size);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		for(int ii=0;ii<(ns/2+1);ii++){
			fprintf(fp,"%lf ",fdata[ii][0]);
			fprintf(fp,"%lf ",fdata[ii][1]);}

		//cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile << endl;
		exit(0);
	}
#endif

}

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

