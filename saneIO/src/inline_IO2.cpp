/*
 * inline_IO2.cpp
 *
 *  Created on: 24 juin 2009
 *      Author: matthieu
 */

#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include "inline_IO2.h"

extern "C" {
#include "nrutil.h"
}

//using namespace std;

// sanePos functions
void write_info_pointing(int nn, string outdir, string termin, int coordsyst, double *tanpix, double *tancoord) {
	FILE *fp;
	string testfile2;


	testfile2 = outdir + "InfoPointing_for_Sanepic_" + termin  + ".txt";
	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fprintf(fp,"%d\n",nn);
		fprintf(fp,"%d\n",coordsyst);
		fprintf(fp,"%lf\n",tanpix[0]);
		fprintf(fp,"%lf\n",tanpix[1]);
		fprintf(fp,"%lf\n",tancoord[0]);
		fprintf(fp,"%lf\n",tancoord[1]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile2 << endl;
		exit(0);
	}
}

void write_samptopix(long ns, long *samptopix, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	// créer un flux de sortie
	std::ostringstream oss;


	oss << outdir + "samptopix_" << iframe << "_" << idet << "_" + termin + ".bi";

	// récupérer une chaîne de caractères
	std::string temp = oss.str();

	if((fp = fopen(temp.c_str(),"w"))!=NULL){
		//fwrite(&ns,sizeof(long),fp);
		fwrite(samptopix,sizeof(long), ns, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << temp << endl;
		exit(0);
	}

	oss.str("");

	// debug
	oss << outdir + "samptopix_" << iframe << "_" << idet << "_" + termin + ".txt";
	temp = oss.str();

	if((fp = fopen(temp.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		for(int ii = 0; ii< ns; ii++)
			fprintf(fp,"%ld ",samptopix[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << temp << endl;
		exit(0);
	}
}


void write_indpix(long ind_size, int npix, long *indpix, string termin, string outdir, int flagon) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad_" + termin + ".bi";

	if((fp = fopen(testfile2.c_str(),"w"))!=NULL){
		fwrite(&flagon,sizeof(int),1,fp); // mat 04/06
		fwrite(&npix,sizeof(int),1,fp);
		fwrite(&ind_size,sizeof(long),1,fp);
		fwrite(indpix,sizeof(long), ind_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile2 << endl;
		exit(0);
	}

	//Debug
	testfile2 = outdir + "Indpix_for_conj_grad_" + termin + ".txt";
	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fprintf(fp,"%d\n",flagon); // mat 04/06
		fprintf(fp,"%d\n",npix);
		fprintf(fp,"%ld\n",ind_size);
		for(int ii =0;ii<ind_size;ii++)
			fprintf(fp,"%ld ",indpix[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}
}


//sanePre functions

void read_info_pointing(int &nn, string outdir, string termin, int &coordsyst2, double *&tanpix, double *&tancoord) {
	FILE *fp;
	string testfile2;


	testfile2 = outdir + "InfoPointing_for_Sanepic_" + termin + ".txt";
	if((fp = fopen(testfile2.c_str(),"r"))){ // doubles parenthèses sinon warning ...
		fscanf(fp,"%d\n",&nn);
		fscanf(fp,"%d\n",&coordsyst2);
		fscanf(fp,"%lf\n",tanpix);
		fscanf(fp,"%lf\n",tanpix+1);
		fscanf(fp,"%lf\n",tancoord);
		fscanf(fp,"%lf\n",tancoord+1);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}
}

void read_indpix(long &ind_size, int &npix, long *&indpix, string termin, string outdir, int &flagon) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad_" + termin + ".bi";
	if ((fp = fopen(testfile2.c_str(),"r"))!=NULL){
		fread(&flagon,sizeof(int),1,fp); // mat 04/06
		fread(&npix,sizeof(int),1,fp);
		fread(&ind_size,sizeof(long),1,fp);
		indpix=new long[ind_size];
		fread(indpix,sizeof(long), ind_size, fp);
		fclose(fp);
	}else{
		cerr << "Error : cannot find Indpix file " << testfile2 << endl;
		exit(0);
	}
}

void write_PNd(double *PNd, int npix, string termin, string outdir) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "PNdCorr_" + termin + ".bi";

	if((fp = fopen(testfile2.c_str(),"w"))!=NULL){
		//fprintf(fp,"%d\n",npix);
		//for(long ii=0;ii<npix;ii++)
		//fprintf(fp,"%lf ",PNd[ii]);
		fwrite(&npix,sizeof(int),1,fp);
		fwrite(PNd,sizeof(double), npix, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}


	//Debug
	testfile2 = outdir + "PNdCorr_" + termin + ".txt";

	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		//fprintf(fp,"%d\n",npix);
		//for(long ii=0;ii<npix;ii++)
		//fprintf(fp,"%lf ",PNd[ii]);
		fprintf(fp,"%d\n",npix);
		for(int ii= 0; ii<npix;ii++)
			fprintf(fp,"%lf ",PNd[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}
}


void read_PNd(double *&PNdtot, int &npix, string termin, string outdir) {
	FILE *fp;
	string testfile2;
	testfile2 = outdir + "PNdCorr_" + termin + ".bi";

	if ((fp = fopen(testfile2.c_str(),"r"))!=NULL){
		fread(&npix,sizeof(int),1,fp);
		PNdtot= new double[npix]; // mat 04/06
		fread(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
	}else{
		cerr << "Error. Unable to find file : " << testfile2 << endl;
		exit(0);
	}

}
void read_samptopix(long ns, long *&samptopix, string termin, string outdir, int idet, long iframe) {
	FILE *fp;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "samptopix_" << iframe << "_" << idet << "_" + termin + ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"r"))){
		//read(&ns,sizeof(long),fp);
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile << endl;
		exit(0);
	}
}


void write_fdata(long ns, fftw_complex *fdata, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	//long data_size;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "fdata_" << iframe << "_" << idet << "_" + termin + ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		//data_size = (ns/2+1)*2;
		//fwrite(&data_size,sizeof(long),fp);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		//cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile << endl;
		exit(0);
	}

	// Debug
	oss.str("");
	oss << outdir + "fdata_" << iframe << "_" << idet << "_" + termin + ".txt";

	// récupérer une chaîne de caractères
	testfile = oss.str();
	if((fp = fopen(testfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		//data_size = (ns/2+1)*2;
		//fwrite(&data_size,sizeof(long),fp);
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

}


void read_noise_file(long &nbins, double *&ell, double **&SpN_all, string nameSpfile, long ndet) {

	FILE *fp;
	double dnbins;


	if ((fp = fopen(nameSpfile.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find noise power spectra file" << nameSpfile << " , check -k or -K in command line. Exiting. \n";
		exit(1);
	}
	fread(&dnbins,sizeof(double), 1, fp);
	nbins = (long)dnbins;
	SpN_all = dmatrix((long)0,ndet-1,(long)0,nbins-1);
	ell = new double[nbins+1];
	fread(ell,sizeof(double), nbins+1, fp);
	fread(*SpN_all,sizeof(double), nbins*ndet, fp);
	fclose(fp);


}


void read_fdata(long ns, fftw_complex *&fdata, string prefixe, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	//long data_size;

	// créer un flux de sortie
	std::ostringstream oss;

	oss << outdir + prefixe << iframe << "_" << idet << "_" + termin + ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"r"))!=NULL){
		//fread(&data_size,sizeof(long),fp);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		fread(fdata,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}
}

void write_fPs(long ns, fftw_complex *fdata, string termin, string outdir, long idet, long iframe) {
	FILE *fp;
	//long data_size;

	// créer un flux de sortie
	std::ostringstream oss;
	oss << outdir + "fPs" << iframe << "_" << idet << "_" + termin + ".bi";

	// récupérer une chaîne de caractères
	std::string testfile = oss.str();

	if((fp = fopen(testfile.c_str(),"w"))!=NULL){
		//data_size = (ns/2+1)*2;
		//fwrite(&data_size,sizeof(long),fp);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}
}

void write_info_for_second_part(string outdir, string termin, int nn, int npix,
		double pixdeg, double *tancoord, double* tanpix, int coordsyst, bool flagon, long* indpix){

	//char testfile[100];
	FILE *fp;
	string testfile;

	testfile = outdir + "InfoFor2ndStep_" + termin + ".txt";
	//sprintf(testfile,"%s%s%s%s",outdir.c_str(),"InfoFor2ndStep_",termin.c_str(),".txt");
	if((fp = fopen(testfile.c_str(),"w"))==NULL){
		cerr << "Cannot open file :" << testfile << "\tExiting." << endl;
		exit(1);
	}
	fprintf(fp,"%d\n",nn);
	fprintf(fp,"%d\n",npix);
	fprintf(fp,"%lf\n",pixdeg);
	fprintf(fp,"%lf\n",tancoord[0]);
	fprintf(fp,"%lf\n",tancoord[1]);
	fprintf(fp,"%lf\n",tanpix[0]);
	fprintf(fp,"%lf\n",tanpix[1]);
	fprintf(fp,"%d\n",coordsyst);
	fprintf(fp,"%i\n",flagon);
	fprintf(fp,"\n");
	for (long ii=0;ii<nn*nn;ii++)
		fprintf(fp,"%ld\n",indpix[ii]);
	fclose(fp);

}

void read_mixmat_txt(string MixMatfile, long ndet, long ncomp2, double **&mixmat)
{
	FILE *fp;
	double dummy1; // used to read mixing matrix

	if ((fp = fopen(MixMatfile.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
		exit(1);
	}
	fscanf(fp,"%ld",&ncomp2); // modified d => ld to avoid warning, mat-27/05

	mixmat = dmatrix(0,ndet-1,0,ncomp2-1); // ajout mat 24/07 pour eviter de declarer 20 comp useless

	for (long ii=0;ii<ndet;ii++){
		for (long jj=0;jj<ncomp2;jj++){
			fscanf(fp,"%lf",&dummy1);
			mixmat[ii][jj] = dummy1;
		}
	}
	fclose(fp);

}
