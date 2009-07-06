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
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}
}

void write_samptopix(long ns, long *samptopix, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	string testfile2;
	char testfile[100];

	sprintf(testfile,"%ld",iframe);
	testfile2 = outdir + "samptopix_";// + char(iframe) + "_" + char(idet) + "_" + termin + ".bi";
	testfile2 +=testfile;
	testfile2 +="_";
	sprintf(testfile,"%d",idet);
	testfile2 +=testfile;
	testfile2 += "_" + termin + ".bi";

	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fwrite(samptopix,sizeof(long), ns, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		exit(0);
	}
}


void write_indpix(long ind_size, int npix, long *indpix, string termin, string outdir, int flagon) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad_" + termin + ".bi";

	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fwrite(&flagon,sizeof(int),1,fp); // mat 04/06
		fwrite(&npix,sizeof(int),1,fp);
		fwrite(indpix,sizeof(long), ind_size, fp);
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

void read_indpix(long ind_size, int &npix, long *&indpix, string termin, string outdir, int &flagon) {
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad_" + termin + ".bi";
	if ((fp = fopen(testfile2.c_str(),"r"))!=NULL){
		fread(&flagon,sizeof(int),1,fp); // mat 04/06
		fread(&npix,sizeof(int),1,fp);
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

	if((fp = fopen(testfile2.c_str(),"w"))){ // doubles parenthèses sinon warning ...
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
}


void read_PNd(double *&PNdtot, int &npix, string termin, string outdir) {
	FILE *fp;
	char testfile[100];

	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"PNdCorr_",termin.c_str(),".bi");
	if ((fp = fopen(testfile,"r"))!=NULL){
		fread(&npix,sizeof(int),1,fp);
		PNdtot= new double[npix]; // mat 04/06
		//PNd= new double[npix]; // mat 04/06
		fread(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
	}else{
		cerr << "erreur : impossible de trouver le fichier " << testfile << endl;
		exit(0);
	}

}
void read_samptopix(long ns, long *&samptopix, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	char testfile[100];

	sprintf(testfile,"%s%s%ld%s%d%s%s%s",outdir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
	if((fp = fopen(testfile,"r"))){
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile << endl;
		exit(0);
	}
}


void write_fdata(long ns, fftw_complex *fdata, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	char testfile[100];

	sprintf(testfile,"%s%s%ld%s%d%s%s%s",outdir.c_str(),"fdata_",iframe,"_",idet,"_",termin.c_str(),".bi");

	if((fp = fopen(testfile,"w"))){ // doubles parenthèses sinon warning ...
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile << endl;
		exit(0);
	}
}


void read_noise_file(long &nbins, double *&ell, double *&SpN, double **&SpN_all, string nameSpfile, long ndet) {

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
	SpN = new double[nbins];
	fread(ell,sizeof(double), nbins+1, fp);
	fread(*SpN_all,sizeof(double), nbins*ndet, fp);
	fclose(fp);

}


void read_fdata(long ns, fftw_complex *&fdata, string prefixe, string termin, string outdir, int idet, long iframe) {
	FILE *fp;
	char testfile[100];

	sprintf(testfile,"%s%s%s%ld%s%d%s%s%s",outdir.c_str(),prefixe.c_str(),"_",iframe,"_",idet,"_",termin.c_str(),".bi");
	if((fp = fopen(testfile,"r"))!=NULL){
		fread(fdata,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}
}

void write_fPs(long ns, fftw_complex *fdata, string termin, string outdir, long idet, long iframe) {
	FILE *fp;
	char testfile[100];

	sprintf(testfile,"%s%s%ld%s%ld%s%s%s",outdir.c_str(),"fPs_",iframe,"_",idet,"_",termin.c_str(),".bi");
	if((fp = fopen(testfile,"w"))!=NULL){
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << testfile << ". Exiting. \n";
		exit(1);
	}
}

void write_info_for_second_part(string outdir, string termin, int nn, int npix,
		double pixdeg, double *tancoord, double* tanpix, int coordsyst, bool flagon, long* indpix){

	char testfile[100];
	FILE *fp;

	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"InfoFor2ndStep_",termin.c_str(),".txt");
	if((fp = fopen(testfile,"w"))==NULL){
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


