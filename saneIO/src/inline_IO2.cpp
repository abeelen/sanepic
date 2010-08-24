

#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "struct_definition.h"
#include "inputFileIO.h"
#include <iostream>

#include "inline_IO2.h"


extern "C" {
#include "nrutil.h"
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
// TODO : change this for multiple suffix (means multiple covariance matrix)
bool compute_dirfile_format_noisePS(std::string outdir, std::vector<string> det, string suffix)
/*!  Create the format file for kst and write the noise file list inside  */
{

	std::ostringstream oss;
	FILE *fp;
	string temp = outdir + "Noise_data/format";
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

bool compute_dirfile_format_fdata(std::string outdir, struct detectors det, long ntotscan, int rank)
/*!  Create the format file for kst and write the fourier data list inside  */
{

	if(rank==0){
		std::ostringstream oss;
		FILE *fp;
		string temp = outdir + "Fourier_data/format";


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


int write_samptopix(long ns, long long *&samptopix, string outdir, string filename, std::string boloname)
/*!  write a sample to pixel vector to disk  */
{
	FILE *fp;
	// créer un flux de sortie
	//	std::ostringstream oss;
	string outfile;



	outfile=outdir + "/Indexes/samptopix_" + FitsBasename(filename) + "_" + boloname + ".bi";
	// Transform into string

	long long ns2 = (long long)ns;

	if((fp = fopen(outfile.c_str(),"w"))!=NULL){
		fwrite(&ns2,sizeof(long long),1,fp);
		fwrite(samptopix,sizeof(long long), ns, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << outfile << endl;
		return 1;
	}


#ifdef DEBUG_PRINT
	//	oss.str("");
	//
	//	// debug
	//	oss << outdir + "samptopix_" << FitsBasename(filename) << "_" << boloname  << ".txt";
	//	temp = oss.str();

	outfile=outdir + "/Indexes/samptopix_" + FitsBasename(filename) + "_" + boloname + ".txt";

	if((fp = fopen(outfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		fprintf(fp,"%ld ",ns);
		for(long ii = 0; ii< ns; ii++)
			fprintf(fp,"%lld ",samptopix[ii]);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << outfile << endl;
		exit(0);
	}
#endif

	return 0;
}


int read_samptopix(long ns, long long *&samptopix, string outdir, string filename, std::string boloname)
/*!  read a sample to pixel vector from disk  */
{
	FILE *fp;
	size_t result;
	long long ns2;

	//	std::ostringstream oss;
	//	oss << outdir + "/Indexes/samptopix_" << iframe << "_" << boloname  << ".bi";
	//	std::string testfile = oss.str();

	string outfile=outdir + "/Indexes/samptopix_" + FitsBasename(filename) + "_" + boloname + ".bi";

	if((fp = fopen(outfile.c_str(),"r"))){

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
		cerr << "ERROR : Could not find " << outfile << endl;
		return 1;
	}

	return 0;
}

int write_indpsrc(long long map_size, long long  npixsrc, long long * indpsrc, std::string outdir)
/*!  write the bright sources index to disk  */
{
	FILE *fp;
	string testfile = outdir+"indpsrc.bin";

	if((fp = fopen(testfile.c_str(),"w"))!=NULL){
		fwrite(&map_size,sizeof(long long),1,fp);
		fwrite(&npixsrc, sizeof(long long),1,fp);
		fwrite(indpsrc,sizeof(long long), map_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << testfile << endl;
		return 1;
	}

	return 0;
}


int  read_indpsrc(long long &map_size, long long &npixsrc, long long *&indpsrc, std::string outdir)
/*!  read the bright sources index from disk  */
{
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
		return 1;
	}

	return 0;
}


int write_indpix(long long ind_size, long long npix, long long *indpix, string outdir, int flagon)
/*!  write the map index to disk  */
{
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
		return 1;
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

	return 0;
}

int read_indpix(long long &ind_size, long long &npix, long long *&indpix, string outdir, int &flagon)
/*!  read the map index from disk  */
{
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
		return 1;
	}

	return 0;
}

int write_PNd(double *PNd, long long npix,  string outdir)
/*!  write the map preconditioner to disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "PNdCorr.bi";

	if((fp = fopen(testfile2.c_str(),"w"))!=NULL){
		fwrite(&npix,sizeof(long long),1,fp);
		fwrite(PNd,sizeof(double), npix, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
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
	return 0;
}

// correlation between npix here and npix in Indpix file is done in sanePic (main.cpp)
int read_PNd(double *&PNdtot, long long &npix,  string outdir)
/*!  read the map preconditioner from disk  */
{
	FILE *fp;
	string testfile2;
	size_t result;

	testfile2 = outdir + "PNdCorr.bi";

	if ((fp = fopen(testfile2.c_str(),"r"))!=NULL){
		result = fread(&npix,sizeof(long long),1,fp);
		PNdtot= new double[npix];
		result = fread(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
	}else{
		cerr << "Error. Unable to read file : " << testfile2 << endl;
		return 1;
	}

	return 0;
}

int write_fdata(long ns, fftw_complex *fdata, string prefixe, string outdir, long idet, string filename, std::vector<std::string> bolonames)
/*! write Fourier data file to disk */
{
	FILE *fp;
	double data_size;

	//	std::ostringstream oss;
	//	oss << outdir + "/Fourier_data/" + prefixe << iframe << "_" << bolonames[idet] << ".bi";

	// Transform into string
	std::string outfile;

	outfile=outdir + "/Fourier_data/" + prefixe + FitsBasename(filename) + "_" + bolonames[idet] + ".bi";

	if((fp = fopen(outfile.c_str(),"w"))){
		data_size = (double)((ns/2+1)*2);
		fwrite(&data_size,sizeof(double),1,fp);
		fwrite(fdata,sizeof(double), (long) data_size, fp);
		fclose(fp);
	}else{
		cerr << "ERROR : Could not write " << outfile << endl;
		return 1;
	}

#ifdef DEBUG_PRINT
	// Debug
	//	oss.str("");
	//	oss << outdir + "fdata_" << iframe << "_" << bolonames[idet] << ".txt";

	// Transform into string
	//	testfile = oss.str();
	outfile=outdir + "/Fourier_data/" + prefixe + FitsBasename(filename) + "_" + bolonames[idet] + ".txt";
	if((fp = fopen(outfile.c_str(),"w"))){ // doubles parenthèses sinon warning ...
		data_size = (ns/2+1)*2;

		fprintf(fp,"%ld ",data_size);
		//if (data_size!=(ns/2+1)*2) cerr << "Error. fdata size does not correspond to expected size\n";
		for(int ii=0;ii<(ns/2+1);ii++){
			fprintf(fp,"%lf ",fdata[ii][0]);
			fprintf(fp,"%lf \n",fdata[ii][1]);}

		//cout << "writing fdata  : "  << (ns/2+1)*2 << " " << sizeof(double) << endl;
		fclose(fp);
	}else{
		cerr << "ERROR : Could not open " << outfile << endl;
		exit(0);
	}
#endif

	return 0;
}


int read_fdata(long ns, fftw_complex *&fdata, string prefixe,  string outdir, long idet, string filename, std::vector<std::string> bolonames)
/*!  read the map preconditioner from disk  */
{
	FILE *fp;
	size_t result;
	double data_size;
	long data_read;

	string outfile=outdir + "/Fourier_data/" + prefixe + FitsBasename(filename) + "_" + bolonames[idet] + ".bi";

	if((fp = fopen(outfile.c_str(),"r"))!=NULL){
		result = fread(&data_size,sizeof(double),1,fp);
		data_read=(long) data_size;
		if (data_read!=(ns/2+1)*2){
			cerr << "Error. fdata size does not correspond to expected size : " << data_read << " != " << (ns/2+1)*2 << endl;
			fclose(fp);
			return 1;
		}
		result = fread(fdata,sizeof(double), data_read, fp);
		fclose(fp);
	}else{
		cerr << "ERROR: Can't find Fourier transform data file" << outfile << ". Exiting. \n";
		return 1;
	}

	return 0;
}


int read_mixmat_txt(string MixMatfile, long ndet, long ncomp, double **&mixmat)
{
	//TODO : CAN NOT WORK !
	FILE *fp;
	int result;
	long ncomp2;
	double dummy1; // used to read mixing matrix

	if ((fp = fopen(MixMatfile.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
		cout << "Advice : verify the file is in your noise directory and that his name is : " << MixMatfile << endl;
		return 1;
	}
	result = fscanf(fp,"%ld",&ncomp2); // modified d => ld to avoid warning, mat-27/05

	mixmat = dmatrix(0,ndet-1,0,ncomp-1);

	for (long ii=0;ii<ndet;ii++){
		for (long jj=0;jj<ncomp2;jj++){
			result = fscanf(fp,"%lf",&dummy1);
			if (jj< ncomp)
				mixmat[ii][jj] = dummy1;
		}
	}
	fclose(fp);
	return 0;

}

