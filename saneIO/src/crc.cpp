#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iostream>

#include "struct_definition.h"
#include "imageIO.h"

#include "crc.h"

extern "C" {
#include "wcslib/wcshdr.h"
}

using namespace std;

unsigned checksum(void *buffer, size_t len, unsigned int seed)
{
	unsigned char *buf = (unsigned char *)buffer;
	size_t i;

	for (i = 0; i < len; ++i)
		seed += (unsigned int)(*buf++);
	return seed;
}


void compute_checksum(std::string ini_file, std::string tmp_dir, long long npix, long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk)
{

	FILE *fp;
	size_t len;
	string file;
	char buf[6144];


	if (NULL == (fp = fopen(ini_file.c_str(), "r")))
	{
		cout << "Unable to open " << ini_file << " for reading\n";
		return;
	}
	len = fread(buf, sizeof(char), sizeof(buf), fp); // TODO ameliorer ca !
	char buf2[len-116];
	for(long hh=0;hh<(long)(len-116);hh++)
		buf2[hh]=buf[hh];
	fclose(fp);

	chk.chk_ini_file=checksum(buf2, len-116, 0);
	//	printf("The checksum of %s is %u\n", ini_file.c_str(), chk.chk_ini_file);

	file= tmp_dir + "mapHeader.keyrec";
	if (NULL == (fp = fopen(file.c_str(), "r")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return;
	}
	len = fread(buf, sizeof(char), sizeof(buf), fp);
	fclose(fp);

	chk.chk_wcs_file=checksum(buf, len, 0);
	//	printf("The checksum of %s is %u\n", file.c_str(), chk.chk_wcs_file);

	chk.chk_indpix=checksum(indpix, (size_t) npix, 0);
	//	printf("The checksum of Indpix is %u\n", chk.chk_indpix);

	chk.chk_indpsrc=checksum(indpsrc, (size_t) indpsrc_size, 0);
	//	printf("The checksum of Indpsrc is %u\n", chk.chk_indpsrc);

}


int write_checksum(std::string tmp_dir, struct checksum chk, std::string projectname)
{

	FILE* fp;
	string file = tmp_dir + projectname + "_checksum.bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "wb")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return 1;
	}
	len+=fwrite(&chk.chk_ini_file,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_wcs_file,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_indpix,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_indpsrc,sizeof(unsigned int),1,fp);

	fclose(fp);

	return 0;
}

void read_checksum(std::string tmp_dir, struct checksum &chk, std::string projectname)
{

	FILE* fp;
	string file = tmp_dir + projectname + "_checksum.bin";
	size_t len=0;


	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return;
	}
	len+=fread(&chk.chk_ini_file,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_wcs_file,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_indpix,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_indpsrc,sizeof(unsigned int),1,fp);

	fclose(fp);

}

bool compare_checksum(struct checksum chk_t, struct checksum chk_t2){

	if((chk_t.chk_ini_file != chk_t2.chk_ini_file) ||
			(chk_t.chk_wcs_file!=chk_t2.chk_wcs_file) ||
			(chk_t.chk_indpix!=chk_t2.chk_indpix) ||
			(chk_t.chk_indpsrc!=chk_t2.chk_indpsrc))
		return 1;

	return 0;
}

int load_idupl(string tmp_dir, string out_dir, int &idupl){
	FILE* fp;
	string file = tmp_dir + "data_sanePic.bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return 1;
	}
	len+=fread(&idupl,sizeof(int),1,fp);

	fclose(fp);

	return 0;
}

int load_from_disk(string tmp_dir, string out_dir, double *S, double *d, double *r, long long npixeff, double & var_0, double &var_n, double &delta_0, double &delta_n, int &iter, double *Mptot){

	FILE* fp;
	string file = tmp_dir + "data_sanePic.bin";
	size_t len=0;
	int idupl;

	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return 1;
	}
	len+=fread(&idupl,sizeof(int),1,fp);
	len+=fread(&var_0,sizeof(double),1,fp);
	len+=fread(&var_n,sizeof(double),1,fp);
	len+=fread(&delta_0,sizeof(double),1,fp);
	len+=fread(&delta_n,sizeof(double),1,fp);
	len+=fread(&iter,sizeof(int),1,fp);
	len+=fread(d,sizeof(double),npixeff,fp);
	len+=fread(r,sizeof(double),npixeff,fp);
	len+=fread(S,sizeof(double),npixeff,fp);
	len+=fread(Mptot,sizeof(double),npixeff,fp);

	iter++;

	fclose(fp);

	return 0;

}


int write_disk(string tmp_dir, double *d, double *r, double *S, long long npixeff, double var_0, double var_n, double delta_0, double delta_n, int iter, int idupl, double *Mptot){

	FILE *fp;
	string file = tmp_dir + "data_sanePic.bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "wb")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return 1;
	}
	len+=fwrite(&idupl,sizeof(int),1,fp);
	len+=fwrite(&var_0,sizeof(double),1,fp);
	len+=fwrite(&var_n,sizeof(double),1,fp);
	len+=fwrite(&delta_0,sizeof(double),1,fp);
	len+=fwrite(&delta_n,sizeof(double),1,fp);
	len+=fwrite(&iter,sizeof(int),1,fp);
	len+=fwrite(d,sizeof(double),npixeff,fp);
	len+=fwrite(r,sizeof(double),npixeff,fp);
	len+=fwrite(S,sizeof(double),npixeff,fp);
	len+=fwrite(Mptot,sizeof(double),npixeff,fp);

	fclose(fp);

	return 0;

}


int restore_session(string tmp_dir, string filename, int &completed_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns)
{

	FILE* fp;
	string file = tmp_dir + "data_saved_sanePS_" + filename + ".bi";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return 0;
	}

	len+=fread(&completed_step,sizeof(int),1,fp);

	switch(completed_step)
	{
	case 2:
		// load commonm2
		for (int ii = 0; ii < ncomp; ii++)
			for(long in=0; in<ns; in++)
				len+=fread(&commonm2[ii][in], sizeof(double), 1, fp);
		break;

	case 3:
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				//				fprintf(fp,"%10.15g\n",N[idet][ibin]);
				len+=fread(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellth[idet][ibin], sizeof(double), 1, fp);
		break;

	case 4:
	case 5:
		// load N, P, Rellexp, Rellth, SPref
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellth[idet][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellexp[idet][ibin], sizeof(double), 1, fp);

		len+=fread(SPref, sizeof(double), nbins, fp);
		break;

	case 6:
		cout << "completed_step is equal to 6, sanePS has already computed : " << file << ". Continue ...\n" << endl;
		return 0;

	default :
		cout << "completed_step has an incorrect value in : " << file << ". Exiting ...\n" << endl;
		return 1;

	}

	fclose(fp);

	return 0;
}


int save_session(string tmp_dir, string filename, int completed_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns){

	FILE* fp;
	string file = tmp_dir + "data_saved_sanePS_" + filename + ".bi"; // TODO : changer en .bin
	size_t len=0;

	cout << "opening " << file << " for writing\n";

	if (NULL == (fp = fopen(file.c_str(), "w+")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return 1;
	}

	len+=fwrite(&completed_step,sizeof(int),1,fp);

	cout << "completed_step : " << completed_step << endl;

	switch(completed_step)
	{
	case 2:
		// load commonm2
		for (int ii = 0; ii < ncomp; ii++)
			for(long in=0; in<ns; in++)
				len+=fwrite(&commonm2[ii][in], sizeof(double), 1, fp);
		break;

	case 3:
		// write N, P, Rellth
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellth[idet][ibin], sizeof(double), 1, fp);
		break;

	case 4:
	case 5:
		// write N, P, Rellexp, SPref
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellth[idet][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellexp[idet][ibin], sizeof(double), 1, fp);

		len+=fwrite(SPref, sizeof(double), nbins, fp);
		break;

	case 6:
		return 0;

	default :
		cout << "completed_step has an incorrect value in : " << file << ". Exiting ...\n" << endl;
		return 1;

	}

	fclose(fp);

	return 0;
}
