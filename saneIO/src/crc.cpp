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


void compute_checksum(std::string ini_file, std::string tmp_dir, double* Pnd, long long npix, long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk)
{

	FILE *fp;
	size_t len;
	string file;
	char buf[6144];
	//	struct checksum chk2;


	if (NULL == (fp = fopen(ini_file.c_str(), "r")))
	{
		cout << "Unable to open " << ini_file << " for reading\n";
		return;
	}
	len = fread(buf, sizeof(char), sizeof(buf), fp);
	char buf2[len-116];
	for(long hh=0;hh<(long)(len-116);hh++)
		buf2[hh]=buf[hh];
	//	printf("%d bytes read\n", len);
	fclose(fp);

	chk.chk_ini_file=checksum(buf2, len-116, 0);
	printf("The checksum of %s is %u\n", ini_file.c_str(), chk.chk_ini_file);

	file= tmp_dir + "mapHeader.keyrec";
	if (NULL == (fp = fopen(file.c_str(), "r")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return;
	}
	len = fread(buf, sizeof(char), sizeof(buf), fp);
	//	printf("%d bytes read\n", len);
	fclose(fp);

	chk.chk_wcs_file=checksum(buf, len, 0);
	//	printf("The checksum of %s is %#x\n", file.c_str(), cksum);
	printf("The checksum of %s is %u\n", file.c_str(), chk.chk_wcs_file);

	chk.chk_pnd=checksum(Pnd, (size_t) npix, 0);
	printf("The checksum of PNd is %u\n", chk.chk_pnd);

	chk.chk_indpix=checksum(indpix, (size_t) npix, 0);
	printf("The checksum of Indpix is %u\n", chk.chk_indpix);

	chk.chk_indpsrc=checksum(indpsrc, (size_t) indpsrc_size, 0);
	printf("The checksum of Indpsrc is %u\n", chk.chk_indpsrc);


	//	file = tmp_dir + "checksum.bin";
	//	if (NULL == (fp = fopen(file.c_str(), "w")))
	//	{
	//		cout << "Unable to open " << file << " for writing\n";
	//		return;
	//	}
	//		for(int ii =0;ii<5;ii++)
	//	fwrite(chk_tabl,sizeof(unsigned int),5,fp);
	//	fprintf(fp,"%u ",chk_tabl[ii]);

	//	fclose(fp);

	//	for (int ii=0;ii<5;ii++){
	//		chk_tabl[ii]=0;
	//		printf("chk_tabl[%d] is %u\n", ii,chk_tabl[ii]);
	//	}

	//	unsigned int *tab;
	//	tab= new unsigned int[5];
	//
	//	file = tmp_dir + "checksum.bin";
	//	if (NULL == (fp = fopen(file.c_str(), "rb")))
	//	{
	//		cout << "Unable to open " << file << " for reading\n";
	//		return;
	//	}
	//	//for(int ii =0;ii<5;ii++)
	//	//fscanf(fp,"%u ",tab+ii);
	//	len=fread(tab,sizeof(unsigned int),5,fp);
	//	fclose(fp);
	//
	//	cout << "len : " << len << endl;
	//	cout << "size : " << sizeof(unsigned int) << endl;
	//
	//	for (int ii=0;ii<5;ii++)
	//		printf("tab[%d] is %u\n", ii,tab[ii]);

}


void write_checksum(std::string tmp_dir, struct checksum chk)
{

	FILE* fp;
	string file = tmp_dir + "checksum.bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "wb")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return;
	}
	len+=fwrite(&chk.chk_ini_file,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_wcs_file,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_pnd,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_indpix,sizeof(unsigned int),1,fp);
	len+=fwrite(&chk.chk_indpsrc,sizeof(unsigned int),1,fp);

	fclose(fp);
}

void read_checksum(std::string tmp_dir, struct checksum &chk)
{

	FILE* fp;
	string file = tmp_dir + "checksum.bin";
	size_t len=0;


	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return;
	}
	len+=fread(&chk.chk_ini_file,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_wcs_file,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_pnd,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_indpix,sizeof(unsigned int),1,fp);
	len+=fread(&chk.chk_indpsrc,sizeof(unsigned int),1,fp);


	//	cout << "chkini : " << chk.chk_ini_file << endl;

	//	cout << "len read : " << len << endl;

	fclose(fp);

}

bool compare_checksum(struct checksum chk_t, struct checksum chk_t2){

	if((chk_t.chk_ini_file != chk_t2.chk_ini_file) ||
			(chk_t.chk_wcs_file!=chk_t2.chk_wcs_file) ||
			(chk_t.chk_pnd!=chk_t2.chk_pnd) ||
			(chk_t.chk_indpix!=chk_t2.chk_indpix) ||
			(chk_t.chk_indpsrc!=chk_t2.chk_indpsrc))
		return 1;

	return 0;
}

void load_from_disk(string tmp_dir, string out_dir, double *S, double *d, double *r, long long *indpix, long long npixeff, double &var_n, double &delta_n, int &iter){

	FILE* fp;
	string file = tmp_dir + "data_sanePic.bin";
	size_t len=0;
	string fits_temp;
	std::ostringstream oss;
	long NAXIS1, NAXIS2;
	struct wcsprm * wcs; // TODO change this to deal with REAL wcs !!
	int nwcs = 1;

	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return;
	}

	len+=fread(&var_n,sizeof(double),1,fp);
	len+=fread(&delta_n,sizeof(double),1,fp);
	len+=fread(&iter,sizeof(int),1,fp);
	len+=fread(d,sizeof(double),npixeff,fp);
	len+=fread(r,sizeof(double),npixeff,fp);

	oss << out_dir + "optimMap_" << iter << "b.fits";
	fits_temp=oss.str();
	read_MapHeader(tmp_dir,wcs,&NAXIS1, &NAXIS2);
	read_fits_signal(fits_temp, S, indpix, NAXIS1, NAXIS2, wcs);

	iter++;



	fclose(fp);

	wcsvfree(&nwcs, &wcs);
}


void write_disk(string tmp_dir, double *d, double *r, long long npixeff, double var_n, double delta_n, int iter){

	FILE *fp;
	string file = tmp_dir + "data_sanePic.bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "wb")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return;
	}
	len+=fwrite(&var_n,sizeof(double),1,fp);
	len+=fwrite(&delta_n,sizeof(double),1,fp);
	len+=fwrite(&iter,sizeof(int),1,fp);
	len+=fwrite(d,sizeof(double),npixeff,fp);
	len+=fwrite(r,sizeof(double),npixeff,fp);

	fclose(fp);


}
