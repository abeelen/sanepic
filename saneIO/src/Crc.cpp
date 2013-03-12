#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iostream>
#include <string>

#include <stdint.h>

#include "StructDefinition.h"
#include "ImageIO.h"

#include "Crc.h"
#include "Utilities.h"

extern "C" {
#include "wcslib/wcshdr.h"
}

using namespace std;

uint32_t checksum(void *buffer, size_t len, uint32_t seed){

	// transform the buffer pointer to a char* pointer
	char *buf = (char *)buffer;

	for (size_t ii = 0; ii < len; ++ii)
		seed += (uint32_t)(*buf++);

	return seed;
}


void compute_checksum(struct param_common dir, struct param_sanePos pos_param, struct param_saneProc proc_param,
		struct param_saneInv inv_param, struct param_sanePS ps_param, struct param_sanePic pic_param, struct samples samples_struct, long long npix,
		long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk)
{

	FILE *fp;
	size_t len;
	string file;
	char *buf_ini;

	/*------------------------------------------------------------------*/

	// fill a buffer with all the struct values !
	string buffer_crc = dir.input_dir + dir.output_dir + dir.tmp_dir + dir.data_dir + dir.fits_filelist + dir.bolo_suffix +
			dir.bolo +
			StringOf(pos_param.pixdeg) + StringOf(pos_param.flgdupl) + pos_param.maskfile + StringOf(pos_param.projgaps) +
			pos_param.projcode + pos_param.axistype + StringOf(pos_param.lon) + StringOf(pos_param.lat) + StringOf(pos_param.eq2gal) + StringOf(pos_param.gal2eq) +
			StringOf(proc_param.CORRon) + StringOf(proc_param.fill_gap) + StringOf(proc_param.remove_linear) + StringOf(proc_param.fhp) +
			proc_param.fcut_file + StringOf(proc_param.fsamp) + StringOf(proc_param.napod) + StringOf(proc_param.poly_order) +
			StringOf(proc_param.remove_polynomia) +
			inv_param.cov_matrix + inv_param.cov_matrix_suffix + inv_param.noise_dir +
			ps_param.cov_matrix + ps_param.cov_matrix_suffix + ps_param.ell + ps_param.ell_suffix + ps_param.mix +
			ps_param.mix_suffix + StringOf(ps_param.ncomp) + ps_param.signame +
			StringOf(pic_param.itermax) + StringOf(pic_param.iterw) + pic_param.map_prefix;

	buf_ini = (char*)buffer_crc.c_str();

	chk.chk_ini_file=checksum(buf_ini, (long)buffer_crc.size(), 0);
//	cout << "ini chk : " << chk.chk_ini_file<< endl;
	/*------------------------------------------------------------------*/

	char buf[19 * (34 + 48)];
	// 19 lines * ( 34 char for key and value (and spaces) + 48 for the longest comment, but you can extend to 80 )

	file= dir.tmp_dir + "mapHeader.keyrec";
	if (NULL == (fp = fopen(file.c_str(), "r")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return;
	}
	len = fread(buf, sizeof(char), sizeof(buf), fp);
	fclose(fp);

	chk.chk_wcs_file=checksum(buf, len, 0);
	//	printf("The checksum of %s is %u\n", file.c_str(), chk.chk_wcs_file);

	/*------------------------------------------------------------------*/

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

int load_idupl(string tmp_dir, int &idupl){
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

int load_from_disk(string tmp_dir, double *S, double *d, double *r, long long npixeff, double & var_0, double &var_n, double &delta_0, double &delta_n, int &iter, double *Mptot){

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
