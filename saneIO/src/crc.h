#ifndef CRC_H_
#define CRC_H_

#include <stdlib.h>
#include "struct_definition.h"

/*! Given a buffer, which size is len, return the checksum of this buffer, added to seed */
unsigned checksum(void *buffer, size_t len, unsigned int seed);

void compute_checksum(struct param_common dir, struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_saneInv inv_param, struct param_sanePS ps_param, struct param_sanePic pic_param, struct samples samples_struct, long long npix,
		long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk);

/*! Given a checksum structure, write this struct to disk in a binary file */
int write_checksum(std::string tmp_dir, struct checksum chk, std::string projectname);

/*! Given a checksum binary file, read from disk and fill the struct checksum */
void read_checksum(std::string tmp_dir, struct checksum &chk, std::string projectname);

/*! Given 2 filled checksum struct, compare each member of both struct and return True if they are the same, False otherwise */
bool compare_checksum(struct checksum chk_t, struct checksum chk_t2);

int load_idupl(std::string tmp_dir, std::string out_dir, int &idupl);

int load_from_disk(std::string tmp_dir, std::string out_dir, double *S, double *d, double *r, long long npixeff, double &var_0, double &var_n, double &delta_0, double &delta_n, int &iter, double *Mptot);

int write_disk(std::string tmp_dir, double *d, double *r, double *S, long long npixeff, double var_0, double var_n, double delta_0, double delta_n, int iter, int idupl, double *Mptot);

int restore_session(std::string tmp_dir, std::string filename, int &completed_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns);

int save_session(std::string tmp_dir, std::string filename, int completed_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns);

#endif /* CRC_H_ */
