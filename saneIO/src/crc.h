#ifndef CRC_H_
#define CRC_H_

#include <stdlib.h>
#include "struct_definition.h"

unsigned checksum(void *buffer, size_t len, unsigned int seed);
void compute_checksum(std::string ini_file, std::string tmp_dir, long long npix, long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk);
int write_checksum(std::string tmp_dir, struct checksum chk);
void read_checksum(std::string tmp_dir, struct checksum &chk);
bool compare_checksum(struct checksum chk_t, struct checksum chk_t2);
void load_idupl(std::string tmp_dir, std::string out_dir, int &idupl);
void load_from_disk(std::string tmp_dir, std::string out_dir, double *S, double *d, double *r, long long npixeff, double &var_0, double &var_n, double &delta_0, double &delta_n, int &iter, double *Mptot, double *PtNPmatStot);
void write_disk(std::string tmp_dir, double *d, double *r, double *S, long long npixeff, double var_0, double var_n, double delta_0, double delta_n, int iter, int idupl, double *Mptot, double *PtNPmatStot);

#endif /* CRC_H_ */
