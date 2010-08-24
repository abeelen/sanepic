#ifndef CRC_H_
#define CRC_H_

#include <stdlib.h>

#include "struct_definition.h"

unsigned checksum(void *buffer, size_t len, unsigned int seed);
void compute_checksum(std::string ini_file, std::string tmp_dir, double* Pnd, long long npix, long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk);
void write_checksum(std::string tmp_dir, struct checksum chk);
void read_checksum(std::string tmp_dir, struct checksum &chk);
bool compare_checksum(struct checksum chk_t, struct checksum chk_t2);
void load_from_disk(std::string tmp_dir, std::string out_dir, double *S, double *d, long long *indpix, long long npixeff, double &var_n, double &delta_n, int &iter);
void write_disk(std::string tmp_dir, double *d, long long npix, double var_n, double delta_n, int iter);

#endif /* CRC_H_ */
