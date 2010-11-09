#ifndef DATAIO_H_
#define DATAIO_H_

#include <string>
#include <vector>
#include "struct_definition.h"

extern "C" {
#include <fitsio.h>
}

using namespace std;

int read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets);
int read_ReferencePosition_from_fits(string filename, double *&RA, double *&DEC, double *&PHI, long &ns);

int read_ra_dec_from_fits(string filename, string field, double *&ra, double *& dec, long & ns);
int read_signal_from_fits(string filename, string field, double *& signal, long & ns);
int read_flag_from_fits(string filename, string field, int *&mask, long & ns);
int read_image_2D_from_fits(string filename, double *&image, string hdu_name, long & ns, long & ndet);

int read_bolo_list(string fname, std::vector<string> &boloname, long &ndet);
int read_channels(fitsfile *fptr, char **& data, long &nBolos);
long find_channel_index(fitsfile *fptr, const char * field);

int read_time_from_fits(string filename, double *& time, long ns);
int test_format(string fitsname);

/*! copy offsets table from this file to output file */
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);

/*! copy channels list from this file to output file */
void copy_channels(fitsfile * fptr, fitsfile *outfptr);

#endif /* DATAIO_H_ */
