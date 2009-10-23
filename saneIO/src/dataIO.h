#ifndef DATAIO_H_
#define DATAIO_H_

#include <string>
#include <vector>

using namespace std;

// int read_data_std(string fname, int frame, int fs, int ns, void* data, string field, char type); // on garde
// int read_data(string fname, int frame, int fs, int ns, void* data, string field, char type); // on garde


void read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets);
void read_flpoint_from_fits(string filename, short *FLAG);
void read_flag_from_fits(string filename, string field, short *& mask, long &ns);
void read_signal_from_fits(string filename, double *signal, string field);
void read_ReferencePosition_from_fits(string filename, double *&RA, double *&DEC, double *&PHI, short *&FLAG, long &ns);
void read_ra_from_fits(string filename, string field, double *& ra, long & ns);
void read_dec_from_fits(string filename, string field, double *& dec, long & ns);


#endif /* DATAIO_H_ */
