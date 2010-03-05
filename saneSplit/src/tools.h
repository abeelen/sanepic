
#ifndef TOOLS_H_
#define TOOLS_H_

#include "struct_definition.h"

//void read_Split_file(std::string fname, std::vector< long > &cut_sample, struct samples sample_struct);
void copy_ref_pos(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_final);
//void read_ReferencePosition_from_pointer(fitsfile * fptr, double *&RA, double *&DEC, double *&PHI, long &ns);
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);
void copy_channels(fitsfile * fptr, fitsfile *outfptr);
void copy_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final);
void copy_signal(fitsfile * fptr, fitsfile *outfptr,   std::string name, long ns_final, struct detectors det);
void copy_mask(fitsfile * fptr, fitsfile *outfptr,   std::string name, long ns_final, struct detectors det);
void copy_RA(fitsfile * fptr, fitsfile *outfptr, string name, long ns_final, struct detectors det);
void copy_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long ns_final, struct detectors det);

#endif /* TOOLS_H_ */
