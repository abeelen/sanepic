
#ifndef TOOLS_H_
#define TOOLS_H_

#include "struct_definition.h"

void copy_ref_pos(fitsfile * fptr, fitsfile *outfptr, std::string name, long min_sample, long max_sample);
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);
void copy_channels(fitsfile * fptr, fitsfile *outfptr);
void copy_time(fitsfile * fptr, fitsfile *outfptr, double *time, long min_sample, long max_sample);
void copy_signal(fitsfile * fptr, fitsfile *outfptr,   std::string name, long min_sample, long max_sample, struct detectors det);
void copy_mask(fitsfile * fptr, fitsfile *outfptr,   std::string name, long min_sample, long max_sample, struct detectors det);
void copy_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long min_sample, long max_sample, struct detectors det);

#endif /* TOOLS_H_ */
