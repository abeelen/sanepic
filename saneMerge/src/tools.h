

#ifndef TOOLS_H_
#define TOOLS_H_

void file_compatibility_verification(struct samples samples_struct);

void copy_ref_pos(fitsfile *outfptr, struct samples samples_struct, long ns_final);
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);
void copy_channels(fitsfile * fptr, fitsfile *outfptr);
void copy_time(fitsfile *outfptr, struct samples samples_struct, long ns_final);
void copy_signal(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final);
void copy_mask(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final);
void copy_RA_DEC(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final);

#endif /* TOOLS_H_ */
