
#ifndef WRITE_MAPS_TO_DISK_H_
#define WRITE_MAPS_TO_DISK_H_

void write_maps_to_disk(double *S, long NAXIS1, long NAXIS2, std::string outdir, long long *indpix, long long *indpsrc,
		double *Mptot, long long addnpix, long long npixsrc, int factdupl, long ntotscan, struct wcsprm *wcs);

#endif /* WRITE_MAPS_TO_DISK_H_ */
