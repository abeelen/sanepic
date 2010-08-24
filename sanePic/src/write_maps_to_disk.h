
#ifndef WRITE_MAPS_TO_DISK_H_
#define WRITE_MAPS_TO_DISK_H_

int write_maps_to_disk(double *S, long NAXIS1, long NAXIS2, string outdir, long long *indpix, long long *indpsrc,
		double *Mptot, long long addnpix, long long npixsrc, int factdupl, long ntotscan,
		struct param_process proc_param, struct param_positions pos_param, struct detectors det,
		struct samples samples_struct, std::vector<double> fcut, struct wcsprm *wcs, string maskfile);

#endif /* WRITE_MAPS_TO_DISK_H_ */
