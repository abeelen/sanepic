#ifndef WRITE_MAPS_TO_DISK_H_
#define WRITE_MAPS_TO_DISK_H_






//! Write final sanePic map to disk
/*!
 * Computes needed informations such as Coverage
 \param S Map signal
 \param NAXIS1 Number of horizontal pixels
 \param NAXIS2 Number of vertical pixels
 \param npix number of pixels
 \param dir The param_common structure
 \param indpix The pixels indices table
 \param indpsrc The masked pixels indices table
 \param Mptot CG preconditioner : Error map
 \param addnpix Number of pix to add to compute the final maps in case of duplication and/or user binary mask
 \param npixsrc Number of filled pixels in the mask
 \param factdupl Map duplication factor
 \param ntotscan Total number of input scans
 \param pos_param The param_sanePos structure
 \param proc_param The param_saneProc structure
 \param samples_struct The samples structure
 \param fcut A vector containing Noise filter cut-off frequency for each scan
 \param wcs A pointer to a wcsprm struct, that contains informations about map's projection and coordinates
 \param maskfile A fits file in which a binary mask of the map is given
 \param structPS The param_sanePS structure
 \param sanePic_struct The param_sanePic structure
 \param saneInv_struct The param_saneInv structure
 \param key A vector containing input fits METADATA keys
 \param datatype A vector containing input fits METADATA datatypes
 \param val A vector containing input fits METADATA values
 \param com A vector containing input fits METADATA commentaries
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */

int writeMapsToFits(string fname, double *S, double *Mptot, long long addnpix,
		long NAXIS1, long NAXIS2, long long *indpix, long long *indpsrc, long long npixsrc, int factdupl, long ntotscan,
		struct wcsprm *wcs, char * subheader, int nsubkeys, bool extend);


int exportExtraToFits(string fname, struct param_common dir, struct samples samples_struct,
		struct param_saneProc Proc_param, struct param_sanePos Pos_param, struct param_sanePS PS_param,
		struct param_sanePic Pic_param, struct param_saneInv Inv_param);

#endif /* WRITE_MAPS_TO_DISK_H_ */
