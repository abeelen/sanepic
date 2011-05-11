#ifndef SANEPOS_PREPROCESS_H_
#define SANEPOS_PREPROCESS_H_

#include <string>
#include <vector>
#include "mpi_architecture_builder.h"


extern "C" {
	#include <wcslib/wcs.h>
}

//! Modify scans' flags on disk (in dirfile tree), to give them a different value if there are contained in user mask
/*!
 * Each rank treats his scans (indice_min and _max are given has inputs)
 \param tmp_dir A string containing the temporary files pathname
 \param samples_struct A samples structure
 \param indpsrc The masked pixels indices table, to be read from disk
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param iframe_min Actual rank first frame indice
 \param iframe_max Actual rank last frame indice
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int modify_mask_flag_in_dirfile(std::string tmp_dir, struct samples samples_struct, long long *indpsrc,
		long NAXIS1, long NAXIS2, long iframe_min, long iframe_max);

//! Get coordinates of pixels that are seen, computes the position to pixel projetcion matrices : One binary file per bolometer and per scan
/*!
 * This routine is used why SANEPIC format only !!
 * Each rank treats his scans (indice_min and _max are given has inputs)
 \param tmp_dir A string containing the temporary files pathname
 \param samples_struct A samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_sanePre structure
 \param iframe_min Actual rank first frame indice
 \param iframe_max Actual rank last frame indice
 \param wcs A pointer to a wcsprm struct
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param mask User binary mask of the map
 \param factdupl map duplication factor
 \param addnpix Number of pix to add to compute the final maps in case of duplication and/or user binary mask
 \param pixon This array is used to store the rules for pixels : seen or not
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \param indpsrc The masked pixels indices table, to be read from disk
 \param npixsrc number of pixels included in Crossing Constraint Removal (binary mask)
 \param flagon If >0 : some pixels are apodized or outside the map
 \param pixout If >0 : Indicates that at least one pixel has been flagged and is out
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int computePixelIndex(std::string tmpdir,
		struct samples samples_struct, struct param_sanePre proc_param, struct param_sanePos pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout, std::vector<std::vector<std::string> > bolo_vect);

//! Get coordinates of pixels that are seen, computes the position to pixel projetcion matrices : One binary file per bolometer and per scan
/*!
 * This routine is used why HIPE format only !!
 * Each rank treats his scans (indice_min and _max are given has inputs)
 \param tmp_dir A string containing the temporary files pathname
 \param samples_struct A samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_sanePre structure
 \param iframe_min Actual rank first frame indice
 \param iframe_max Actual rank last frame indice
 \param wcs A pointer to a wcsprm struct
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param mask User binary mask of the map
 \param factdupl map duplication factor
 \param addnpix Number of pix to add to compute the final maps in case of duplication and/or user binary mask
 \param pixon This array is used to store the rules for pixels : seen or not
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \param indpsrc The masked pixels indices table, to be read from disk
 \param npixsrc number of pixels included in Crossing Constraint Removal (binary mask)
 \param flagon If >0 : some pixels are apodized or outside the map
 \param pixout If >0 : Indicates that at least one pixel has been flagged and is out
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int computePixelIndex_HIPE(std::string tmpdir,
		struct samples samples_struct, struct param_sanePre proc_param, struct param_sanePos pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout, std::vector<std::vector<std::string> > bolo_vect);

#endif /* SANEPOS_PREPROCESS_H_ */
