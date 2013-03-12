#ifndef CRC_H_
#define CRC_H_

#include <stdint.h>
#include <stdlib.h>
#include "StructDefinition.h"

//! Given a buffer, which size is len, computes the checksum of this buffer, added to seed
/*!
  \param buffer A buffer containing characters which checksum needs to be computed
  \param len The length of the "buffer"
  \param seed A seed that can be add to the checksum (not used)
  \return The checksum value (as an unsigned int)
 */
uint32_t checksum(void *buffer, size_t len, uint32_t seed);

//! Computes the checksum of whole structure, indpix and indpsrc table, stacking them into a buffer
/*!
 \param dir The param_common structure
 \param samples_struct The samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_saneProc structure
 \param ps_param The param_sanePS structure
 \param pic_param The param_sanePic structure
 \param inv_param The param_saneInv structure
 \param chk The checksum structure that has to be filled with all the inputs
 \param npix Number of pixels in the map
 \param indpix The pixels indices table
 \param indpsrc The masked pixels indices table
 \param indpsrc_size indpsrc Size
 */
void compute_checksum(struct param_common dir, struct param_sanePos pos_param, struct param_saneProc proc_param,
		struct param_saneInv inv_param, struct param_sanePS ps_param, struct param_sanePic pic_param, struct samples samples_struct, long long npix,
		long long* indpix, long long* indpsrc, long long indpsrc_size, struct checksum &chk);

//! Given a checksum structure, write this struct to disk in a binary file
/*!
  \param tmp_dir The temporary files directory
  \param chk The checksum structure that has to be wrote on disk
  \param projectname The name of the project calling this function (sanePic or sanePS)
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_checksum(std::string tmp_dir, struct checksum chk, std::string projectname);

//! Given a checksum binary file, read from disk and fill the struct checksum */
/*!
  \param tmp_dir The temporary files directory
  \param chk The checksum structure that has to be filled with the file
  \param projectname The name of the project calling this function (sanePic or sanePS)
 */
void read_checksum(std::string tmp_dir, struct checksum &chk, std::string projectname);

//! Given 2 filled checksum struct, compare each member of both struct and return True if they are the same, False otherwise
/*!
  \param chk_t The checksum structure that has to be compared to chk_t2
  \param chk_t2 The checksum structure that has to be compared to chk_t
  \return True if chk_t and chk_t2 are different, False if they are the same
 */
bool compare_checksum(struct checksum chk_t, struct checksum chk_t2);

//! Load idupl from "data_sanePic.bin" file (used to restore a lost sanePic session)
/*!
 \param tmp_dir The temporary files directory
 \param idupl An integer : Duplication map iteration marker
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int load_idupl(std::string tmp_dir, int &idupl);

//! Load the informations stored in "data_sanePic.bin" file to restore a lost sanePic session
/*!
 \param tmp_dir The temporary files directory
 \param S Map signal
 \param d
 \param r
 \param npixeff number of pixel (effective value)
 \param var_0 initial conjuguent criteria
 \param var_n actual loop conjuguent gradient criteria
 \param delta_0 initial conjuguent criteria
 \param delta_n actual loop conjuguent gradient criteria
 \param iter The number of iterations that has already been done
 \param Mptot The CG preconditioner
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int load_from_disk(std::string tmp_dir, double *S, double *d, double *r, long long npixeff, double &var_0, double &var_n, double &delta_0, double &delta_n, int &iter, double *Mptot);

//! Save a sanePic running session by storing informations on disk in "data_sanePic.bin" file
/*!
 \param tmp_dir The temporary files directory
 \param idupl An integer : Duplication map iteration marker
 \param S Map signal
 \param d
 \param r
 \param npixeff number of pixel (effective value)
 \param var_0 initial conjuguent gradient convergence criteria
 \param var_n actual loop conjuguent gradient convergence criteria
 \param delta_0 initial conjuguent gradient convergence criteria
 \param delta_n actual loop conjuguent gradient convergence criteria
 \param iter The number of iterations that has already been done
 \param Mptot The CG preconditioner
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_disk(std::string tmp_dir, double *d, double *r, double *S, long long npixeff, double var_0, double var_n, double delta_0, double delta_n, int iter, int idupl, double *Mptot);

#endif /* CRC_H_ */
