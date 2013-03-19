#ifndef DATAIO_H_
#define DATAIO_H_

#include <string>
#include <vector>
#include "StructDefinition.h"

extern "C" {
#include <fitsio.h>
}

#define EXT_SIGNAL 0x01
#define EXT_MASK   0x02
#define EXT_CHAN   0x04
#define EXT_TIME   0x08

#define EXT_LON    0x10
#define EXT_LAT    0x20
#define EXT_REFPOS 0x40
#define EXT_OFFSET 0x80

#define EXT_NECESSARY  (EXT_SIGNAL | EXT_MASK | EXT_CHAN)
#define EXT_POS        (EXT_LON | EXT_LAT | EXT_REFPOS | EXT_OFFSET)
#define HIPE_FORMAT    (EXT_LON | EXT_LAT)
#define SANEPIC_FORMAT (EXT_REFPOS | EXT_OFFSET)
#define BOTH_FORMAT    (HIPE_FORMAT | SANEPIC_FORMAT)


using namespace std;

//! Read the channels offsets in regards to reference position, for a given fits file
/*!
 * \param filename The fits filename that contain the offsets table
 * \param bolonames A vector containing the channel list for this scan
 * \param offsets A double matrices, in which are stored the offsets
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets);

//! Read the reference position (RA, DEC, PHI and ns) for a given fits file
/*!
 * \param filename The fits filename that contain the positions table
 * \param RA RA table of the reference channel positions
 * \param DEC DEC table of the reference channel positions
 * \param PHI PHI table of the reference channel positions
 * \param ns This scan number of samples
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_ReferencePosition_from_fits(string filename, double *&RA, double *&DEC, double *&PHI, long &ns);

//! Read RA/DEC tables of a given channel in a fits file (also fill ns value)
/*!
 * \param filename The input fits filename
 * \param field A string containing the name of channel which positions has to be read
 * \param ra RA table of the considered channel (field)
 * \param dec DEC table of the considered channel (field)
 * \param ns This scan number of samples
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_LON_LAT_from_fits(string filename, string field, double *&lon, double *& lat, long & ns);

//! Read the "signal" table of a given channel in a fits file (also return ns value)
/*!
 * \param filename The input fits filename
 * \param field A string containing the name of channel which signal has to be read
 * \param signal A table containing the signal for the considered "field"
 * \param ns This scan number of samples
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_signal_from_fits(string filename, string field, double *& signal, long & ns);

//! Read the flag table (mask) of a given channel in a fits file (also return ns value)
/*!
 * \param filename The input fits filename
 * \param field A string containing the name of channel which flag table has to be read
 * \param mask A table containing the flags for the considered "field"
 * \param ns This scan number of samples
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_flag_from_fits(string filename, string field, int *&mask, long & ns);

//! Read a fits file channels list
/*!
 * \param fname The input fits filename
 * \param boloname A vector containing the name of the channels
 * \param ndet The number of channels contained in boloname
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_bolo_list(string fname, std::vector<string> &boloname, long &ndet);

//! Read the channel table of an opened fits file (called by read_bolo_list and find_channel_index)
/*!
 * \param fptr A pointer to an opened fits file
 * \param data A table containing the name of the channels
 * \param nBolos The number of channels contained in data
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_channels(fitsfile *fptr, char **& data, long &nBolos);

//! Find the channel index in an open fits file, given the bolometer name
/*!
 * \param fptr A pointer to an opened fits file
 * \param field A table containing the name of the channels
 * \return The channel index
 */
long find_channel_index(fitsfile *fptr, const char * field);

//! Read time table from a given fits file
/*!
 * \param filename The input fits filename
 * \param time A table containing the scan times
 * \param ns The number of sample in this scan
 * \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_time_from_fits(string filename, double *& time, long &ns);

int getImageExtensionSize(fitsfile *fptr,  string extname, long &nx, long &ny);
int checkImageExtension(fitsfile *fptr,    string extname, long nx, long ny);
int checkBintableExtension(fitsfile *fptr, string extname, long nx);

//! Check whether the fits files have HIPE or SANEPIC format
/*!
 * \param fitsname The input fits filename
 * \return An integer : 1 if HIPE format has been found. 2 if SANEPIC. 0 if NONE of them
 */
int testExtensions(string fitsname);

//! copy offsets table from this file to output file
/*!
 * \param fptr A pointer to the input opened fits file
 * \param outfptr A pointer to the output opened fits file in which the offsets table will be copied
 */
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);

//! copy channels list from this file to output file
/*!
 * \param fptr A pointer to the input opened fits file
 * \param outfptr A pointer to the output opened fits file in which the channels table will be copied
 */
void copy_channels(fitsfile * fptr, fitsfile *outfptr);

#endif /* DATAIO_H_ */
