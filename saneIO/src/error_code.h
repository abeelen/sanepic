/*
 * error_code.h
 *
 *  Created on: Oct 28, 2011
 *      Author: abeelen
 */

#ifndef ERROR_CODE_H_
#define ERROR_CODE_H_

#define OK 0x0000

/*! ini file was not found : incorrect name, or name was not given by user */
#define INI_NOT_FOUND 0x0001

/*! input file related problems */
#define DIR_PROBLEM 0x0002
#define FILE_PROBLEM 0x0003
#define FILE_SIZE_PROBLEM 0x0004

/*! data or/and input path(s) are incorrects or were not correctly filled in ini file */
#define DATA_INPUT_PATHS_PROBLEM 0x0005

/*! output path is incorrect or was not correctly filled in ini file */
#define OUPUT_PATH_PROBLEM 0x0006

/*! temporary files path is incorrect or was not correctly filled in ini file */
#define TMP_PATH_PROBLEM 0x0008

/*! channels files were not found or not correctly filled in ini file */
#define BOLOFILE_NOT_FOUND 0x0010

/*! pixel size was not correctly filled in ini file or the value is wrong */
#define PIXDEG_WRONG_VALUE 0x0020

/*! file format was not correctly filled in ini file or the value is absent */
#define FILEFORMAT_NOT_FOUND 0x0040

/*! napod was not correctly filled in ini file or the value is < 0 */
#define NAPOD_WRONG_VALUE 0x0080

/*! fsamp was not correctly filled in ini file or the value is < 0 */
#define FSAMP_PROBLEM 0x0100

/*! fhp was not correctly filled in ini file or the value is < 0 */
#define FHP_PROBLEM 0x0101

/*! fcut was not correctly filled in ini file or the value is < 0 */
#define FCUT_PROBLEM 0x0102

/*! ncomp was not correctly filled in ini file or the value is < 0 */
#define NCOMP_WRONG_VALUE 0x0400

//! Parser functions binary return flag
/*! ell file(s) were not found or the name was not correctly given in ini file */
#define ELL_FILE_NOT_FOUND 0x0800

/*! mixing matrix file(s) were not found or not correctly filled in ini file */
#define MIX_FILE_NOT_FOUND 0x1000

//! Parser functions binary return flag
/*! Noise path or covariance matrices was/were not found or not correctly filled in ini file */
#define SANEINV_INPUT_ERROR 0x2000

/*! fits filelist file was not found or not correctly filled in ini file */
#define FITS_FILELIST_NOT_FOUND 0x4000


#endif /* ERROR_CODE_H_ */
