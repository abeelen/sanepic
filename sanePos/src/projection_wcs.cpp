/*
 * projection_wcs.cpp
 *
 *  Created on: 7 août 2009
 *      Author: matthieu
 */




#include "projection_wcs.h"

#define NELEM 2

using namespace std;


// TODO : remove, this has been replaced by computeMapHeader()
void fits_header_generation(string outdir, const char *fits_file,double pixdeg, bool default_projection,double *tanpix,double *tancoord)
{

	FILE *fout;
	fitsfile *fptr;
	int i, status = 0;
	struct wcsprm *wcs;
	int nkeyrec;
	char *header, *hptr;


//TODO : CRVAL different entre header et header2 (optimap_flux), why ??
	if(default_projection){

		//int nkeyrec;

		//#define NELEM 2
		//int i, status = 0;
		//double pixcrd[10][NELEM];
		//double imgcrd[10][NELEM],phi[10],theta[10],world[10][NELEM];

		//int stat[NWCSFIX];
		//char *header, *hptr;

		//struct wcsprm *wcs;

		/* The following routine simulates the actions of a FITS header parser. */
		wcs = (wcsprm *)malloc(sizeof(struct wcsprm));
		wcs->flag = -1;

		//	struct wcsprm *wcs2;




		/* List status return messages. */
		/*printf("\nList of wcs status return values:\n");
		for (status = 1; status <= 13; status++) {
			printf("%4d: %s.\n", status, wcs_errmsg[status]);
		}*/

		//		printf("before\n");

		parser(wcs,tanpix,tancoord, pixdeg);
		//printf("after\n");


		/* print the header */

		//	wcsprt(wcs);

		/* print the header */

		if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
			printf("%4d: %s.\n", status, wcs_errmsg[status]);
			exit(0);
		}

		hptr = header;
		printf("\n\n Projection Header :\n");
		outdir=outdir + "projection_header.keyrec";
		fout = fopen(outdir.c_str(),"w");
		for (i = 0; i < nkeyrec; i++, hptr += 80) {
			printf("%.80s\n", hptr);
			fprintf(fout,"%.80s\n", hptr);
		}
		fclose(fout);
		free(header);

		//printf("good\n");

	}else{


		double pixcrd[10][NELEM];
		double imgcrd[10][NELEM],phi[10],theta[10],world[10][NELEM];
		int stat[NWCSFIX];
		int nreject, nwcs;

		fits_open_file(&fptr, fits_file, READONLY, &status);

		if ((status = fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec,
				&status))) {
			fits_report_error(stderr, status);
			exit(0);
		}

		/*-----------------------------------------------------------------------*/
		/* Basic steps required to interpret a FITS WCS header, including -TAB.  */
		/*-----------------------------------------------------------------------*/

		/* Parse the primary header of the FITS file. */
		if ((status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs,
				&wcs))) {
			fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
		}

		/* Read coordinate arrays from the binary table extension. */
		if ((status = fits_read_wcstab(fptr, wcs->nwtb, (wtbarr *)wcs->wtb,
				&status))) {
			fits_report_error(stderr, status);
			exit(0);
		}

		/* Translate non-standard WCS keyvalues. */
		if ((status = wcsfix(7, 0, wcs, stat))) {
			for (i = 0; i < NWCSFIX; i++) {
				if (stat[i] > 0) {
					fprintf(stderr, "wcsfix ERROR %d: %s.\n", status,
							wcsfix_errmsg[stat[i]]);
				}
			}

			exit(0);
		}


		//printf the header
		//wcsprt(wcs);

		/* print the header */
		if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
			exit(0);
		}

		hptr = header;
		printf("\n\n Projection Header :\n");
		outdir=outdir + "projection_header2.keyrec";
		fout = fopen(outdir.c_str(),"w");
		for (i = 0; i < nkeyrec; i++, hptr += 80) {
			printf("%.80s\n", hptr);
			fprintf(fout,"%.80s\n", hptr);
		}
		fclose(fout);
		free(header);

		//printf("good2b\n");


		//printf("Testing %s;\n",wcs->cel.prj.code);

		/*for (i=0; i < 10; i++){
			pixcrd[i][wcs->lng] = (double) i;
			pixcrd[i][wcs->lat] = (double) i;
		}*/


		if ((status = wcsp2s(wcs, 10, NELEM, pixcrd[0], imgcrd[0], phi, theta, world[0], stat))) {
			printf("   wcsp2s(1) ERROR %2d (lat1 = %f)\n", status, phi[1] );
		}
		//printf("good2bb\n");

		//wcsprt(wcs);

		/*for (i = 0; i < 10; i++) {

			printf("  (%5.1f,%6.1f) -> (%10.6f,%10.6f)",
					pixcrd[i][0], pixcrd[i][1],
					world[i][0],  world[i][1]);

			if (stat[i]) printf("  (BAD)");
			printf("\n");
		}*/


		wcsfree(wcs);

		free(wcs);


	}


}
// TODO: what is the parser doing here ?

void parser(struct wcsprm *wcs,double *tanpix,double *tancoord, double pixdeg)

//struct wcsprm *wcs;

{
	int i, j, status;
	double *pcij;

	const int NAXIS = 2; // number of dimensions
	double CRPIX[2] =  {tanpix[0], tanpix[0]}; // reference pixel coordinates
	const double PC[2][2] = {{    1.0,  0.0}, // linear transformation matrix
			{    0.0,  1.0}};
	double CDELT[2] =  {-pixdeg, pixdeg}; // -pixdeg, pixdeg (coordinate scales)

	char CUNIT[2][18] = {"deg", "deg"}; // en degre // units of CRVAL and CDELT
	char CTYPE[2][18] = {"RA---TAN", "DEC--TAN"}; // X, Y, projection type
	char CNAME[2][18] = {"Right Ascension", "Declination"};

	double CRVAL[2] = {tancoord[0], tancoord[1]}; // ra mean, dec mean
	// => Celestial longitude and latitude of the fiducial point

	const double LONPOLE  = 150.0; // reference pole longitude
	const double LATPOLE  = 150.0; // reference pole latitude

	int NPV = 0;
	struct pvcard PV[2]; // native latitude, longitude of the fiducial point /* Projection parameters are set in main(). */



	//printf("inside : %i\n",NAXIS);

	/* In practice a parser would read the FITS header until it encountered  */
	/* the NAXIS keyword which must occur near the start, before any of the  */
	/* WCS keywords.  It would then use wcsini() to allocate memory for      */
	/* arrays in the wcsprm struct and set default values.  In this          */
	/* simulation the header keyvalues are set as global variables.          */
	wcsnpv(0);
	wcsnps(0);
	wcsini(1, NAXIS, wcs);

	//printf("inside : %i\n",NAXIS);

	/* Now the parser scans the FITS header, identifying WCS keywords and    */
	/* loading their values into the appropriate elements of the wcsprm      */
	/* struct.                                                               */

	for (j = 0; j < NAXIS; j++) {
		wcs->crpix[j] = CRPIX[j];
	}

	pcij = wcs->pc;
	for (i = 0; i < NAXIS; i++) {
		for (j = 0; j < NAXIS; j++) {
			*(pcij++) = PC[i][j];
		}
	}

	for (i = 0; i < NAXIS; i++) {
		wcs->cdelt[i] = CDELT[i];
	}

	for (i = 0; i < NAXIS; i++) {
		strcpy(wcs->cunit[i], &CUNIT[i][0]);
		strcpy(wcs->ctype[i], &CTYPE[i][0]);
		strcpy(wcs->cname[i], &CNAME[i][0]);
	}

	for (i = 0; i < NAXIS; i++) {
		wcs->crval[i] = CRVAL[i];
	}

	wcs->lonpole = LONPOLE;
	wcs->latpole = LATPOLE;

	/*    wcs->restfrq = RESTFRQ; */
	/*    wcs->restwav = RESTWAV; */

	wcs->npv = NPV;
	for (i = 0; i < NPV; i++) {
		wcs->pv[i] = PV[i];
	}

	/* Extract information from the FITS header. */
	if ((status = wcsset(wcs))) {
		printf("wcsset ERROR%3d\n", status);
	}

	return;
}
