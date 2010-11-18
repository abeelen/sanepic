#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <sysexits.h>
#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "Corr_preprocess.h"
#include "dataIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "todprocess.h"
#include "map_making.h"
#include "struct_definition.h"
#include "covMatrix_IO.h"

#include <time.h>

#include <gsl/gsl_math.h>
#include <fftw3.h>
#include "mpi_architecture_builder.h"
#include "struct_definition.h"


#include "dataIO.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "parser_functions.h"



using namespace std;

int test_map_mak(struct samples samples_struct, long ns,struct param_sanePre proc_param,
		std::string outdir,	std::vector<std::string> det,long ndet, double f_lppix, long iframe, int para_bolo_indice, int para_bolo_size);

int write_to_fits_data_lp(std::string fits_name, double *data_lp, string outdir,string field, long ns_total);
int copy_fits(std::string fits_name,  std::string outdir, long ns_total);

template <class T>
void insert_array_in_image(fitsfile *fptr, fitsfile *outfptr, string field, T *array, long ns_total)
/*! insert a mask row in the output mask image */
{

	int status=0; // fits error status


	long fpixel[2]={1,1};
	long rowIndex = find_channel_index(fptr, field.c_str());
	fpixel[1] = rowIndex;

	string type = (typeid (T).name());

	if(type=="i")
		fits_write_pix(outfptr, TINT, fpixel, ns_total, array, &status);
	else
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_total, array, &status);
}

