#include <string>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "dataIO.h"

using namespace std;

//TODO: remove C++ dependancies (string fname)
//TODO: inline/extern inline if needed (not likely)

#include "dataIO.h"

// #include "getdata.h"

using namespace std;


int read_data_std(string fname, int frame, int fs, int ns,
		void* data, string field, char type)
{

	int sizetype;
	char test[2];
	size_t result;

	test[0] = type;
	test[1] = '\0';
	string typestr = string(test);
	//  printf("type = %s\n",test);
	free(test);

	FILE *fp;

	if (typestr == "d") sizetype = 8;
	if (typestr == "c") sizetype = 1;

	string filename = fname + field;

	if((fp = fopen(filename.c_str(),"r"))!=NULL){
		fseek(fp,(20*frame+fs)*sizetype,SEEK_SET);
		result = fread(data,sizetype,ns,fp);
		fclose(fp);
	}else{
		printf("Error. Could not open %s. Exiting...\n",filename.c_str());
		exit(0);
	}

		return 1;
}




// int read_data(string fname, int frame, int fs, int ns,
//               void* data, string field, char type)
// {
//   int error_code, nread;

//   char ffname[100];
//   strcpy(ffname,fname.c_str());

//   nread = GetData(ffname, field.c_str(),
//                     frame,fs, /* 1st sframe, 1st samp */
//                     0, ns, /* num sframes, num samps */
//                     type, data,
//                     &error_code);

//   if (error_code != GD_E_OK) {
//     cerr << "    GetData Error while reading "<< field
// 	 << " from " << fname <<":\n";
//     cerr << GD_ERROR_CODES[error_code] << "\n";
//     cerr << " Frame: " << frame << "\n";

//     //exit(0);
//   }

//   if(nread == 0) {
//     cerr << "Warning: nread = 0\n";
//     //exit(0);
//   }

//   return nread;
// }
