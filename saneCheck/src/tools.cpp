

#include "covMatrixIO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}



using namespace std;


int check_flag(string fname,struct detectors det,long ns, string outname){

	short *flag;
	short sum=0;
	FILE * fp;

	int marge = 20;
	long ii=1;
	long tt=0;
	long rr=0;


	fp=fopen(outname.c_str(),"w");
	//	flag = new short[ns];

	//cout << "ndet : "<< det.ndet << endl;

	for(int jj=0;jj<det.ndet;jj++){

		ii=1;
//		cout << "det : " << det.boloname[jj] << endl;
		sum=0;
		//		cout << "jj : " << jj << endl;
		read_flag_from_fits(fname, det.boloname[jj], flag, ns);
//		cout << "ns " << ns << endl;
		//		cout << flag[0] << " " << flag[1] << " " << flag[2] << " " << flag[3] << endl;
		for(int kk=0;kk<ns;kk++){
			if(flag[kk]==NAN)
				cout << "Warning ! there is a NaN in the flag field of " << det.boloname[jj] << endl;
		}

		while(ii<ns-1){

			if((flag[ii]==0)&&(flag[ii+1]==1)&&(flag[ii-1]==1)){
				tt=ii-1;
				rr=ii+1;

				while((tt>0)&&(flag[tt]==1)&&((ii-tt)<marge+2))
					tt--;

				while((rr<ns)&&(flag[rr]==1)&&((rr-ii)<marge+2))
					rr++;
			}

			if(((rr-ii)>=marge)&&((ii-tt)>=marge)){
				cout << ii << " " << rr-ii << " " << ii-tt << endl;
				//getchar();
				flag[ii]=1;
				cout << "singleton trouvé en " << det.boloname[jj] << " sample n° " << ii << endl;
			}

			sum+=flag[ii];
			ii++;

		}

		if(sum==ns){
			cout << "Warning ! " << det.boloname[jj] << " is totally flagged" << endl;
			fprintf(fp, "%s\n", (char*)(det.boloname[jj].c_str()));
		}else{
			if(sum>80*ns/100)
				cout << "Warning ! " << det.boloname[jj] << " is more than 80% flagged" << endl;
			else if(sum>50*ns/100)
				cout << "Warning ! " << det.boloname[jj] << " is more than 50% flagged" << endl;
		}

//		cout << "fin bolo : " << det.boloname[jj] << endl;
//		getchar();

		delete [] flag;
	}

	fclose(fp);

	return 0;
}
