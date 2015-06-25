#include <iostream>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include "Moose.h"
using namespace std;

int main(int argc, char **argv)
{	
	if (argc!=12) {
		cout << "delta Ggx Gtx Gsx Msx musx mugx Gginvx IMx Ggdomx Lx Qx" << endl;
		cout << "IM=0 (UPI) or 1 (BPI)"<<endl;
		return 1;
	}
	int Ggx = atoi(argv[1]);
	int Gtx = atoi(argv[2]);
	int Gsx = atoi(argv[3]);
	int Msx = atoi(argv[4]);
	double musx = atof(argv[5]);
	double mugx = atof(argv[6]);
	int Gginvx = atoi(argv[7]);
	int IMx = atoi(argv[8]);
	int Ggdomx = atoi(argv[9]);
	int Lx = atoi(argv[10]);
	int Qx = atoi(argv[11]);
	assert(Ggdomx==Ggx or Ggdomx==Gginvx);
	Moose moosy;
	int ntot = 0;
	int nfixed = 0;
	for (int i=0;i<10;i++){
		moosy = Moose(100,Ggx,Gtx,Gsx,Msx,musx,mugx,Gginvx,IMx,Ggdomx,Lx,Qx);
    		nfixed += moosy.retval;
		ntot += 1;
    	}
	cout<<ntot<<' '<<nfixed<<endl;
	return 0;
}
