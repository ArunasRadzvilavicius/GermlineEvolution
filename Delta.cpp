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
	int Ggx = atoi(argv[1]);	// age of germline sequestration (in cell divs), 1...10
	int Gtx = atoi(argv[2]);	// age of tissue fate acquisition (in cell divs), 1...10
	int Gsx = atoi(argv[3]);	// tissue size (in cell divs)
	int Msx = atoi(argv[4]);	// number of mitochondria per somatic cell, 50...200
	double musx = atof(argv[5]);	// mtDNA replication error rate
	double mugx = atof(argv[6]);	// background mutation rate
	int Gginvx = atoi(argv[7]);	// age germline seq determined by the invading allele
	int IMx = atoi(argv[8]);	// mitochondrial inheritance: 0=UPI, 1=BPI
	int Ggdomx = atoi(argv[9]);	// dominant age of germline seq, either Ggx or Gginvx
	int Lx = atoi(argv[10]);	// total lifespan L=Gtx+Gsx!
	int Qx = atoi(argv[11]);	// log egg size in oogamy 0...Ggx
	assert(Ggdomx==Ggx or Ggdomx==Gginvx);
	Moose moosy;
	int ntot = 0;
	int nfixed = 0;
	// Calculate the fixation probability of the invading germline allele
	for (int i=0;i<1;i++){
		moosy = Moose(500,Ggx,Gtx,Gsx,Msx,musx,mugx,Gginvx,IMx,Ggdomx,Lx,Qx);
    		nfixed += moosy.retval;
		ntot += 1;
    	}
	cout<<ntot<<' '<<nfixed<<endl;
	return 0;
}
