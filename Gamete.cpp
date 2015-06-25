#include "Organism.h"
#include <assert.h>
#include <iostream>
using namespace std;

Gamete::Gamete(){};

Gamete::Gamete(int mx,int Mx,int nAx,int WZx,int Ggx,int Gtx,int Gsx,int Qx){
	m = mx;
	M = Mx;
	WZ = WZx;
	nA = nAx;
	Gg = Ggx;
	Gt = Gtx;
	Gs = Gsx;
	Q = Qx;
}

void Gamete::PrintGam(){
	cout<<" Gamete "<< m<<'-'<<M<<" WZ="<<WZ<<" Gg="<<Gg;
}
