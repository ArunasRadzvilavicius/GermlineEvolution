#ifndef GUARD_GAMETE_H
#define GUARD_GAMETE_H
class Gamete
{
public:
	Gamete();
	Gamete(int mx,int Mx,int nAx, int WZx, int Ggx,int Gtx,int Gsx,int Qx);
	int m;
	int M;
	int nA;
	int WZ;
	int Gg;
	int Gt;
	int Gs;
	int Q;
	void PrintGam();
};
#endif
