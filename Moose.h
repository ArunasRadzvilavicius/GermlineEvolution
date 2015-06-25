#ifndef GUARD_MOOSE_H
#define GUARD_MOOSE_H
#include <vector>
#include "Organism.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using std::vector;

class Moose
{
public:
	Moose();
	Moose(int Nx,int Ggx,int Gtx,int Gsx,int Msx,double musx,double mugx,int Gginvx,int IMx,int Ggdomx,int Lx,int Qx);
	int Gginv;
	int N;
	int IM;
	int Ggdom;
	int L;
	int Q;
	vector<Organism> Population;
	gsl_rng * rngm;
	void GrowAll();
	void Selection();
	void Fusion();
	vector<double> Sample(vector<double> list, int n);
	void PrintPop();
	vector<Organism> Fuse(Gamete MGamx, Gamete FGamx);
	double AvgFitness();
	double AvgInvFitness();
	double AvgWTFitness();
	double GetGginvFreq();
	double mus;
	double mug;
	int Ms;
	double retval;
	double retval1;
	long GenSeed();
};
#endif
