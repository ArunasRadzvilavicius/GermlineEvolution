#ifndef GUARD_ORGANISM_H
#define GUARD_ORGANISM_H
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Gamete.h"
using std::vector;
typedef vector<int> cell;
typedef vector<int> gtype;
class Organism
{
public:
	Organism();
	Organism(int mx,int Mx,gtype WZgx,gtype Gggx,int Gtx,int Gsx,gsl_rng * rngox,double musx,double mugx,int Ggdomx,int Lx,gtype Qgx,int Msx);
	vector< vector<int> > soma;
	vector< Gamete > gline;
	gtype WZg;
	gtype Qg;
	int L;
	gtype Ggg;
	int GgP;
	int Gt;
	int Gs;
	int Ms;
	double mus;
	double mug;
	void Grow();
	void NoDivMitosis();
	void Mitosis();
	void TissueMitosis();
	void Mutation(double mux);
	void TissueMutation();
	void Gametogenesis();
	void SomaticGametogenesis();
	vector<vector<cell> > Tissues;
	gsl_rng * rngo;
	double AdultFitness;
	double TissueFitness(vector<cell> tissuex);
	vector<int> Mix(int,int);
	vector<int> Mix4(int,int,int,int);
	int Amplify(int mx, int Qx);
};
#endif
