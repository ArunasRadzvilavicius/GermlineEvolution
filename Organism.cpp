#include "Organism.h"
#include <assert.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <numeric>
using namespace std;

Organism::Organism(){}

Organism::Organism(int mx,int Mx,gtype WZgx,gtype Gggx,int Gtx,int Gsx,gsl_rng * rngox,double musx,double mugx,int Ggdomx,int Lx,gtype Qgx,int Msx){
	WZg = WZgx;
	Qg = Qgx;
	L = Lx;
	Gt = Gtx;
	Gs = Gsx;
	Ggg = Gggx;
	if (Ggg[0]==Ggg[1]){GgP = Ggg[0];}
	else if (Ggg[0]==Ggdomx or Ggg[1]==Ggdomx) {GgP = Ggdomx;}
	if (GgP<0) {GgP = -1*GgP;}
	Ms = Msx;
	mus = musx;
	mug = mugx;
	int cstate[2] = {mx,Mx};
	vector<int> cell (cstate, cstate+2);
	soma.push_back(cell);
	assert(soma.size()==1);
	assert(gline.size()==0);
	rngo = rngox;
	AdultFitness = 1.0 - pow(1.0*mx/Mx,2);
}

void Organism::Grow(){
	// go through all stages of development
	assert(int(soma.size())==1);
	int Qp = 0;
	if (soma[0][1]>Ms){Qp = log2(soma[0][1]/Ms);}
	int divcount = 0;
	assert(Qp<=Gt);
	for (int i=0;i<Qp;i++){
		Mutation(mug);
		NoDivMitosis();
		divcount+=1;
	}
	assert(soma[0][1]==Ms);
	if (GgP<=Gt){
		while (divcount<GgP){
			Mutation(mug+mus);
			Mitosis();
			divcount += 1;
		}
		assert(pregline.size()==0);
		Gametogenesis();
		assert(pregline.size()==1);
	}

	while (divcount<Gt){
		Mutation(mug+mus);
		Mitosis();
		divcount+=1;
	}
	
	assert( int(soma.size() )==pow(2.0, Gt));

	// tissue development
	Tissues = vector<vector<cell> > (0);
	// initiate tissues
	vector<cell> Tissue;
	for (int l=0;l<int(soma.size());l++){
		Tissue = vector<cell> (1);
		Tissue[0] = soma[l];
		Tissues.push_back(Tissue);
	}

	if (GgP>Gt){
		//assert(1==2);
		assert(GgP<=Gt+Gs);
		while  (divcount < GgP){
			TissueMutation();
			TissueMitosis();
			divcount += 1;
		}
		SomaticGametogenesis();
		while  (divcount < Gs+Gt){
			TissueMutation();
			TissueMitosis();
			divcount += 1;
		}
	}
	else{
		assert(pregline.size()!=0);
		while (divcount < Gs+Gt){
			TissueMutation();
			TissueMitosis();
			divcount += 1;
		}		
	}
	// apply mu_B to all
	for (int k=0;k<(L-Gt-Gs);k++){
		assert(L>Gt+Gs);
		for (int i=0;i<(int)Tissues.size();i++){
			for (int j=0;j<(int)Tissues[i].size();j++){
					cell cx = Tissues[i][j];
					int m = cx[0];
					int M = cx[1];
					int dm = gsl_ran_binomial(rngo, mug, M-m);
					Tissues[i][j][0] += dm;
			}
		}
	}
	assert(Tissues.size()==pow(2.0,Gt));
	vector<double> fits;
	fits = vector<double> (0);
	for (int i=0;i<int(Tissues.size());i++){
		fits.push_back(TissueFitness(Tissues[i]));
	}
	if(Gs>0){AdultFitness = *min_element(fits.begin(), fits.end());}
	else {AdultFitness = accumulate(fits.begin(), fits.end(), 0.0) / fits.size();}
	//AdultFitness = *min_element(fits.begin(), fits.end());
}

void Organism::NoDivMitosis(){
	// divide mitotically every cell in org
	// without prior duplication of mitochondria
	// use only in anysogamous species
	vector< vector<int> > newsoma;
	int s1 = soma.size();
	for (int i=0;i<int(soma.size());i++){
		vector<int> cell = soma[i];
		int m = cell[0];
		int M = cell[1];
		int m1, m2;
		m1 = gsl_ran_hypergeometric (rngo, m, M-m, M/2);
		m2 = m - m1;
		M = M/2;
		int cstate1[2] = {m1,M};
		int cstate2[2] = {m2,M};
		vector<int> newcell1 (cstate1, cstate1+2);
		vector<int> newcell2 (cstate2, cstate2+2);
		newsoma.push_back (newcell1);
		newsoma.push_back (newcell2);		
	}
	soma = newsoma;
	int s2 = soma.size();
	assert(s2==2*s1);
}

void Organism::Mitosis(){
	vector< vector<int> > newsoma;
	int s1 = soma.size();
	for (int i=0;i<int(soma.size());i++){
		vector <int> cell = soma[i];
		int m = cell[0];
		int M = cell[1];
		int m1, m2;
		m1 = gsl_ran_hypergeometric (rngo, 2*m, 2*M-2*m, M);
		m2 = 2*m - m1;
		int cstate1[2] = {m1,M};
		int cstate2[2] = {m2,M};
		vector<int> newcell1 (cstate1, cstate1+2);
		vector<int> newcell2 (cstate2, cstate2+2);
		newsoma.push_back (newcell1);
		newsoma.push_back (newcell2);		
	}
	soma = newsoma;
	int s2 = soma.size();
	assert(s2==2*s1);
}

void Organism::TissueMitosis(){
	vector<vector<cell> > newTissues;
	newTissues = vector<vector<cell> > (0);
	for (int i=0;i<int(Tissues.size());i++){
		vector<cell> tissue = Tissues[i];
		int s1 = tissue.size();
		vector<cell> newtissue;
		newtissue = vector<cell> (0);
		for (int j=0;j<tissue.size();j++){
			vector <int> cell = tissue[j];
			int m = cell[0];
			int M = cell[1];
			int m1, m2;
			m1 = gsl_ran_hypergeometric (rngo, 2*m, 2*M-2*m, M);
			m2 = 2*m - m1;
			int cstate1[2] = {m1,M};
			int cstate2[2] = {m2,M};
			vector<int> newcell1 (cstate1, cstate1+2);
			vector<int> newcell2 (cstate2, cstate2+2);
			newtissue.push_back (newcell1);
			newtissue.push_back (newcell2);
		}
		int s2 = newtissue.size();
		assert(s2 = 2*s1);
		newTissues.push_back(newtissue);
	}
	Tissues = newTissues;
}

void Organism::Mutation(double mux){
	for (int i=0;i<(int)soma.size();i++){
		vector<int> cell= soma[i];
		int m = cell[0];
		int M = cell[1];
		int dm = gsl_ran_binomial(rngo, mux, M-m);
		soma[i][0] += dm;			
	}
}

void Organism::TissueMutation(){
	for (int i=0;i<(int)Tissues.size();i++){
		for (int j=0;j<(int)Tissues[i].size();j++){
			vector<int> cell= Tissues[i][j];
			int m = cell[0];
			int M = cell[1];
			int dm = gsl_ran_binomial(rngo, mus+mug, M-m);
			Tissues[i][j][0] += dm;
		}			
	}
}

double Organism::TissueFitness(vector<cell> tissuex){
	vector<double> fits;
	fits = vector<double> (0);
	for (int i=0;i<(int)tissuex.size();i++){
		fits.push_back(1.0 - pow(1.0*tissuex[i][0]/tissuex[i][1],2));//change
	}
	return accumulate(fits.begin(), fits.end(), 0.0) / fits.size();
}

void Organism::CompleteGametogenesis(){
        //meiosis
	int Qp;
        if (Qg[0]==Qg[1]){assert(Qg[0]==0);}
        else {assert(Qg[0]>Qg[1]);}
        Qp = Qg[0];
	assert((int)pregline.size()==1);
	int m = pregline[0][0];
	m = Amplify(m,Qp);
	int M = pregline[0][1];
	assert(M==Ms);
	M = M*pow(2.0,Qp); 
        int m1 = gsl_ran_hypergeometric (rngo, 2*m, 2*M-2*m, M);
        int m2 = 2*m - m1;
        int m11, m21;
        m11 = gsl_ran_hypergeometric (rngo, m1, M-m1, M/2);
        m21 = gsl_ran_hypergeometric (rngo, m2, M-m2, M/2);
        Gamete gam1;
        Gamete gam2;
	Ggg = Mix(Ggg[0],Ggg[1]);
        int nA = 1; // no UPI invasion, so can be anything
        gam1 = Gamete(m11,M/2,nA,WZg[0],Ggg[0],Gt,Gs,Qg[0]);

        assert((int)pregline.size()==1);
        if (Qg[0]==Qg[1]){assert(Qg[0]==0);}
        else {assert(Qg[0]>Qg[1]);}
        Qp = Qg[0];
        m = pregline[0][0];
	m = Amplify(m,Qp);
        M = pregline[0][1];
	M = M*pow(2.0,Qp);
        m1 = gsl_ran_hypergeometric (rngo, 2*m, 2*M-2*m, M);
        m2 = 2*m - m1;
        m11 = gsl_ran_hypergeometric (rngo, m1, M-m1, M/2);
        m21 = gsl_ran_hypergeometric (rngo, m2, M-m2, M/2);
        nA = 1; // no UPI invasion, so can be anything

        gam2 = Gamete(m21,M/2,nA,WZg[1],Ggg[1],Gt,Gs,Qg[1]);
        gline.push_back(gam1);
        gline.push_back(gam2);

}

void Organism::Gametogenesis(){
	assert(pregline.size()==0);
	int ncells = soma.size();
	int r = gsl_rng_uniform_int (rngo, ncells);
	int m = soma[r][0];
	int M = soma[r][1];	
	for (int i=0;i<(L-GgP);i++){
		int dm = gsl_ran_binomial(rngo, mug, M-m);
		m += dm;
	}
	assert(m<=M);
	vector<int> pregam (2,0);
	pregam[0] = m;
	pregam[1] = M;
	pregline.push_back(pregam);
}

void Organism::SomaticGametogenesis(){
	assert(pregline.size()==0);
	int ntissues = Tissues.size();
	int rt = gsl_rng_uniform_int (rngo, ntissues);
	vector<cell> tissue = Tissues[rt];
	int ncells = tissue.size();
	int r = gsl_rng_uniform_int (rngo, ncells);
	int m = tissue[r][0];
	int M = tissue[r][1];
	for (int i=0;i<(L-GgP);i++){
		int dm = gsl_ran_binomial(rngo, mug, M-m);
		m += dm;
	}
	vector<int> pregam (2,0);
	pregam[0] = m;
	pregam[1] = M;
	pregline.push_back(pregam);
}

int Organism::Amplify(int mx, int Qx){
	int m = mx;
	int Ma = Ms;
	for (int i=0;i<Qx;i++){
		int dm = gsl_ran_binomial(rngo, mus, Ma-m);
		m += dm;
		//int mnew = gsl_ran_binomial(rngo, 1.0*m/Ma, 2*Ma); // this generates extra variance
		int mnew = 2*m; // not to generate any extra variance
		m = mnew;
		Ma = 2*Ma;
	}
	assert(Ma==Ms*pow(2.0,Qx));
	return m;
}

vector<int> Organism::Mix(int a, int b){
	int prevec[2] = {a,b};
	vector<int> vec (prevec,prevec+2);
	random_shuffle(vec.begin(), vec.end());
	return vec;
}

vector<int> Organism::Mix4(int a,int b,int c,int d){
	int prevec[4] = {a,b,c,d};
	vector<int> vec (prevec,prevec+4);
	random_shuffle(vec.begin(), vec.end());
	return vec;
}
