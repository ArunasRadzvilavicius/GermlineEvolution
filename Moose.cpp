#include "Moose.h"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <math.h> 
using namespace std;

Moose::Moose(){
}

Moose::Moose(int Nx,int Ggx,int Gtx,int Gsx,int Msx,double musx,double mugx,int Gginvx,int IMx,int Ggdomx,int Lx,int Qx){
	srand( time(0) );
	N = Nx;
	Ms = Msx;
	int Gg = Ggx; // ancestral homozygote
	int Gt = Gtx; // gamma_T
	int Gs = Gsx; // gamma_S
	L = Lx;
	Q = Qx;
	Gginv = Gginvx;
	Ggdom = Ggdomx;
	IM = IMx;
	mus = musx;
	mug = mugx;
	rngm = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rngm, GenSeed());

	// initialize population
	Organism org;
	int m0;
	gtype WZg;
	gtype Ggg;
	Ggg = gtype (2,Gg);
	gtype Qgx; Qgx = gtype (2);
	for (int i=0;i<N;i++){
		WZg = gtype (2,0);
		if (i<N/2) {
			// females WZ
			WZg[0] = 1;
			WZg[1] = 0;
			Qgx[0] = Q;
			Qgx[1] = 0;
		}
		else {
			// males ZZ
			WZg[0] = 0;
			WZg[1] = 0;
			Qgx[0] = 0;
			Qgx[1] = 0;
		}
		
		m0 = gsl_rng_uniform_int(rngm, Ms);
		org = Organism(m0,Ms,WZg,Ggg,Gt,Gs,rngm,mus,mug,Ggdom,L,Qgx,Ms);
		Population.push_back( org );
	}

	int gen = 0;
	double FU = 0;
	double FA = 0;
	double Fwt = 0;
	double Fmut = 0;
	// equilibrate
	for (int i=0; i<50; i++){
		FU = AvgFitness();
		GrowAll();
		FA = AvgFitness();
		//Fwt = AvgWTFitness();
		//Fmut = AvgInvFitness();
		Selection();
		Fusion();
		gen += 1;
		//cout<<gen<<' '<<FU<<' '<<FA<<' '<<GetGginvFreq()<<endl;
		//cout<<gen<<' '<<Fwt<<' '<<Fmut<<' '<<GetGginvFreq()<<endl;
	}

	//insert some Gg mutants
	int nmut = 0;
	int i = 0;
	while (nmut < 0.05*2*N){
		if (Population[i].Ggg[0]==Gg){
			Population[i].Ggg[0]=Gginv;
			nmut+=1;
		}
		i+=1;
	}
	// evolve further, calc avg
	double FUsum = 0;
	double FAsum = 0;
	double ntot = 0;
	for (int i=0; i<1000000; i++){
		FU = AvgFitness();
		//FUsum += FU;
		GrowAll();
		//FA = AvgFitness();
		//Fwt = AvgWTFitness();
		//Fmut = AvgInvFitness();
		//FAsum += FA;
		Selection();
		Fusion();
		gen += 1;
		ntot += 1;
		double Ggfreq = GetGginvFreq();
		if (Ggfreq==0 or Ggfreq==1){
			if (Ggfreq==0){ retval = 0;}
			else if (Ggfreq==1){ retval = 1;}
			break;
		}
		//cout<<gen<<' '<<FU<<' '<<FA<<' '<<Ggfreq<<endl;
		//cout<<gen<<' '<<Fwt<<' '<<Fmut<<' '<<GetGginvFreq()<<endl;
	}
}

void Moose::GrowAll(){
	for (int i=0;i<N;i++){
		Population[i].Grow();
	}
}

void Moose::Selection(){
	// separate MTypes
	vector<Organism> Mpop;
	vector<Organism> Fpop;
	for (int i=0;i<N;i++){
		assert(Population[i].gline.size()==2);
		if (Population[i].WZg[0]!=Population[i].WZg[1]) {Fpop.push_back(Population[i]);}
		else {assert(Population[i].WZg[0]==0);Mpop.push_back(Population[i]);}
	}
	assert(int(Fpop.size())==N/2);
	assert(int(Mpop.size())==N/2);
	vector<double> fitsM (N/2);
	vector<double> fitsF (N/2);
	for (int i=0;i<N/2;i++){
		fitsM[i] = Mpop[i].AdultFitness;
		fitsF[i] = Fpop[i].AdultFitness;
	}
	// now sample with replacement according to weights
	vector <double> sampleM = Sample(fitsM, fitsM.size());
	vector <double> sampleF = Sample(fitsF, fitsF.size());
	// PrintVec(sample);
	vector<Organism> newPopulation;
	for (int i=0;i<N/2;i++){
		newPopulation.push_back(Mpop[sampleM[i]]);
		newPopulation.push_back(Fpop[sampleF[i]]);
	}
	Population = newPopulation;
}

void Moose::Fusion(){
	vector<Gamete> Mgams;
	vector<Gamete > Fgams;
	for (int i=0;i<N;i++){
		assert( int(Population[i].gline.size())==2 );
		if (Population[i].WZg[0]!=Population[i].WZg[1]){
			Fgams.push_back(Population[i].gline[0]);
			Fgams.push_back(Population[i].gline[1]);
		}
		else{
			assert(Population[i].WZg[0]==0);
			Mgams.push_back(Population[i].gline[0]);
			Mgams.push_back(Population[i].gline[1]);
		}
	}
	assert(int(Fgams.size())==N);
	assert(int(Mgams.size())==N);
	random_shuffle(Fgams.begin(), Fgams.end());
	random_shuffle(Mgams.begin(), Mgams.end());
	vector<Organism> newPopulation(N);
	for (int j=0;j<N;j++){
		Gamete GamM = Mgams[j];
		Gamete GamF = Fgams[j];
		vector<Organism> newOrgs = Fuse(GamM,GamF);
		assert(int(newOrgs.size())==1);
		newPopulation[j]=newOrgs[0];
	}
	Population = newPopulation;
}

vector<double> Moose::Sample(vector<double> list, int n){
	// sample n items from the list of weights with replacement
	// returns indices of selected items from the list of weights
	double total = accumulate(list.begin(),list.end(), 0.0);
	vector<double> selection;
	int j = 0;
	double w = list[0];
	while (n>0){
		// draw a single weighted sample
		double r = gsl_rng_uniform (rngm);
		double x = total*(1.0 - pow(r, 1.0/n));
		total = total - x;
		while (x > w){
			x = x-w;
			j+=1;
			w = list[j];
		}
		w-=x;
		n-=1;
		selection.push_back(j);
	}
	return selection;
}

vector<Organism> Moose::Fuse(Gamete MGamx, Gamete FGamx){
	int mn = 0;
	int Mn = 0;
	assert(MGamx.M==Ms/2);
	//if (FGamx.nA==1){
	if (IM==1){
		//BPI
		mn = MGamx.m + FGamx.m;
		Mn = MGamx.M + FGamx.M;
	}
	//else if (FGamx.nA==0){
	else if (IM==0){
		//UPI
		mn = FGamx.m;
		Mn = FGamx.M;
	}
	//resample to M*2^Q
	int q = log2(1.0*FGamx.M/Ms) + 1;
	assert(q==Q);
	//assert(q==0);
	assert(FGamx.M == Ms/2 * pow(2.0,q));
	mn = gsl_ran_binomial(rngm, 1.0*mn/Mn, Ms*pow(2.0,q));
	Mn = Ms*pow(2.0,q);
	Organism Org;
	gtype Gggx; Gggx = gtype (2,0);
	Gggx[0] = FGamx.Gg; Gggx[1] = MGamx.Gg;
	gtype WZgf; WZgf = gtype (2,0);
	WZgf[0] = FGamx.WZ; WZgf[1] = MGamx.WZ;
	gtype Qgx; Qgx = gtype (2,0);
	Qgx[0] = FGamx.Q; Qgx[1] = MGamx.Q;
	Org = Organism(mn,Mn,WZgf,Gggx,MGamx.Gt,MGamx.Gs,rngm,mus,mug,Ggdom,L,Qgx,Ms);
	vector<Organism> Orgs;
	Orgs.push_back(Org);
	return Orgs;
}

double Moose::AvgFitness(){
	double fits = 0;
	for (int i=0; i<N; i++){
		fits += Population[i].AdultFitness;
	}
	return fits/N;
}

double Moose::AvgInvFitness(){
	double fits = 0;
	int ntot = 0;
	for (int i=0; i<N; i++){
		if (Population[i].GgP==Gginv){
			fits += Population[i].AdultFitness;
			ntot += 1;
		}
	}
	return fits/ntot;
}

double Moose::AvgWTFitness(){
	double fits = 0;
	int ntot = 0;
	for (int i=0; i<N; i++){
		if (Population[i].GgP!=Gginv){
			fits += Population[i].AdultFitness;
			ntot += 1;
		}
	}
	return fits/ntot;
}


double Moose::GetGginvFreq(){
	double c = 0;
	for (int i=0;i<N;i++){
		assert(int( Population[i].soma.size() )==1); 
		if (Population[i].Ggg[0] == Gginv) {c+=1;}
		if (Population[i].Ggg[1] == Gginv) {c+=1;}
	}
	return c/(2.0*N);
}

void Moose::PrintPop(){
	for (int i=0;i<N;i++){
		cout<<Population[i].AdultFitness<<endl;
	}
	cout << endl;
}

long Moose::GenSeed(){
	long s, seed, pid;
	int sr;
	time_t seconds;
	pid = getpid();
	s = time(&seconds);
	seed = abs(((s*181)*((pid-83)*359))%104729);
	ifstream file("/dev/urandom", std::ios::binary);
	if (file.is_open()){
	char * memblock;
	int size = sizeof(int);
	memblock = new char [size];
	file.read(memblock, size);
	file.close();
	sr = *reinterpret_cast<int*>(memblock);
	delete[] memblock;
	}
	else{sr = 1 ;}
	//cout << seed << endl;
	seed += sr;
	//cout << seed << endl;
	return seed;
}
