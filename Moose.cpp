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
	N = Nx; // population size
	Ms = Msx;
	int Gg = Ggx;
	int Gt = Gtx;
	int Gs = Gsx;
	L = Lx;
	Q = Qx;
	Gginv = Gginvx;
	Ggdom = Ggdomx;
	IM = IMx;
	mus = musx;
	mug = mugx;
	rngm = gsl_rng_alloc (gsl_rng_mt19937); // init random number generator
	gsl_rng_set (rngm, GenSeed());		// seed rng

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
	// run until fixtion or extinction of the invader allele
	for (int i=0; i<1000000; i++){
		FU = AvgFitness();
		GrowAll();
		Selection();
		Fusion();
		gen += 1;
		double Ggfreq = GetGginvFreq();
		if (Ggfreq==0 or Ggfreq==1){
			if (Ggfreq==0){ retval = 0;}
			else if (Ggfreq==1){ retval = 1;}
			break;
		}
	}
}

void Moose::GrowAll(){
	// Run the Organism life cycle one by one
	for (int i=0;i<N;i++){
		Population[i].Grow();
	}
}

void Moose::Selection(){
	// Applies selection to both mating types as a weighted sampling with replacement
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
	vector <double> sampleM = Sample(fitsM, fitsM.size());
	vector <double> sampleF = Sample(fitsF, fitsF.size());
	vector<Organism> newPopulation;
	for (int i=0;i<N/2;i++){
		newPopulation.push_back(Mpop[sampleM[i]]);
		newPopulation.push_back(Fpop[sampleF[i]]);
	}
	Population = newPopulation;
}

void Moose::Fusion(){
	vector<Gamete> Mgams;
	vector<Gamete> Fgams;
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
	// Sample n items from the list of weights with replacement
	// Returns indices of selected items from the list of weights
	double total = accumulate(list.begin(),list.end(), 0.0);
	vector<double> selection;
	int j = 0;
	double w = list[0];
	while (n>0){
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
	// Fuse two gametes, return a zygote as an instance of Organism
	int mn = 0;
	int Mn = 0;
	assert(MGamx.M==Ms/2);
	if (IM==1){
                //BPI
                //int Fmn = gsl_ran_binomial(rngm, 1.0*FGamx.m/FGamx.M, 2*FGamx.M-Ms/2);
                int Fmn = gsl_ran_hypergeometric (rngm, 2*FGamx.m, 2*FGamx.M-2*FGamx.m, 2*FGamx.M-Ms/2);
                int FMn = 2*FGamx.M-Ms/2;
                int Mmn = MGamx.m; int MMn = MGamx.M;
                mn = Fmn + Mmn;
                Mn = FMn + MMn;
	}
	else if (IM==0){
		//UPI fusion
		mn = 2*FGamx.m;
		Mn = 2*FGamx.M;
	}
	int q = log2(1.0*FGamx.M/Ms) + 1;
	assert(q==Q);
	assert(FGamx.M == Ms/2 * pow(2.0,q));
	//mn = gsl_ran_binomial(rngm, 1.0*mn/Mn, Ms*pow(2.0,q));
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
	// Get average fitnes of the whole population
	double fits = 0;
	for (int i=0; i<N; i++){
		fits += Population[i].AdultFitness;
	}
	return fits/N;
}

double Moose::AvgInvFitness(){
	// Get average fitness of the invading allele
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
	// Get the average fitness of the wild type allele
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
	// Get current frequency of the invader allele
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
	// Generates the seed for the random number generator
	// based on current time and process number 
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
	seed += sr;
	return seed;
}
