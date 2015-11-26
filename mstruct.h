#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<set>
#include<ctime>
#include<cstdlib>
#include<map>
#include<climits>
#include<algorithm>

using namespace std;

int PRINTCT=0,SATPRINTCT=0;
#define DEBUG 1
//#define DEBUG_0 1
int MAXALLELES=50;

struct datamatrix{
  int*** data;
  int ninds,nloci,ploidy;
};

struct model{
  int nloci,npops,ploidy;
  int** centroids;
  int* numcentroids;
  double* alpha;
  double*** allelefreq;
  map<int,int>* uniqalleles; //stores unique alleles at each locus to reduce computation
};

struct satellite{
  int ninds,nloci,npops,ploidy;
  int* numcentroids;
  double** loggamma;
  double**** logrho;
  double** theta;
};

const int MAX_CENTROIDS=20;
const int MAX_ITER=100;
const double THRESH=1e-5;
string OUTDIR;

//data is ninds*nloci*ploidy
double log1p(double);

model* initialise_model(int*** data,int nloci,int npops,int ploidy,int ninds);
void putinitialvalues(model* m,int*** data,int ninds);
void infermodel(model* m,int*** data,int ininds,double**** varz,double** vartheta);
void updatemodel(model* m,int*** data,int ninds,double**** varz,double** vartheta);
double computelogf(int x, int mu);
void updateparams(model* m,satellite* info);
void learnmodel(string filename,int npops);

satellite* init_from_model(model* m,int numinds,int*** data);
satellite* makecopy(satellite* current);
double sumlogs(double a, double b);
void printmodel(model* m,string s);
void dumpmodel(model* m,string s);
void printsatellite(satellite* info,string s);
void dumpsatellite(satellite* info,string s);

double zeroonerand();
double logit(double x);
double sigmoid(double x);
double psi(double x);
datamatrix* readdata(string filename);
