#include "varstruct.h"
//#define DEBUG_2_2 1

model* initialise_model(int*** data,int nloci,int npops,int ploidy,int ninds){
  model* m=new model();
  //	INIT  no of pops, no of loci, ploidy
  m->npops=npops;
  m->nloci=nloci;
  m->ploidy=ploidy;
  cout<<"Model description is: npops="<<m->npops<<" nloci="<<m->nloci<<" ploidy="<<m->ploidy<<endl;
  //	INIT alpha array
  m->alpha=new double[m->npops];
  for(int i=0;i<m->npops;i++)
    m->alpha[i]=1.0;
  //	INIT centroid arrays
  m->numcentroids=new int[nloci];
  m->centroids=new int*[nloci];
  m->allelefreq=new double**[nloci];
  for(int i=0;i<nloci;i++){
    m->centroids[i]=new int[MAX_CENTROIDS];
    m->allelefreq[i]=new double*[npops];
    for(int j=0;j<npops;j++){
      m->allelefreq[i][j]=new double[MAX_CENTROIDS];
    }
  }
  m->uniqalleles=new map<int,int>[nloci];
  //	INIT centroid parameters
  putinitialvalues(m,data,ninds);
#ifdef DEBUG
  cout<<"Returning brand new model\n";
#endif
  return m;
}

double logit(double x){
  return log(x)-log(1-x);
}

double sigmoid(double x){
  return 1/(1+exp(-1*x));
}

int signum(int x){
  if(x==0)
    return 0;
  else return x/abs(x);
}

int signum(double x){
  if(abs(x)<1e-20)
    return 0;
  else if(x<0)
    return -1;
  else return 1;
}

void putinitialvalues(model* m,int*** data,int ninds){
  int itercount;
  int* thislocus=new int[ninds*m->ploidy];
  int max,min,ncentroids;
  int bestncentroids;
  int allelecount=0;
  for(int i=0;i<m->nloci;i++){
    m->uniqalleles[i] =map<int,int>(); // inited allele unique map for this locus
#ifdef DEBUG_0
    cout<<"Doing locus "<<i<<endl;
#endif
    //		Put locus data into an array.
    allelecount=0;
    for(int j=0;j<ninds*m->ploidy;j++){
      thislocus[j]=data[j%ninds][i][j/ninds];
      if(thislocus[j]>=0){
	if(m->uniqalleles[i].find(thislocus[j])==m->uniqalleles[i].end()){
		m->uniqalleles[i][thislocus[j]]=allelecount; //added new unique allele to list for locus i
		m->centroids[i][allelecount]=thislocus[j];
		for(int k=0;k<m->npops;k++)
			m->allelefreq[i][k][allelecount]=log(zeroonerand());
		allelecount++;
	}
	data[j%ninds][i][j/ninds]=m->uniqalleles[i][thislocus[j]]; // Replace allele with an index into the map
      }
    }
    bestncentroids=allelecount;
    m->numcentroids[i]=bestncentroids;
    for(int k=0;k<m->npops;k++){
      double logsum=-100;
      for(int l=0;l<bestncentroids;l++){
	logsum=sumlogs(logsum,m->allelefreq[i][k][l]);
      }
      for(int l=0;l<bestncentroids;l++){
	m->allelefreq[i][k][l]-=logsum;
      }
    }
  }
  delete[] thislocus;
}

double sumlogs(double a, double b){
  //	Given a=(log A), and b=(log B) -> return log(A+B)
  double ans;
  if(a>=b)
    ans=a+log1p(exp(b-a));
  else
    ans=b+log1p(exp(a-b));
#ifdef DEBUG
  //	cout<<"In sumlog "<<a<<" "<<b<<" "<<ans<<endl;
#endif
  return ans;
}

double psi(double y){
  double x=y;
  //	Digamma function- from http://www.cog.brown.edu/~mj/code/digamma.c
  double result = 0, xx, xx2, xx4;
  //	assert(x > 0); assume x to be >= zero
  for ( ; x < 7; ++x)
    result -= 1/x;
  x -= 1.0/2.0;
  xx = 1.0/x;
  xx2 = xx*xx;
  xx4 = xx2*xx2;
  result += log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
  if(isnan(result)){
    cout<<"Problem computing psi function for y="<<y<<"\n";
  }
  return result;
}

double trigamma(double y){
  double x=y;
  x=x+6;
  double ans=1/(x*x);
  ans=(((((0.075757575757576*ans-0.033333333333333)*ans+0.0238095238095238)*ans-0.033333333333333)*ans+0.166666666666667)*ans+1)/x+0.5*ans;
  for(int i=0;i<=5;i++){
    x=x-1;
    ans=1/(x*x)+ans;
  }
  return ans;
}

void printsatellite(satellite* info, string s){
  string s1=s+".general";
  ofstream ofile;
  ofile.open(s1.c_str());	
  ofile<<info->npops<<" "<<info->nloci<<" "<<info->ploidy<<endl;
  for (int i=0; i<info->nloci; i++) {
    ofile<<info->numcentroids[i]<<endl;
  }
  ofile.close();
  string s2=s+".gamma";
  ofile.open(s2.c_str());	
  for(int i=0;i<info->ninds;i++){
    for(int k=0;k<info->npops;k++){
      ofile<<info->loggamma[i][k]<<" ";
    }
    ofile<<endl;
  }
  ofile.close();
  string s3=s+".theta";
  ofile.open(s3.c_str());	
  for(int i=0;i<info->ninds;i++){
    for(int k=0;k<info->npops;k++){
      ofile<<info->theta[i][k]<<" ";
    }
    ofile<<endl;
  }
  ofile.close();
  string s4=s+".rhos";
  ofile.open(s4.c_str());
  for(int i=0;i<info->ninds;i++){
    for(int p=0;p<info->ploidy;p++){
      for(int j=0;j<info->nloci;j++){
	for(int k=0;k<info->npops;k++){
	  ofile<<info->logrho[i][j][k][p]<<" ";
	}
	ofile<<endl;
      }
    }
  }
  ofile.close();
}

void dumpsatellite(satellite* info,string sin){
  /*char buffer[100];
    sprintf(buffer,"satellite_dump.%s",sin);
    string s(buffer);*/
  string s=sin;
  s=OUTDIR+"/"+s;
  printsatellite(info,s);	
}

double zeroonerand(){
  double num=rand()%INT_MAX;
  double bound=INT_MAX;
  double ans=num/bound;
  //	cout<<"Number is ="<<ans<<endl;
  return ans;
}

satellite* init_from_model(model* m,int numinds,int*** data){
  satellite* info=new satellite();
  info->nloci=m->nloci;
  info->ninds=numinds;
  info->ploidy=m->ploidy;
  info->npops=m->npops;
  info->numcentroids=new int[info->nloci];
  //Initialize number of centroids
  for(int i=0;i<m->nloci;i++){
    info->numcentroids[i]=m->numcentroids[i];
  }
  //Initialize log-gamma, log-rho 
  info->theta=new double*[numinds];
  info->loggamma=new double*[numinds];
  info->logrho=new double***[numinds];
  for(int i=0;i<numinds;i++){
    info->theta[i]=new double[info->npops];
    info->loggamma[i]=new double[info->npops];
    info->logrho[i]=new double**[info->nloci];
    double thetasum=0;
    for(int j=0;j<info->npops;j++){
      //			info->loggamma[i][j]=m->alpha[j]+info->nloci/info->npops;
      double ctemp=zeroonerand();
      info->loggamma[i][j]=ctemp;
      info->theta[i][j]=ctemp;
      thetasum+=ctemp;
    }
    for(int j=0;j<info->npops;j++){
      info->theta[i][j]/=thetasum;
    }
    for(int j=0;j<info->nloci;j++){
      info->logrho[i][j]=new double*[m->npops];
      for(int k=0;k<info->npops;k++){
	info->logrho[i][j][k]=new double[info->ploidy];
	for(int p=0;p<info->ploidy;p++){
	  info->logrho[i][j][k][p]=log(zeroonerand());
	}
      }			
      for(int p=0;p<info->ploidy;p++){
	double logsum=-100;
	for(int k=0;k<info->npops;k++){
	  logsum=sumlogs(logsum,info->logrho[i][j][k][p]);
	}
	for(int k=0;k<info->npops;k++){
	  info->logrho[i][j][k][p]-=logsum;
	}
      }
    }
  }
  return info;
}

double log1p(double x){
  if(abs(x)<1e-15){
    //		cout<<"Had to use approximation!\n";
    return x;
  }
  else
    return log(1+x);
}

double computelogf(int x,int mu,double delta){
  if(x<-1)
    return -1000;
  double ans=(x==mu ? 0 :-1000);
  return ans;
}

double maxdiff(double* x,double* y,int n){
  // returns the max fractional diff between arrays x and y of size n each
  double best=-1;
  for(int i=0;i<n;i++){
    if(best<2*abs(x[i]-y[i])/(x[i]+y[i]))
      best=2*abs(x[i]-y[i])/(x[i]+y[i]);

  }
  return best; 
}

void clearsatellite(satellite* info){
  for(int i=0;i<info->ninds;i++){
    for(int j=0;j<info->nloci;j++){
      for(int k=0;k<info->npops;k++){
	delete[] info->logrho[i][j][k];
      }		
      delete[] info->logrho[i][j];
    }
    delete[] info->logrho[i];
    delete[] info->loggamma[i];
    delete[] info->theta[i];
  }
  delete[] info->logrho;
  delete[] info->loggamma;
  delete[] info->theta;
  delete[] info->numcentroids;
  delete info;
}

satellite* copysatellite(satellite* oldinfo){
  satellite* info=new satellite();
  info->nloci=oldinfo->nloci;
  info->ninds=oldinfo->ninds;
  info->ploidy=oldinfo->ploidy;
  info->npops=oldinfo->npops;
  info->numcentroids=new int[info->nloci];
  //copy number of centroids
  for(int i=0;i<oldinfo->nloci;i++){
    info->numcentroids[i]=oldinfo->numcentroids[i];
  }
  //copy log-gamma, log-rho 
  info->theta=new double*[oldinfo->ninds];
  info->loggamma=new double*[oldinfo->ninds];
  info->logrho=new double***[oldinfo->ninds];
  for(int i=0;i<oldinfo->ninds;i++){
    info->theta[i]=new double[info->npops];
    info->loggamma[i]=new double[info->npops];
    info->logrho[i]=new double**[info->nloci];
    //double thetasum=0;
    for(int j=0;j<info->npops;j++){
      info->loggamma[i][j]=oldinfo->loggamma[i][j];
      info->theta[i][j]=oldinfo->theta[i][j];
    }
    for(int j=0;j<info->nloci;j++){
      info->logrho[i][j]=new double*[oldinfo->npops];
      for(int k=0;k<info->npops;k++){
	info->logrho[i][j][k]=new double[info->ploidy];
	for(int p=0;p<info->ploidy;p++){
	  info->logrho[i][j][k][p]=oldinfo->logrho[i][j][k][p];
	}
      }			
      /*    for(int p=0;p<info->ploidy;p++){
	double logsum=-100;
	for(int k=0;k<info->npops;k++){
	  logsum=sumlogs(logsum,info->logrho[i][j][k][p]);
	}
	for(int k=0;k<info->npops;k++){
	  info->logrho[i][j][k][p]-=logsum;
	}
	}*/
    }
  }
  return info;
}

satellite* infer(model* m,int*** data,int numinds,satellite* oldinfo){
  //	Given a model and the data, return the inferred variational parameters
  //	satellite* oldinfo=init_from_model(m,numinds);//Might be unnecessary
  cout<<"Performing inference\n";
  satellite* currinfo=copysatellite(oldinfo); 
  int itercount;
  map<int,int>::iterator it;
  double gammasum,logsum;
  double* oldgamma=new double[m->npops];
  for(int i=0;i<numinds;i++){
    //		In a loop, do stuff
    itercount=0;
    for(int k=0;k<m->npops;k++)
      oldgamma[k]=1;
    while(itercount++<=10 && maxdiff(oldgamma,currinfo->loggamma[i],m->npops)>1e-2){
      //			update rho
      for(int k=0;k<m->npops;k++)
	oldgamma[k]=currinfo->loggamma[i][k];
#ifdef DEBUG_0
      cout<<"Start rho update\n";
#endif
      gammasum=0;
      for(int k=0;k<m->npops;k++){
	gammasum+=currinfo->loggamma[i][k];
      }			
      for(int j=0;j<m->nloci;j++){
	for(int k=0;k<m->npops;k++){
	  for(int p=0;p<m->ploidy;p++){
		if(data[i][j][p]<0)
			continue;
	    currinfo->logrho[i][j][k][p]=psi(currinfo->loggamma[i][k])-psi(gammasum);
	    //						cout<<currinfo->loggamma[i][k]<<" "<<gammasum<<endl;
	    if(isnan(currinfo->logrho[i][j][k][p])){
	      cout<<"problem computing log rho zero up "<<i<<" "<<j<<" "<<k<<" "<<p<<endl;
	    }
	    if(data[i][j][p]>=0)
		    currinfo->logrho[i][j][k][p]+=m->allelefreq[j][k][data[i][j][p]];
	    if(isnan(currinfo->logrho[i][j][k][p])){
	      cout<<"problem computing log rho first up "<<i<<" "<<j<<" "<<k<<" "<<p<<endl;
	    }
	  }
	}
      }
      //			normalize
      for(int j=0;j<m->nloci;j++){
	for(int p=0;p<m->ploidy;p++){
	  logsum=-1000;
	  if(data[i][j][p]<0)
	    continue;
	  for(int k=0;k<m->npops;k++){
	    logsum=sumlogs(logsum,currinfo->logrho[i][j][k][p]);
	  }
	  for(int k=0;k<m->npops;k++){
	    currinfo->logrho[i][j][k][p]-=logsum;
	    if(isnan(currinfo->logrho[i][j][k][p])){
	      cout<<"problem computing log rho "<<i<<" "<<j<<" "<<k<<" "<<p<<endl;
	    }
	  }
	}
      }
#ifdef DEBUG_0
      cout<<"Start gamma update\n";
#endif
      //			update gamma	
      for(int k=0;k<m->npops;k++){
	currinfo->loggamma[i][k]=m->alpha[k];
	for(int j=0;j<m->nloci;j++){
	  for(int p=0;p<m->ploidy;p++){
	    if(data[i][j][p]>=0){
	      currinfo->loggamma[i][k]+=exp(currinfo->logrho[i][j][k][p]);
	    }
	  }
	}
	if(isnan(currinfo->loggamma[i][k]) || currinfo->loggamma[i][k]<0){
	  cout<<"problem computing loggamma "<<i<<" "<<k<<endl;
	}
	//				cout<<"Gamma is for i="<<i<<",k="<<k<<", val="<<currinfo->loggamma[i][k]<<endl;
      }
    }
    //cout<<"Finished in "<<itercount<<" iterations\n";
    gammasum=0;
    for(int k=0;k<m->npops;k++){
      gammasum+=currinfo->loggamma[i][k];
    }			
    for(int k=0;k<m->npops;k++){
      currinfo->theta[i][k]=currinfo->loggamma[i][k]/gammasum;
    }			
#ifdef DEBUG_0
    if(i%50==0){
      cout<<"Computed var. params for "<<i<<" individuals\n";
    }
#endif
  }
  delete[] oldgamma;
  cout<<"Finished inference \n";
  return currinfo;
}

double lngamma(double y){
  double x=y;
  double z=1/(x*x);
  x=x+6;
  z=(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z+0.083333333333333)/x;
  z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
  if(isnan(z)){
    cout<<"Problem computing log gamma, y="<<y<<endl;
  }
  return z;
}

double computelikelihood(model* m,satellite* info,int*** data){
  double llkhd=0;
  double alphasum=0,gammasum,temp;
  /*cout<<"Printing alpha for fun for npops="<<m->npops<<" populations\n";
    for(int k=0;k<m->npops;k++){
    cout<<"*"<<m->alpha[k]<<":";
    }
    cout<<"Printed alphas\n";*/
  for(int k=0;k<m->npops;k++){
    alphasum+=m->alpha[k];
    llkhd-=lngamma(m->alpha[k]);
  }
  llkhd+=lngamma(alphasum);
  llkhd*=info->ninds;
  if(isnan(llkhd)){
    cout<<"Lkhd nan at alpha term computation\n";
    while(1){}
  }
  //	cout<<"Computed alpha terms\n";
  for(int i=0;i<info->ninds;i++){
    gammasum=0;
    for(int k=0;k<m->npops;k++){
      gammasum+=info->loggamma[i][k];
      //			cout<<"Computing ln gamma for i="<<i<<",k="<<k<<", gamma val= "<<info->loggamma[i][k]<<endl;
      llkhd+=lngamma(info->loggamma[i][k]);
    }
    if(isnan(llkhd)){
      cout<<"Lkhd nan at one gamma term computation\n";
      while(1){}
    }
    //		cout<<"Computed one gamma term\n";
    llkhd-=lngamma(gammasum);
    for(int k=0;k<m->npops;k++){
      temp=(psi(info->loggamma[i][k])-psi(gammasum));
      llkhd+=(m->alpha[k]-info->loggamma[i][k])*temp;
      for(int j=0;j<m->nloci;j++){
	for(int p=0;p<m->ploidy;p++){
	if(data[i][j][p]>=0)
	  llkhd+=(temp-info->logrho[i][j][k][p])*exp(info->logrho[i][j][k][p]);
	}
      }
    }		
    if(isnan(llkhd)){
      cout<<"Lkhd nan at rho term computation\n";
      while(1){}
    }
    //		cout<<"Computed a rho term\n";
    for(int j=0;j<m->nloci;j++){
      //			cout<<"Num centroids for this locus is "<<m->numcentroids[j]<<endl;
      for(int k=0;k<m->npops;k++){
	  for(int p=0;p<m->ploidy;p++){
	    if(data[i][j][p]>=0){
	      llkhd+=exp(info->logrho[i][j][k][p])*m->allelefreq[j][k][data[i][j][p]];
	    }
	}
      }
    }
    if(isnan(llkhd)){
      cout<<"Lkhd nan at final beta term computation\n";
      while(1){}
    }
    //		cout<<" Done individual "<<i<<endl;
  }
  //	cout<<"Computed rest alpha terms\n";
  return llkhd;
}

void dumpmodel(model* m,string sin){
  /*char buffer[100];
    sprintf(buffer,"model_dump.%s",sin);
    string s(buffer);*/
  string s=sin;
  s=OUTDIR+"/"+s;
  printmodel(m,s);
}

void printmodel(model* m,string s){
#ifdef DEBUG_0
  cout<<"Ready to print model now\n";
#endif
  string s1=s+".general";
  ofstream ofile;
  ofile.open(s1.c_str());
  ofile<<m->npops<<" "<<m->nloci<<" "<<m->ploidy<<endl;
  for (int i=0; i<m->nloci; i++) {
    ofile<<m->numcentroids[i]<<endl;
  }
  ofile.close();
  string s2=s+".centroids";
  ofile.open(s2.c_str());
  for(int i=0;i<m->nloci; i++) {
    for(int j=0;j<m->numcentroids[i];j++){
      ofile<<m->centroids[i][j]<<" ";
    }
    ofile<<endl;
  }
  ofile.close();
  string s4=s+".alphas";
  ofile.open(s4.c_str());
  for(int j=0;j<m->npops;j++){
    ofile<<m->alpha[j]<<" ";
  }
  ofile.close();
  string s5=s+".betas";
  ofile.open(s5.c_str());
  for(int j=0;j<m->nloci;j++)
    for(int k=0;k<m->npops;k++){
      ofile<<j<<":"<<k<<endl;
      for(int l=0;l<m->numcentroids[j];l++)
	ofile<<exp(m->allelefreq[j][k][l])<<" ";
      ofile<<endl;
    }
  ofile.close();
			
}

void updateparams(model* m,satellite* info,int*** data){
  //	update beta matrix
  double logsum,temp;
  for(int j=0;j<m->nloci;j++)
    for(int k=0;k<m->npops;k++)
      for(int l=0;l<m->numcentroids[j];l++)
	m->allelefreq[j][k][l]=-1000;
  for(int i=0;i<info->ninds;i++){
  	for(int j=0;j<m->nloci;j++){
	  for(int p=0;p<m->ploidy;p++){
		  if(data[i][j][p]>=0){
    			for(int k=0;k<m->npops;k++)
			  m->allelefreq[j][k][data[i][j][p]]=sumlogs(m->allelefreq[j][k][data[i][j][p]],info->logrho[i][j][k][p]);
		  }
	  }
	}
  }
  for(int j=0;j<m->nloci;j++){
    for(int k=0;k<m->npops;k++){
      logsum=-1000;
      for(int l=0;l<m->numcentroids[j];l++){
	logsum=sumlogs(logsum,m->allelefreq[j][k][l]);
      }
      for(int l=0;l<m->numcentroids[j];l++){
	m->allelefreq[j][k][l]-=logsum;
      }
    }
  }
  //	update alpha vector
  double max_relchange=1;
  double stepsize=0.1;
  double alphasum,psisum;
  double* gvec=new double[m->npops];
  double* qvec=new double[m->npops];
  double z,b,bnum,bden,newalpha;
  double g_constterm=0,i_gammasum;
  for(int i=0;i<info->ninds;i++){
    i_gammasum=0;
    for(int k=0;k<m->npops;k++){
      i_gammasum+=info->loggamma[i][k];
      g_constterm+=psi(info->loggamma[i][k]);
    }
    g_constterm-=psi(i_gammasum);
  }
  int iter=0;
  while(max_relchange>1e-4 && iter++<100){
    alphasum=0;
    psisum=0;
    for(int k=0;k<m->npops;k++){
      alphasum+=m->alpha[k];
      psisum+=psi(m->alpha[k]);
      qvec[k]=-1*info->ninds*trigamma(m->alpha[k]);
    }
    z=info->ninds*trigamma(alphasum);
    for(int k=0;k<m->npops;k++){
      gvec[k]=info->ninds*(psi(alphasum)-psisum)+g_constterm;
    }
    bnum=0,bden=1/z;
    for(int k=0;k<m->npops;k++){
      bnum+=(gvec[k]/qvec[k]);
      bden-=(1/qvec[k]);
    }
    b=bnum/bden;
    for(int k=0;k<m->npops;k++){
      newalpha=m->alpha[k]-stepsize*(gvec[k]-b)/qvec[k];
      if(abs(newalpha-m->alpha[k])/m->alpha[k]>max_relchange)
	max_relchange=abs(newalpha-m->alpha[k])/m->alpha[k];
      m->alpha[k]=newalpha;
    }
    alphasum=0;
    for(int k=0;k<m->npops;k++){
      alphasum+=m->alpha[k];
    }
    if(alphasum>m->npops){
      for(int k=0;k<m->npops;k++){
	m->alpha[k]/=alphasum;
      }	
    }
  }
  delete[] qvec;
  delete[] gvec;
  return;
}

datamatrix* readdata(string filename){
  //	 Format- A line having ninds, nloci, ploidy
  //	Then, for each ind, there are <ploidy> lines, each line containing <nloci> integers
  //	 Therefore the file has <ninds>*<ploidy> lines
  ifstream ifile;
  int temp;
  cout<<"Filename is "<<filename<<endl;
  ifile.open(filename.c_str());
  datamatrix* d=new datamatrix();
  d->ninds=NUMINDS;
  d->nloci=NUMLOCI;
  d->ploidy=PLOIDY;
  d->data=new int**[d->ninds];
  for(int i=0;i<d->ninds;i++){
    d->data[i]=new int*[d->nloci];
    for(int j=0;j<d->nloci;j++){
      d->data[i][j]=new int[d->ploidy];
    }
  }
  for(int i=0;i<d->ninds;i++){
    for(int p=0;p<d->ploidy;p++){
      for(int j=0;j<d->nloci;j++){
	ifile>>temp;
	//cout<<"Reading "<<i<<" "<<p<<" "<<j<<":"<<temp<<endl;
	d->data[i][j][p]=temp;
      }
    }
  }
  ifile.close();
  return d;
}

void clearmatrix(datamatrix* d){
  for(int i=0;i<d->ninds;i++){
    for(int p=0;p<d->ploidy;p++){
      delete[] d->data[i][p];
      }
    delete[] d->data[i];
    }
  delete[] d->data;
  delete d;
}

void learnmodel(string filename,int npops){
#ifdef DEBUG
  cout<<"Starting file read\n";
#endif
  datamatrix* d=readdata(filename);
#ifdef DEBUG
  cout<<"Read file\n";
#endif
  model* currmodel=initialise_model(d->data,d->nloci,npops,d->ploidy,d->ninds);
  //	putinitialvalues(currmodel,d->data,d->ninds);
#ifdef DEBUG
  cout<<"Inited model\n";
#endif
  double currlkhd,oldlkhd;
  satellite* oldinfo=init_from_model(currmodel,d->ninds,d->data);
  satellite* info;
  //	currlkhd=computelikelihood(currmodel,info,d->data);
  currlkhd=-1e20;
  oldlkhd=-1e25;
  int iter=0;
  string iterstring="";
  char buffer[10];
  while( (currlkhd-oldlkhd)/abs(oldlkhd)>1e-4 && iter<150){
    iter++;
    cout<<"Start iteration "<<iter<<" with llkhd= "<<currlkhd<<endl;
    oldlkhd=currlkhd;
    info=infer(currmodel,d->data,d->ninds,oldinfo);
    clearsatellite(oldinfo);
    currlkhd=computelikelihood(currmodel,info,d->data);
    cout<<"After inference, lkhd = "<<currlkhd<<endl;
    updateparams(currmodel,info,d->data);
    currlkhd=computelikelihood(currmodel,info,d->data);
    cout<<"After parameter update, lkhd = "<<currlkhd<<endl;
    oldinfo=info;
    //clearsatellite(info);
    if(iter%3==1){
      sprintf(buffer,"%d",iter);
      iterstring.assign(buffer);
      cout<<"Iterstring is "<<iterstring<<endl;
      cout<<"Intermediate print for iteration "<<iter<<endl;
      dumpmodel(currmodel,iterstring);
      dumpsatellite(info,iterstring);
    }
    cout<<"Finished iteration\n";
  }
  dumpmodel(currmodel,"final");
  dumpsatellite(info,"final");
  clearmatrix(d);
}

void readglobals(char* filename){
  FILE* f=fopen(filename,"r");
  float pt1;
  //fscanf(f,"pseudostrength=%f\n", &pt1);
  fclose(f);
  //pseudostrength=pt1;
  //cout<<"Pseudocount strength is "<<pseudostrength<<endl;
}

void valuegiven(int count,int len,char* arr[]){
	if(count+1==len){
		cout<<"No value provided for switch "<<arr[count]<<endl;
		exit(1);
	}
	else if(arr[count+1][0]=='-'){
		cout<<"No value provided for switch "<<arr[count]<<endl;
		exit(1);
	}
}

void verifyoption(int option,string message,string optswitch){
	if(option==-1){
		cout << "Must provide "<<message<<" with switch "<<optswitch<<endl;
		cout << "For help, use the switch \" -h\"\n";
		exit(1);
	}		
}

void usage(){
	cout<<"Usage: varStruct requires the following parameters in any order:\n";
	cout<<"-d <datafile>\n";
	cout<<"-o <output-directory>\n";
	cout<<"-k <number of populations>\n";
	cout<<"-n <number of individuals>\n";
	cout<<"-m <number of loci>\n";
	cout<<"-p <ploidy>\n";
	cout<<"You can provide a random seed (optional) with \"-r <integer seed>\" \n";
}

void readopts(int len,char* arr[]){
	int count=1,len1;
	string temp;
	int dgiven=-1,ngiven=-1,ogiven=-1,pgiven=-1,rgiven=-1,ggiven=-1,kgiven=-1,mgiven=-1;
	if(len==1){
		usage();
		exit(0);
	}
	while(count!=len){
		temp=string(arr[count]);
		len1=temp.size();
		if(len1!=2){
			cout<<"Invalid switch "<<arr[count]<<endl;
			exit(1);
		}
		if(arr[count][1]!='h')
			valuegiven(count,len,arr);
		switch (arr[count][1]) {
			case 'h':
				//help option
				usage();
				exit(0);				
			case 'd':
				// Data file
				DATAFILE=string(arr[count+1]); dgiven=1;
				break;
			case 'o':
				//Output directory
				OUTDIR=string(arr[count+1]); ogiven=1;
				break;
			case 'k':
				// No. of pops
				NPOPS=atoi(arr[count+1]); kgiven=1;
				break;
			case 'n':
				// No. of individuals
				NUMINDS=atoi(arr[count+1]); ngiven=1;
				break;
			case 'm':
				// No. of loci
				NUMLOCI=atoi(arr[count+1]); mgiven=1;
				break;
			case 'r':
				// Random seed
				SEED=atoi(arr[count+1]); rgiven=1;
				break;
			case 'p':
				// Ploidy
				PLOIDY=atoi(arr[count+1]); pgiven=1;
				break;
			default:
				cout<<"Unknown switch "<<arr[count]<<endl;
				exit(1);
		}
		count+=2;
	}	
	verifyoption(dgiven,"Datafile","-d");
	verifyoption(ngiven,"Number of individuals","-n");
	verifyoption(mgiven,"Number of loci","-m");
	verifyoption(pgiven,"Ploidy","-p");
	verifyoption(kgiven,"Number of populations","-k");
	verifyoption(ogiven,"Output directory","-o");
	if(rgiven==-1)
		cout<<"No seed provided, using default value: "<<SEED<<endl;
}

int main(int argc,char** argv){
  readopts(argc,argv);
  srand(SEED);
  learnmodel(DATAFILE,NPOPS);
  return 0;
}
