//Swift GRB population study code.
//Written by John G Baker NASA-GSFC (2016)

//#include <valarray>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "omp.h"
//#include "options.hh"
#include "bayesian.hh"
//#include "proposal_distribution.hh"

//This is a very hacky way to access the contents of Phil's codes, which require a runargs struct that we don't really use
extern "C" {
  #include "utils.h"
}

#define ZMIN		0.0
#define ZMAX		10.0

struct {
	double tobs=0.8;
	MLmethod method=RANDOMFOREST;
	bool vary_z1;
}runargs;

extern "C" {
  #include "mock_sample_functions.h"     
  #include "poputils.c"
}

using namespace std;

#define Z1DATA			3.60

///Set up a bayes_data component
///
///In particular, this object just reads in a list of z-values for observed values from a text file.
class GRBpop_z_only_data : public bayes_data {
  int Nevents;
  int ndetdata=0;
  int ndetpop;
public:
  GRBpop_z_only_data():bayes_data(),Nevents(0){};
  int size()const{return Nevents;};
  virtual double getFocusLabel(bool original=false)const{return 0;};
  void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    haveWorkingStateSpace();
  };
  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
    addOption("file","Filename of z-only data.");
  };
  virtual void setup(){
    haveSetup();
    ///Set up the output stateSpace for this object
    stateSpace space(0);
    nativeSpace=space;
    setPrior(new uniform_dist_product(&nativeSpace,0));//we don't need a prior, since space is zero-D, but we need to have set something
    string filename;
    *optValue("file")>>filename;
    ifstream file(filename.c_str());
    if(file.good()){
      string line;
      while(getline(file,line)){
	if(line[0]=='#')continue;//skip comment lines
	double z;
	stringstream(line)>>z;
	labels.push_back(z);
      }
    } else {
      if(filename.size()>0){//empty path signifies go forward without data
	cerr << "Error: " << strerror(errno)<<endl;
	cout<<"GRBpop_z_only_data::setup: Could not open file '"<<filename<<"'."<<endl;
	exit(1);
      }
    }
    Nevents=values.size();
    cout<<"Data read in from file with "<<ndetdata<<" detected GRBs.\n"<<endl;
    //logpois0 = -0.5 * log(2.0 * M_PI * (double) ndetdata);
    haveData(LABELS);
  };
};

///Define one-break GRB population redshift dependence model
///
class GRBpop_one_break_z_signal : public bayes_signal{
  int idx_n0, idx_n1, idx_n2, idx_z1;
  bool varying_z1;
  bool using_flat_n0=false;
  stateSpace localSpace;
  shared_ptr<const sampleable_probability_function> localPrior;
public:
  GRBpop_one_break_z_signal(){
    idx_n0=idx_n0=idx_n2=idx_z1=-1;
    localPrior=nullptr;
    varying_z1=false;
    //This is part of the set-up for Phil's codes.
    load_splines();
  };
  ~GRBpop_one_break_z_signal(){};
  ///From bayes_signal
  virtual std::vector<double> get_model_signal(const state &st, const std::vector<double> &zvals)const{
    double n0,n1,n2,z1;
    n0=st.get_param(idx_n0);
    n1=st.get_param(idx_n1);
    n2=st.get_param(idx_n2);
    if(varying_z1)z1=st.get_param(idx_n2);
    else z1=Z1DATA;
    double lnew = -1.0 * GRBRateIntegral(n0, n1, n2, z1);
    int i;
    std::vector<double> rates;
    for (i = 0; i < zvals.size(); i++)rates.push_back(GRBRate(zvals[i], n0, n1, n2, z1));
    rates.push_back(GRBRateIntegral(n0, n1, n2, z1));//For Poisson data, last entry is the expected total number
    return rates;
  };
  ///From StateSpaceInterface (via bayes_signal)
  void defWorkingStateSpace(const stateSpace &sp){
    checkSetup();//Call this assert whenever we need options to have been processed.
    idx_n0=sp.requireIndex("n0");
    idx_n1=sp.requireIndex("n1");
    idx_n2=sp.requireIndex("n2");
    if(varying_z1)idx_z1=sp.requireIndex("z1");
    haveWorkingStateSpace();
  };
  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
    opt.add(Option("varyz1","Also vary z1."));
    opt.add(Option("nstar","Not yet implemented."));
    opt.add(Option("ntotal","Not yet implemented."));
    opt.add(Option("flatn0","Set uniform prior for n0, otherwise uniform in log(n0)."));
  };
  void setup(){
    if(optSet("flatn0"))using_flat_n0=true;
    cout<<"     flatn0 = "<<(using_flat_n0?"true":"false")<<endl;
    if(optSet("varyz1")){
      varying_z1=true;
      cout<<"Varying z1."<<endl;
    }
    if(optSet("nstar"))cout<<"Ignoring --nstar flag.  Not yet implemented."<<endl;
    if(optSet("ntotal"))cout<<"Ignoring --ntotal flag.  Not yet implemented."<<endl;
    haveSetup();
    //set nativeSpace
    int npar=3;
    if(varying_z1)npar=4;
    cout<<"signal: npar="<<npar<<endl;
    stateSpace space(npar);
    string names[]={"n0","n1","n2","z1"};
    space.set_names(names);  
    nativeSpace=space;
    const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, pol=mixed_dist_product::polar, log=mixed_dist_product::log; 
    //Priors from Phil's code
    //n0 = CubeToLogPrior(Cube[0], 0.01, 2.00);  
    //n1 = CubeToFlatPrior(Cube[1], 0.00, 4.00);
    //n2 = CubeToFlatPrior(Cube[2], -6.00, 0.00);
    //if(--varyz1):z1 = CubeToFlatPrior(Cube[3], 0.00, 10.0);
    //double CubeToFlatPrior(double r, double xmin, double xmax){return r*(xmax-xmin)+xmin;}
    //double CubeToLogPrior(double r, double xmin, double xmax){double lmin=log(xmin),lmax=log(xmax);return exp(r*(lmax-lmin)+lmin);}
    //nstar = n0 * pow(1.0 + z1, n1);
    //ntotal = GRBNumberIntegral(n0, n1, n2, z1);
    valarray<double>    centers((initializer_list<double>){0.01*sqrt(200.0),   2.0,  -3.0,  5.0});
    valarray<double> halfwidths((initializer_list<double>){     sqrt(200.0),   2.0,   3.0,  5.0});
    valarray<int>         types((initializer_list<int>)   {             log,   uni,   uni,  uni});
    if(using_flat_n0){
      centers[0]=1.005;
      halfwidths[0]=0.995;
      types[0]=uni;
    }
    //The rest of these options (from Phil's code) are not yet implemented...
    //if(--nstar):
    //  nstar = CubeToLogPrior(Cube[0], 0.10, 10000.0);
    //  n0 = nstar * pow(1.0 + z1, -n1);
    //  ntotal = GRBNumberIntegral(n0, n1, n2, z1);
    //else if(--ntotal):
    //  ntotal = CubeToLogPrior(Cube[0], 1.00, 1e5);
    //  double ntmp = GRBNumberIntegral(1.0, n1, n2, z1);
    //  n0 = ntotal / ntmp;
    //  nstar = n0 * pow(1.0 + z1, n1);
    //else: [default explicit n0 prior]
    //  if(--flatn0):
    //    n0 = CubeToFlatPrior(Cube[0], 0.01, 2.00);
    setPrior(new mixed_dist_product(&nativeSpace,types,centers,halfwidths));
  };
};

///Bayes class for a likelihood function object
///
///The interface here is probably not everything we want, but is enough for what was already in the main function.
class GRBpop_likelihood : public bayes_likelihood {
  bool zeroLogLike=false;
public:
  GRBpop_likelihood(bayes_data *data,bayes_signal *signal):bayes_likelihood(nullptr,data,signal){};
  void addOptions(Options &opt,const string &prefix=""){
    Optioned::addOptions(opt,prefix);
    opt.add(Option("zeroLogLike","Set the log-likelihood to zero (thus purely sampling from prior)."));
  };
  void setup(){
    if(optSet("zeroLogLike"))zeroLogLike=true;
    cout<<"zeroLogLike = "<<(zeroLogLike?"true":"false")<<endl;
    bayes_likelihood::setup();
  };
  double evaluate_log(state &s){
    valarray<double>params=s.get_params();
    double result=0;
    if(not zeroLogLike)result=log_poisson(s);
    double post=result;
    post+=nativePrior->evaluate_log(s);//May need a mechanism to check that Prior is set
    #pragma omp critical
    {     
      if(post>best_post){
        best_post=post;
        best=state(s);
      }
      if(!isfinite(result)){
        cout<<"Whoa dude, loglike is NAN! What's up with that?"<<endl;
        cout<<"params="<<s.get_string()<<endl;
	result=-INFINITY;
      }
    }
    return result;
  };
};
  



  

  
