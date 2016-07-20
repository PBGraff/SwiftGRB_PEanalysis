//Swift GRB code...
//Written by John G Baker NASA-GSFC (2016)

//#include "mlfit.hh"
#include <valarray>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include "omp.h"
#include "options.hh"
//#include <mcheck.h>
#include "bayesian.hh"
#include "proposal_distribution.hh"
#include "ptmcmc.hh"
#include "GRB_bayesian_components.hh"

using namespace std;

shared_ptr<Random> globalRNG;//used for some debugging... 

//***************************************************************************************8
//main test program
int main(int argc, char*argv[]){

  Options opt(true);
  //Create the sampler
  ptmcmc_sampler mcmc;
  bayes_sampler *s0=&mcmc;
  //Create the model components and likelihood;
  bayes_data *data=new GRBpop_z_only_data();
  bayes_signal *signal=new GRBpop_one_break_z_signal();
  bayes_likelihood *like=new GRBpop_likelihood(data,signal);
  
  //prep command-line options
  s0->addOptions(opt);
  data->addOptions(opt);
  signal->addOptions(opt);
  like->addOptions(opt);

  //Add some command more line options
  opt.add(Option("nchains","Number of consequtive chain runs. Default 1","1"));
  opt.add(Option("seed","Pseudo random number grenerator seed in [0,1). (Default=-1, use clock to seed.)","-1"));
  opt.add(Option("precision","Set output precision digits. (Default 13).","13"));
  opt.add(Option("outname","Base name for output files (Default 'mcmc_output').","mcmc_output"));
  
  int Nlead_args=1;

  bool parseBAD=opt.parse(argc,argv);
  if(parseBAD) {
    cout << "Usage:\n mcmc [-options=vals] " << endl;
    cout <<opt.print_usage()<<endl;
    return 1;
  }
    
  cout<<"flags=\n"<<opt.report()<<endl;

  data->setup();  
  signal->setup();  
  like->setup();

  double seed;
  int Nchain,output_precision;
  int Nsigma=1;
  int Nbest=10;
  string outname;
  ostringstream ss("");
  istringstream(opt.value("nchains"))>>Nchain;
  istringstream(opt.value("seed"))>>seed;
  //if seed<0 set seed from clock
  if(seed<0)seed=fmod(time(NULL)/3.0e7,1);
  istringstream(opt.value("precision"))>>output_precision;
  istringstream(opt.value("outname"))>>outname;

  //report
  cout.precision(output_precision);
  cout<<"\noutname = '"<<outname<<"'"<<endl;
  cout<<"seed="<<seed<<endl; 
  cout<<"Running on "<<omp_get_max_threads()<<" thread"<<(omp_get_max_threads()>1?"s":"")<<"."<<endl;

  //Should probably move this to ptmcmc/bayesian
  ProbabilityDist::setSeed(seed);
  globalRNG.reset(ProbabilityDist::getPRNG());//just for safety to keep us from deleting main RNG in debugging.
  
  //Get the space/prior for use here
  stateSpace space;
  shared_ptr<const sampleable_probability_function> prior;  
  space=*like->getObjectStateSpace();
  cout<<"like.nativeSpace=\n"<<space.show()<<endl;
  prior=like->getObjectPrior();
  cout<<"Prior is:\n"<<prior->show()<<endl;
  valarray<double> halfw;prior->getHalfwidths(halfw);

  //Read Params
  int Npar=space.size();
  cout<<"Npar="<<Npar<<endl;
  
  //Bayesian sampling [assuming mcmc]:
  //Set the proposal distribution 
  int Ninit;
  proposal_distribution *prop=ptmcmc_sampler::new_proposal_distribution(Npar,Ninit,opt,prior.get(),&halfw);
  cout<<"Proposal distribution is:\n"<<prop->show()<<endl;
  //set up the mcmc sampler (assuming mcmc)
  mcmc.setup(Ninit,*like,*prior,*prop,output_precision);


  //Prepare for chain output
  ss<<outname;
  string base=ss.str();

  //Loop over Nchains
  for(int ic=0;ic<Nchain;ic++){
    bayes_sampler *s=s0->clone();
    s->initialize();
    s->run(base,ic);
    s->analyze(base,ic,Nsigma,Nbest,*like);
    delete s;
  }
  
  //Dump summary info
  cout<<"best_post "<<like->bestPost()<<", state="<<like->bestState().get_string()<<endl;
  delete data;
  delete signal;
  delete like;
}

