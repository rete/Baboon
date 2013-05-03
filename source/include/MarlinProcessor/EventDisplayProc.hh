
#ifndef EVENTDISPLAYPROC_HH
#define EVENTDISPLAYPROC_HH

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TBranch.h>
#include <EVENT/CalorimeterHit.h>

using namespace lcio ;
using namespace marlin ;

const size_t tetamax=100;
typedef struct{
  int I;
  int J;
  int K;
  int TH;
}MyHit_t;

typedef struct{
  int dens2D;
  int dens3D;
  int pos[3];
}Density_t;

typedef struct{
  std::vector<EVENT::CalorimeterHit*> hits;
  int size;
  int begin;
  int end;
  std::vector<float> pos;
  float radius;
}MyCore_t;

typedef struct{
  std::vector<EVENT::CalorimeterHit*> cluster;
  float pos[3];
}Cluster;

typedef struct{
  float pos[3];
  int index;
  int rhox[tetamax];
  int rhoy[tetamax];
}IsATrack_t;

class Track
{
 public:
  Track(){;}
  ~Track();
  void getTrackStartingPoint();
  void getLastTrackPoint();
  //  void getTrackHit();
  //std::vector<EVENT::CalorimeterHit*> TrackHit;
  float getTrackAngle();
  std::vector<Cluster> TrackClus;
  int begin;
  int end;
};

class EventDisplayProc : public Processor {

 public:

  virtual Processor*  newProcessor() { return new EventDisplayProc ; }

  
  EventDisplayProc() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  std::vector<int> Nhit();
  virtual int Nlayer();
  virtual int Begin();
  virtual int End();
  virtual float CenterOfXGravity();
  virtual float CenterOfYGravity();
  virtual float CenterOfZGravity(int begin, int end);
  virtual float Radius(int Zbegin, float* cog);
  virtual int NhitLastPlan();
  virtual void fillTree();
  virtual void ClearVector();
  std::vector<Density_t> Density();
  virtual int Edge();
  virtual int holeFinder(int begin, float *cog);
  virtual int hole1Finder(int begin, float *cog);
  std::vector<MyHit_t> getIJK();
  std::vector<MyCore_t> CoreFinder(std::vector<Density_t> dens);
  virtual int NumberOfCores(std::vector<MyCore_t> vec);
  virtual int CoreSize(MyCore_t core);
  virtual int CoreBegin(MyCore_t core);
  virtual int CoreEnd(MyCore_t core);
  std::vector<float> getCorePosition(MyCore_t core);
  virtual float getCoreRadius(MyCore_t core);
  std::vector<Cluster> GetClusters();
  std::vector<bool> TagTracks(std::vector<Cluster> clusters);
  std::vector<Track*> getRealTrack(std::vector<bool> Tag, std::vector<Cluster> clusters);
  std::vector<int> HitInTrack(std::vector<bool> Tag, std::vector<Cluster> clusters);
  
 protected:
  
  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hcalCollections;

  std::string _outputHcalCollection;
  std::string _outputRelCollection;

  //float _thresholdHcal;
  std::vector<float> _thresholdHcal;
  std::string treeFileName_;
  int _digitalHcal;

  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _hcalLayers;
 private:
  int numElements;
  LCCollection * col;

  TFile *file;
  TTree *tree;
  TBranch *N;
  TBranch *NLay;
  TBranch *ZBegin;
  TBranch *RAdius;
  TBranch *E;
  TBranch *I;
  TBranch *J;
  TBranch *K;
  TBranch *TH;
  TBranch *DEns;
  TBranch *EH;
  TBranch *COG;
  TBranch *HOle;
  TBranch *HOle1;
  TBranch *NCores;
  TBranch *COreSIze;
  TBranch *COreBegin;
  TBranch *COreEnd;
  TBranch *COrePosK;
  TBranch *COreRadius;
  TBranch *TRackTag;
  TBranch *TRackBegin;
  TBranch *TRackEnd;
  TBranch *TRackLength;
  TBranch *TRackMulti;

  
  std::vector<bool> tagtrack;
  int energy_;
  int energy;
  int nhit;
  int begin;
  int nlayer;
  float radius;
  float EHratio;
  int hole;
  int hole1;
  float cog[3];
  int ncores;
  int zend;
  int nlastplan;
  int edge;
  int TrackMultiplicity;
  std::vector<int> Ivec;
  std::vector<int> Jvec;
  std::vector<int> Kvec;
  std::vector<int> THvec;
  std::vector<int> vecdensity;
  std::vector<int> coresize;
  std::vector<int> corebegin;
  std::vector<int> coreend;
  std::vector<double> coreposK;
  std::vector<double> coreradius;
  std::vector<int> TrackTag;
  std::vector<int> TrackBegin;
  std::vector<int> TrackEnd;
  std::vector<int> TrackLength;
} ;

#endif  //  EVENTDISPLAYPROC_HH
