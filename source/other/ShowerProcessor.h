
#ifndef ShowerProcessor_h
#define ShowerProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <map>
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

class ShowerProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new ShowerProcessor ; }
  
  
  ShowerProcessor() ;
  
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
  
  virtual void makeHitMap();
  std::vector<int> Nhit();
  virtual int Nlayer();
  virtual int getShowerStartingLayer();
  virtual int End();
  virtual float CenterOfXGravity();
  virtual float CenterOfYGravity();
  virtual float CenterOfZGravity(int begin);
  virtual float Radius(int Zbegin, float* cog);
  virtual int NhitLastPlan();
  virtual void fillTree();
  virtual void ClearVector();
  std::vector<Density_t> Density();
  virtual int Edge();
  virtual int holeFinder(int begin, float *cog);
  virtual int hole1Finder(int begin, float *cog);
  virtual float HitinFirstPlates();
  virtual float HitinCentralCells(float *cog);
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
  std::map<int,EVENT::CalorimeterHit*> hitmap;
  std::vector<EVENT::CalorimeterHit*> calohit;
  int numElements;
  LCCollection * col;
  TFile *file;
  TTree *tree;
  TBranch *N;
  TBranch *N1;
  TBranch *N2;
  TBranch *N3;
  TBranch *NLayer;
  TBranch *ZBegin;
  TBranch *ZEnd;
  TBranch *RAdius;
  TBranch *CoGX;
  TBranch *CoGY;
  TBranch *CoGZ;
  TBranch *E;
  TBranch *NLast;
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
  TBranch *TRackBegin;
  TBranch *TRackEnd;
  TBranch *TRackLength;
  TBranch *TRackMulti;
  TBranch *LOngitudinal;
  TBranch *TRansversal;

  int nhit;
  int nhit1;
  int nhit2;
  int nhit3;
  int nlayer;
  int begin;
  int zend;
  float radius;
  int energy;
  int energy_;
  float EHratio;
  int hole;
  int hole1;
  float cog[3];
  int ncores;
  int nlastplan;
  int edge;
  float transversal;
  float longitudinal;
  int TrackMultiplicity;
  std::vector<int> vecdensity;
  std::vector<int> coresize;
  std::vector<int> corebegin;
  std::vector<int> coreend;
  std::vector<double> coreposK;
  std::vector<double> coreradius;
  std::vector<int> TrackBegin;
  std::vector<int> TrackEnd;
  std::vector<int> TrackLength;

} ;

#endif



