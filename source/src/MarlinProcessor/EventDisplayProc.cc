
#include "MarlinProcessor/EventDisplayProc.hh"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
#include <TROOT.h>
#include <TMath.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;
using namespace std;

#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

EventDisplayProc aEventDisplayProc ;


EventDisplayProc::EventDisplayProc() : Processor("EventDisplayProc") {

  _description = "EventDisplayProc calculates shower variable" ;
  
  // register steering parameters: name, description, class-variable, default value
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALBarrel"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			   "HCALCollections" , 
			   "HCAL Collection Names"  ,
			   _hcalCollections  ,
			   hcalCollections);

  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationCalorimeterHit")) ; 
    
  registerProcessorParameter( "Energy" ,
			      "Pion energy",
			      energy_,
			      (int) 0 ); 

  
  std::vector<float> thresholdHcal;
  thresholdHcal.push_back(0.114);
  thresholdHcal.push_back(1.45);
  thresholdHcal.push_back(3.80);
  registerProcessorParameter("HCALThreshold" , 
  			       "Threshold for HCAL Hits in GeV" ,
  			       _thresholdHcal,
  			       thresholdHcal);

  registerProcessorParameter( "RootFileName" ,
			      "Name of the ROOT file where tree is stored",
			      treeFileName_,
			      std::string("showers.root") ); 

  std::vector<int> hcalLayers;
  hcalLayers.push_back(48);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);
}

void EventDisplayProc::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    //fg: need to set default encoding in for reading old files...
    UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");
    file = new TFile(treeFileName_.c_str(),"recreate");
    
    tree = (TTree*)file->Get("tree");
    if(!tree){
      std::cout << "tree creation" << std::endl; 
      tree = new TTree("tree","Showers");
    }
    N = tree->Branch("Nhit",&nhit);
    ZBegin = tree->Branch("Zbegin",&begin);
    RAdius = tree->Branch("Radius",&radius);
    NLay = tree->Branch("Nlayer",&nlayer);
    E = tree->Branch("Ebeam",&energy);
    I = tree->Branch("I","std::vector<int>",&Ivec);
    J = tree->Branch("J","std::vector<int>",&Jvec);
    K = tree->Branch("K","std::vector<int>",&Kvec);
    TH = tree->Branch("fThr","std::vector<int>",&THvec);
    DEns = tree->Branch("Density","std::vector<int>",&vecdensity);
    EH = tree->Branch("EHratio",&EHratio);
    HOle = tree->Branch("Hole",&hole);
    HOle1=tree->Branch("Hole1",&hole1);
    COG = tree->Branch("CoG",&cog,"CoG[3]/F");
    NCores=tree->Branch("Ncores",&ncores);
    COreSIze=tree->Branch("CoreSize","std::vector<int>",&coresize);
    COreBegin=tree->Branch("CoreBegin","std::vector<int>",&corebegin);
    COreEnd=tree->Branch("CoreEnd","std::vector<int>",&coreend);
    COrePosK=tree->Branch("CorePos","std::vector<double>",&coreposK);
    COreRadius=tree->Branch("CoreRadius","std::vector<double>",&coreradius);
    TRackTag = tree->Branch("TrackTag","std::vector<int>",&TrackTag);
    TRackBegin = tree->Branch("TrackBegin","std::vector<int>",&TrackBegin);
    TRackEnd = tree->Branch("TrackEnd","std::vector<int>",&TrackEnd);
    TRackLength = tree->Branch("TrackLength","std::vector<int>",&TrackLength);
    TRackMulti = tree->Branch("TrackMultiplicity",&TrackMultiplicity);
}

void EventDisplayProc::ClearVector()
{
  Ivec.clear(); 
  Jvec.clear(); 
  Kvec.clear(); 
  THvec.clear();
  vecdensity.clear();
  coresize.clear();
  corebegin.clear();
  coreend.clear();
  coreposK.clear();
  coreradius.clear();
  TrackTag.clear();
  TrackBegin.clear(); 
  TrackEnd.clear();
  TrackLength.clear();
}

void EventDisplayProc::fillTree()
{
  file->cd();
  tree->Fill();
}

std::vector<int> EventDisplayProc::Nhit()
{
  std::vector<int> _nhit;
  int Nhit1 = 0; 
  int Nhit2 = 0; 
  int Nhit3 = 0;
  for (int j=0; j < numElements; ++j) {
    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
    float fThr = hit->getEnergy();
    if(fThr>0.0 && fThr<1.2) Nhit2++;
    if(fThr>1.1 && fThr<2.2) Nhit1++;
    if(fThr>2.2) Nhit3++;
    
  }
  _nhit.push_back(numElements);
  _nhit.push_back(Nhit1);
  _nhit.push_back(Nhit2);
  _nhit.push_back(Nhit3);
  return _nhit;
}

int EventDisplayProc::Nlayer()
{
  int nlayer = 0;
  UTIL::CellIDDecoder<IMPL::CalorimeterHitImpl> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int iK=0; iK<48; iK++){
    for (int j=0; j<numElements; ++j){
      IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
      int layer=decoder(hit)["K-1"];
      if(iK==layer) { nlayer++; break; }
    }
  }
  return nlayer;
}

int EventDisplayProc::Begin()
{
  UTIL::CellIDDecoder<IMPL::CalorimeterHitImpl> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int Zbegin = -1;
  for(Int_t iK=1; iK<=48; iK++){
    int count = 0;
    for (int j=0; j<numElements; ++j){
      IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
      int Layer=decoder(hit)["K-1"]+1;
      if(iK==Layer) count++;
    }
    if(count<5) continue;
    int count2 = 0;
    for (int j=0; j<numElements; ++j){
      IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
      int Layer=decoder(hit)["K-1"]+1;
      if(iK+1==Layer || iK+2==Layer || iK+3==Layer) count2++;
    }
    if(count2>15) {Zbegin=iK; break;}
  }
  return Zbegin;
}

int EventDisplayProc::End()
{
  UTIL::CellIDDecoder<IMPL::CalorimeterHitImpl> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int end = 49;
  for(Int_t iK=1; iK<=48; iK++){
    int count = 0;
    for (int j=0; j<numElements; ++j){
      IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
      int Layer=decoder(hit)["K-1"]+1;
      if((49-iK)==Layer) count++;
    }
    if(count>15) {end=49-iK; break;}
    if(count<5) continue;
    int count2 = 0;
    for (int j=0; j<numElements; ++j){
      IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
      int Layer=decoder(hit)["K-1"]+1;
      if((49-iK-1)==Layer || (49-iK-2)==Layer || (49-iK-3)==Layer) count2++;
    }
    if(count2>15) {end=49-iK; break;}
  }
  return end;
}

float EventDisplayProc::CenterOfXGravity()
{
  float weight = 0;
  float cog = 0;
  float sumweight = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
    int Icell=idDecoder(hit)["I"];
    float fThr = hit->getEnergy();
    if(fThr>=2.5)      weight = 1;
    else if(fThr>=1.5) weight = 1;
    else if(fThr>=0.5) weight = 1;
    cog += weight*Icell;
    sumweight += weight;
  }
  return cog/sumweight;
}


float EventDisplayProc::CenterOfYGravity()
{
  float weight = 0;
  float cog = 0;
  float sumweight = 0;
  UTIL::CellIDDecoder<IMPL::CalorimeterHitImpl> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
    int Jcell=decoder(hit)["J"];
    float fThr = hit->getEnergy();
    if(fThr>=2.5)      weight = 1;
    else if(fThr>=1.5) weight = 1;
    else if(fThr>=0.5) weight = 1;
    cog += weight*Jcell;
    sumweight += weight;
  }   
  return cog/sumweight;
}

float EventDisplayProc::CenterOfZGravity(int begin, int end)
{
  float weight = 0;
  float cog = 0;
  float sumweight=0;
  UTIL::CellIDDecoder<IMPL::CalorimeterHitImpl> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
    int Layer=decoder(hit)["K-1"]+1;
    if(Layer<begin || Layer>end) continue;
    float fThr = hit->getEnergy();
    if(fThr>=2.5)      weight = 3;
    else if(fThr>=1.5) weight = 2;
    else if(fThr>=0.5) weight = 1;
    cog += weight*Layer;
    sumweight += weight;
  }
  return cog/sumweight;
}

float EventDisplayProc::Radius(int Zbegin, float* cog)
{
  float radius = 0;
  int compteur = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
    if(idDecoder(hit)["K-1"]>=Zbegin){
      int I = idDecoder(hit)["I"];
      int J = idDecoder(hit)["J"];
      radius=radius+((cog[0]-I)*(cog[0]-I)+
		     (cog[1]-J)*(cog[1]-J));
      compteur++;
    }
  }
  if(compteur==0) return 0;
  return TMath::Sqrt(radius/compteur);
}

int EventDisplayProc::NhitLastPlan()
{
  int nlastplan = 0;
  UTIL::CellIDDecoder<IMPL::CalorimeterHitImpl> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int iK=42; iK<48; iK++){
    for (int j=0; j<numElements; ++j){
      IMPL::CalorimeterHitImpl *hit = dynamic_cast<IMPL::CalorimeterHitImpl*>( col->getElementAt( j ) );
      int layer=decoder(hit)["K-1"]+1;
      if(iK==layer) nlastplan++;
    }
  }
  return nlastplan;
}

float getSimClusterPosition(int i,Cluster clus)
{
  size_t size=clus.cluster.size();
  float pos=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=clus.cluster.begin(); it!=clus.cluster.end(); ++it){
    if(i==0) pos+=idDecoder(*it)["I"];
    if(i==1) pos+=idDecoder(*it)["J"];
    if(i==2) pos+=idDecoder(*it)["K-1"];
  }
  return pos/size;
}

bool SimNeighbourhood(Cluster clus,std::vector<Cluster> clVec)
{
  bool isol=false;
  int neighbour=0;
  int neighbourbis=0;
  int neighbourbisbis=0;
  float X=clus.pos[0];
  float Y=clus.pos[1];
  float Z=clus.pos[2];
  for(std::vector<Cluster>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
    float x=(*jt).pos[0];
    float y=(*jt).pos[1];
    float z=(*jt).pos[2];
    if(abs(X-x)>1 && abs(X-x)<10 &&
       abs(Y-y)>1 && abs(Y-y)<10 &&
       Z==z) 
      neighbour++;
    //    if((*jt).cluster.size()>10) continue;
    if(Z!=z && abs(Z-z)<3 &&
       abs(X-x)<5 &&
       abs(Y-y)<5)
      neighbourbis++;
    if(Z!=z && abs(Z-z)<3 &&
       abs(X-x)<5 &&
       abs(Y-y)<5
       && (*jt).cluster.size()>8)
      neighbourbisbis++;
  }
  if(neighbour<2&&neighbourbis>0&&neighbourbis<6&&neighbourbisbis<2) isol=true;
  return isol;
}

std::vector<Cluster> EventDisplayProc::GetClusters()
{
  std::vector<Cluster> vec;
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    std::vector<EVENT::CalorimeterHit*> _hit;
    EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
    if(std::find(_temp.begin(),_temp.end(), hit)!=_temp.end()) continue;
    bool append=false;
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=_temp.begin(); it!=_temp.end(); it++){
      if(hit==*it){
	append=true; break; 
      } 
    }
    if(append==true) continue;
    _temp.push_back(hit);
    _hit.push_back(hit);
    for(int i=0; i<numElements; ++i){
      EVENT::CalorimeterHit *hit2=dynamic_cast<EVENT::CalorimeterHit*>(col->getElementAt(i));
      if(idDecoder(hit2)["K-1"]!=idDecoder(hit)["K-1"]) continue;
      append=false;
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=_temp.begin(); it!=_temp.end(); ++it){
	if(hit2==*it){ 
	  append=true; 
	  break;
	}
      }
      if(append==true) continue;
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=_hit.begin(); it!=_hit.end(); ++it){
	if(abs(idDecoder(*it)["I"]-idDecoder(hit2)["I"])<2 &&
	   abs(idDecoder(*it)["J"]-idDecoder(hit2)["J"])<2){ 
	  _hit.push_back(hit2);
	  _temp.push_back(hit2);
	  break;
	}
      }
    }
    Cluster cl;
    cl.cluster=_hit;
    for(int i=0; i<3; i++) cl.pos[i]=getSimClusterPosition(i,cl);
    vec.push_back(cl);
  }
  return vec;
}

std::vector<bool> EventDisplayProc::TagTracks(std::vector<Cluster> clVec)
{  
  std::vector<bool> tag;

  int houghx[tetamax][140];
  int houghy[tetamax][140];
  for(size_t i=0; i<tetamax; i++){
    for(size_t j=0; j<140; j++){
      houghx[i][j]=0;
      houghy[i][j]=0;
    }
  }

  float pi=TMath::Pi();
  int ClusterVecSize=clVec.size();
  tag.reserve(ClusterVecSize);
  std::vector<bool> tagx;
  std::vector<bool> tagy;
  tagx.reserve(ClusterVecSize);
  tagy.reserve(ClusterVecSize);
  for(int i=0; i<ClusterVecSize; i++){tagx.push_back(false);tagy.push_back(false);}
  std::vector<IsATrack_t> isTrack;
  for(std::vector<Cluster>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if((*it).cluster.size()>10) continue;
    bool isol=SimNeighbourhood((*it),clVec);
    if(isol==true){
      IsATrack_t track;
      for(size_t ith=0; ith<tetamax; ith++){
	track.rhox[ith]=abs((int)((*it).pos[2]*cos(-pi/2+ith*pi/tetamax)+(*it).pos[0]*sin(-pi/2+ith*pi/tetamax)));
  	track.rhoy[ith]=abs((int)((*it).pos[2]*cos(-pi/2+ith*pi/tetamax)+(*it).pos[1]*sin(-pi/2+ith*pi/tetamax)));
	track.index=std::distance(clVec.begin(),it);
	houghx[ith][track.rhox[ith]]+=1;
      }
      track.pos[0]=(*it).pos[0];
      track.pos[1]=(*it).pos[1];
      track.pos[2]=(*it).pos[2];
      isTrack.push_back(track);
    }
  }
  for(size_t ith=0; ith<tetamax; ith++){
    for(std::vector<IsATrack_t>::iterator it=isTrack.begin(); it!=isTrack.end(); ++it){
      if(houghx[ith][(*it).rhox[ith]]<6) continue;
      tagx[(*it).index]=true;
      int nunu=0;
      for(std::vector<IsATrack_t>::iterator jt=isTrack.begin(); jt!=isTrack.end(); ++jt){
	if(abs((*jt).rhox[ith]-(*it).rhox[ith])>1 || 
	   (*jt).pos[2]==(*it).pos[2] ) continue;
	if(abs((*it).pos[0]-(*jt).pos[0])<15 &&
	   abs((*it).pos[1]-(*jt).pos[1])<15 &&
	   abs((*it).pos[2]-(*jt).pos[2])<15){
	  nunu++;
	}
      }
      if(nunu<2) continue;
      for(size_t jth=0; jth<tetamax; jth++)
	houghy[jth][(*it).rhoy[jth]]+=1;
    }
  }
  for(size_t jth=0; jth<tetamax; jth++){
    for(std::vector<IsATrack_t>::iterator jt=isTrack.begin(); jt!=isTrack.end(); ++jt){
      if(houghy[jth][(*jt).rhoy[jth]]<6) continue;
      tagy[(*jt).index]=true;
    }
  }
  for(std::vector<Cluster>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if(tagx[std::distance(clVec.begin(),it)]==true &&
       tagy[std::distance(clVec.begin(),it)]==true) tag.push_back(true);
    else tag.push_back(false);
  }
  //  std::vector<bool> tag2=getRealTrack(tag1,clVec);
  return tag;
}

std::vector<Track*> EventDisplayProc::getRealTrack(std::vector<bool> Tag,std::vector<Cluster> clVec)
{
  tagtrack.clear();
  std::vector<Track*> TrackVec;
//  std::vector<bool> tag;
  tagtrack.reserve(clVec.size());
  for(size_t i=0; i<clVec.size(); i++) tagtrack.push_back(false);
  std::vector<Cluster> _temp;
  for(std::vector<Cluster>::iterator it=clVec.begin(); it!=clVec.end(); ++it){
    if(Tag[std::distance(clVec.begin(),it)]==false) continue;
        bool append=false;
    for(std::vector<Cluster>::iterator kt=_temp.begin(); kt!=_temp.end(); ++kt){
      if((*kt).cluster==(*it).cluster){
	append=true;
	break;
      }
    }
    if(append==true) continue;
    std::vector<Cluster> aTrack;
    std::vector<size_t> index;
    _temp.push_back(*it);
    aTrack.push_back(*it);
    index.push_back(std::distance(clVec.begin(), it));
    for(std::vector<Cluster>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt){
      if(Tag[std::distance(clVec.begin(),jt)]==false) continue;
      append=false;
      for(std::vector<Cluster>::iterator kt=_temp.begin(); kt!=_temp.end(); ++kt){
      	if((*jt).cluster==(*kt).cluster){
	  append=true;
	  break;
	}
      }
      if(append==true) continue;
      for(std::vector<Cluster>::iterator kt=aTrack.begin(); kt!=aTrack.end(); ++kt){
	if(abs((*kt).pos[2]-(*jt).pos[2])<=2 &&
	   abs((*kt).pos[1]-(*jt).pos[1])<=2 &&
	   abs((*kt).pos[0]-(*jt).pos[0])<=2){
	  aTrack.push_back(*jt);
	  _temp.push_back(*jt);
	  index.push_back(std::distance(clVec.begin(), jt));
	  break;
	}
      }
    }
    Track *track=new Track();
    track->TrackClus=aTrack;
    track->getTrackStartingPoint();
    track->getLastTrackPoint();
    int length=track->end - track->begin;
//    float Xangle=track->getTrackAngle(0);
    if(int(track->TrackClus.size())<4) continue;
    if(length<=2) continue;
    float angle=track->getTrackAngle();
    std::cout << "Find one track with " << track->TrackClus.size() << " cluster, starting at " << track->begin << ", and ending at " << track->end << ", with an angle of: " << angle*360/(2*TMath::Pi()) << std::endl;
    for(std::vector<size_t>::iterator kt=index.begin(); kt!=index.end(); ++kt){
      tagtrack[*kt]=true;
    }
    TrackVec.push_back(track);
  }
  return TrackVec;
}

void Track::getTrackStartingPoint()
{
  std::vector<int> K;
  for(std::vector<Cluster>::iterator it=this->TrackClus.begin(); it!=this->TrackClus.end(); ++it){
    K.push_back(int((*it).pos[2]));
  }
  begin=*std::min_element(K.begin(),K.end());
}

void Track::getLastTrackPoint()
{
  std::vector<int> K;
  for(std::vector<Cluster>::iterator it=this->TrackClus.begin(); it!=this->TrackClus.end(); ++it){
    K.push_back(int((*it).pos[2]));
  }
  end=*std::max_element(K.begin(),K.end());
}

float Track::getTrackAngle()
{
  std::vector<float> I;
  std::vector<float> J;
  std::vector<float> K;
  for(std::vector<Cluster>::iterator it=this->TrackClus.begin(); it!=this->TrackClus.end(); ++it){
    K.push_back((*it).pos[2]);
    J.push_back((*it).pos[1]);
    I.push_back((*it).pos[0]);
  }
  
  float deltaX=I[std::distance(K.begin(),std::max_element(K.begin(),K.end()))]-I[std::distance(K.begin(),std::min_element(K.begin(),K.end()))];
  float deltaY=J[std::distance(K.begin(),std::max_element(K.begin(),K.end()))]-J[std::distance(K.begin(),std::min_element(K.begin(),K.end()))];
  float deltaZ=*std::max_element(K.begin(),K.end())-*std::min_element(K.begin(),K.end());
  
  float angle=sqrt(deltaX*deltaX+deltaY*deltaY)/deltaZ;
  return angle;
}

//float Track::getTrackAngle(int i)
////i=0,1 for (z,x) and (z,y) plan
//{
//  TH2* h=new TH2D("h","",this->end,this->begin,this->end,96,0,96);
//  for(std::vector<Cluster>::iterator it=this->TrackClus.begin(); it!=this->TrackClus.end(); ++it){
//    h->Fill((*it).pos[2],(*it).pos[i]);
//  }
//  TF1* func = new TF1("func","pol1",this->begin,this->end);
//  h->Fit(func);
//  float a= func->GetParameter(1);
//  float angle=-1/atan(a);
//  delete h;
//  delete func;
//  return angle*360/(2*TMath::Pi());
//}

std::vector<MyHit_t> EventDisplayProc::getIJK()
{
  std::vector<MyHit_t> vec;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
    MyHit_t myhit;
    myhit.I=idDecoder(hit)["I"];
    myhit.J=idDecoder(hit)["J"];
    myhit.K=idDecoder(hit)["K-1"];
    myhit.TH=int(hit->getEnergy());
    vec.push_back(myhit);
  }
  return vec;
}

std::vector<int> EventDisplayProc::HitInTrack(std::vector<bool> Tag,std::vector<Cluster> clVec)
{
  std::vector<int> taghit;
  int size=clVec.size();
  for(int i=0; i<size; i++){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=clVec[i].cluster.begin(); it!=clVec[i].cluster.end(); ++it){
      if(Tag[i]==true) taghit.push_back(1);
      else taghit.push_back(0);
    }
  }
  return taghit;
}

std::vector<Density_t> EventDisplayProc::Density()
{
  std::vector<Density_t> vec;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    Density_t dens;
    EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
    int Icell=idDecoder(hit)["I"];
    int Jcell=idDecoder(hit)["J"];
    int Kcell=idDecoder(hit)["K-1"];
    dens.pos[0]=Icell;
    dens.pos[1]=Jcell;
    dens.pos[2]=Kcell;
    dens.dens2D=0;
    dens.dens3D=0;
    for (int j2=0; j2<numElements; ++j2){
      EVENT::CalorimeterHit *hit2 = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j2 ) );
      if(j2==j) continue;
      int Icell2=idDecoder(hit2)["I"];
      int Jcell2=idDecoder(hit2)["J"];
      int Kcell2=idDecoder(hit2)["K-1"];
      //it can become useful
      int weight=0;
      float fThr = hit2->getEnergy();
      if(fThr>=2.5)      weight = 10;
      else if(fThr>=1.5) weight = 1;
      else if(fThr>=0.5) weight = 5;
      if(abs(Icell-Icell2)<2 && 
	 abs(Jcell-Jcell2)<2 &&
	 abs(Kcell-Kcell2)<2){
	dens.dens3D += weight;
      }
      if(abs(Icell-Icell2)<3 && 
	 abs(Jcell-Jcell2)<3 &&
	 Kcell==Kcell2){
	dens.dens2D += 1;
      }
    }
    vec.push_back(dens);
  }
  return vec;
}

int EventDisplayProc::Edge()
{
  int edge=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for (int j=0; j<numElements; ++j){
    EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
    if( idDecoder(hit)["I"]<10 || 
	idDecoder(hit)["I"]>90 ||
	idDecoder(hit)["J"]<10 || 
	idDecoder(hit)["J"]>90  ){
      edge+=1;
    }
  }
  //  std::cout << edge << std::endl;
  return edge;
}

int EventDisplayProc::hole1Finder(int begin, float *cog)
{
  if(begin<0) return 0;
  int hole=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int k=0; k<begin-1; k++){
    int nhit=0;
    for (int j=0; j<numElements; ++j){
      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
      if(idDecoder(hit)["K-1"]==begin-1-k){
	nhit++;
	if(abs(idDecoder(hit)["I"]-cog[0])>10 ||
	   abs(idDecoder(hit)["J"]-cog[1])>10 ) nhit--;
      }
    }
    if(nhit==0) hole++;
  }
  return hole;
}

int EventDisplayProc::holeFinder(int begin, float *cog)
{
  if(begin<0) return 0;
  int hole=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int k=begin+1; k<begin+8; k++){
    int nhit=0;
    for (int j=0; j<numElements; ++j){
      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
      if(idDecoder(hit)["K-1"]==k){
	nhit++;
	if(abs(idDecoder(hit)["I"]-cog[0])>10 ||
	   abs(idDecoder(hit)["J"]-cog[1])>10 ) nhit--;
      }
    }
    if(nhit==0) hole++;
  }
  return hole;
}

std::vector<MyCore_t> EventDisplayProc::CoreFinder(std::vector<Density_t> density)
{
  std::vector<MyCore_t> vec;
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int k=0; k<50; k++){
    for (int j=0; j<numElements; ++j){
      if(density[j].dens3D<30) continue;
      std::vector<EVENT::CalorimeterHit*> _hit;
      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
      if(idDecoder(hit)["K-1"]!=k) continue;
      bool append=false;
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=_temp.begin(); it!=_temp.end(); it++){
	if(hit==*it){
	  append=true; break; 
	} 
      }
      if(append==true) continue;
      _temp.push_back(hit);
      _hit.push_back(hit);
      for(int K=0; K<50; K++){
	for(int i=0; i<numElements; ++i){
	  if(density[i].dens3D<30) continue;
	  EVENT::CalorimeterHit *hit2=dynamic_cast<EVENT::CalorimeterHit*>(col->getElementAt(i));
	  if(idDecoder(hit2)["K-1"]!=K) continue;
	  append=false;
	  for(std::vector<EVENT::CalorimeterHit*>::iterator it=_temp.begin(); it!=_temp.end(); ++it){
	    if(hit2==*it){ 
	      append=true; 
	      break;
	    }
	  }
	  if(append==true) continue;
	  for(std::vector<EVENT::CalorimeterHit*>::iterator it=_hit.begin(); it!=_hit.end(); ++it){
	    if(abs(idDecoder(*it)["K-1"]-idDecoder(hit2)["K-1"])<=1){ 
	      _hit.push_back(hit2);
	      _temp.push_back(hit2);
	      break;
	    }
	  }
	}
      }
      MyCore_t core;
      core.hits=_hit;
      core.size=CoreSize(core);
      core.begin=CoreBegin(core);
      core.end=CoreEnd(core);
      core.pos=getCorePosition(core);
      core.radius=getCoreRadius(core);
      vec.push_back(core);
    }
  }
  return vec;
}

int EventDisplayProc::NumberOfCores(std::vector<MyCore_t> vec)
{
  return vec.size();
}

int EventDisplayProc::CoreSize(MyCore_t core)
{
  return core.hits.size();
}

int EventDisplayProc::CoreBegin(MyCore_t core)
{
  int min=50;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    if(idDecoder(*it)["K-1"]<min) min=idDecoder(*it)["K-1"];
  }
  return min;
}

int EventDisplayProc::CoreEnd(MyCore_t core)
{
  int max=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    if(idDecoder(*it)["K-1"]>max) max=idDecoder(*it)["K-1"];
  }
  return max;
}

std::vector<float> EventDisplayProc::getCorePosition(MyCore_t core)
{
  float pos[3];
  for(size_t i=0; i<3; i++) pos[i]=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    pos[0]+=idDecoder(*it)["I"];
    pos[1]+=idDecoder(*it)["J"];
    pos[2]+=idDecoder(*it)["K-1"];
  }
  std::vector<float> vec;
  for(size_t i=0; i<3; i++) vec.push_back(pos[i]/core.hits.size());
  return vec;
}

float EventDisplayProc::getCoreRadius(MyCore_t core)
{
  float radius=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    radius+=(idDecoder(*it)["I"]-core.pos[0])*(idDecoder(*it)["I"]-core.pos[0])+
            (idDecoder(*it)["J"]-core.pos[1])*(idDecoder(*it)["J"]-core.pos[1]);
  }
  return sqrt(radius/CoreSize(core));
}

void EventDisplayProc::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void EventDisplayProc::processEvent( LCEvent * evt )
{   
  //
  // * Reading HCAL Collections of CalorimeterHits* 
  //
  ClearVector();
  std::string initString;
  for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
    std::string colName =  _hcalCollections[i] ;
    try{
      col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      nlayer = Nlayer();
      begin = Begin();
      zend = End();
      cog[0]=CenterOfXGravity();
      cog[1]=CenterOfYGravity();
      cog[2]=CenterOfZGravity(begin, zend);
      radius = Radius(begin,cog);
      nlastplan=NhitLastPlan();
      hole=holeFinder(begin,cog);

      edge=Edge();
      if(zend<100 /*&& begin>0 && 
	 zend-begin>2 &&
	 float(nlastplan)/numElements<15 &&
	 radius/cog[2]>0.30 &&
	 float(edge)/numElements<0.1*/){
	nhit=numElements;
	hole1=hole1Finder(begin,cog);
	std::vector<MyHit_t> vec=getIJK();
	int size=vec.size();
	Ivec.reserve(size);
	Jvec.reserve(size);
	Kvec.reserve(size);
	THvec.reserve(size);
	for(int j=0; j<size; j++){
	  Ivec.push_back(vec[j].I);
	  Jvec.push_back(vec[j].J);
	  Kvec.push_back(vec[j].K);
	  THvec.push_back(vec[j].TH);
	}
	std::vector<Density_t> _density=Density();
	size=_density.size();
	vecdensity.reserve(size);
	float ratio=0;
	for(int i=0; i<size; i++){
	  vecdensity.push_back(_density[i].dens3D);
	  if(_density[i].dens3D>30) ratio++;
	}
	std::vector<MyCore_t> core=CoreFinder(_density);
	ncores=NumberOfCores(core);
	for(int i=0; i<ncores; i++){
	  coresize.push_back(core[i].size);
	  corebegin.push_back(core[i].begin);
	  coreend.push_back(core[i].end);
	  coreposK.push_back(core[i].pos[2]);
	  coreradius.push_back(core[i].radius);
	}
	std::vector<Cluster> clVec=GetClusters();
	TrackTag.reserve(size);
	std::vector<bool> Tag=TagTracks(clVec);
	std::vector<Track*> TrackVec=getRealTrack(Tag,clVec);
	TrackBegin.reserve(TrackVec.size());
	TrackEnd.reserve(TrackVec.size());
	TrackLength.reserve(TrackVec.size());
	TrackMultiplicity=int(TrackVec.size());
	for(size_t i=0; i<TrackVec.size(); i++){
	  TrackBegin.push_back(TrackVec[i]->begin);
	  TrackEnd.push_back(TrackVec[i]->end);
	  TrackLength.push_back(TrackVec[i]->end-TrackVec[i]->begin);
	}
	TrackTag=HitInTrack(tagtrack, clVec);
	if(int(TrackTag.size())!=numElements) std::cout << "WARNING too much hits in track vector" << std::endl; 
	EHratio=ratio/numElements;
	fillTree();
      }
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void EventDisplayProc::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EventDisplayProc::end(){ 
  //  file->cd();

  file->Write();
  file->Close();
  std::cout << "EventDisplayProc::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}

