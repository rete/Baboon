
#include "ShowerProcessor.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
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

#include <map>
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<double> +;
#endif

ShowerProcessor aShowerProcessor ;


ShowerProcessor::ShowerProcessor() : Processor("ShowerProcessor") {

  // modify processor description
  _description = "ShowerProcessor calculates shower variable" ;
  
  
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("HCALBarrel"));

  // register steering parameters: name, description, class-variable, default value
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
  
  registerProcessorParameter( "RootFileName" ,
			      "Name of the ROOT file where tree is stored",
			      treeFileName_,
			      std::string("showers.root") ); 
  
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


  std::vector<int> hcalLayers;
  hcalLayers.push_back(48);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);
}



void ShowerProcessor::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    //fg: need to set default encoding in for reading old files...
    UTIL::CellIDDecoder<CalorimeterHit*>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

    file = new TFile(treeFileName_.c_str(),"RECREATE");
    
    tree = (TTree*)file->Get("tree");
    if(!tree){
      std::cout << "tree creation" << std::endl; 
      tree = new TTree("tree","Shower variables");
    }
    N = tree->Branch("Nhit",&nhit);
    N1 = tree->Branch("Nhit1",&nhit1);
    N2 = tree->Branch("Nhit2",&nhit2);
    N3 = tree->Branch("Nhit3",&nhit3);
    NLayer = tree->Branch("Nlayer",&nlayer);
    ZBegin = tree->Branch("Begin",&begin);
    ZEnd = tree->Branch("End",&zend);
    RAdius = tree->Branch("Radius",&radius);
    E = tree->Branch("Ebeam",&energy);
    NLast = tree->Branch("Lasthit",&nlastplan);
    DEns = tree->Branch("Density","std::vector<int>",&vecdensity);
    EH = tree->Branch("EHratio",&EHratio);
    HOle = tree->Branch("Hole",&hole);
    HOle1=tree->Branch("Hole1",&hole1);
    COG = tree->Branch("CoG",&cog,"CoG[3]/F");
    LOngitudinal=tree->Branch("LongitudinalCut",&longitudinal);
    TRansversal=tree->Branch("TransversalCut",&transversal);
    NCores=tree->Branch("Ncores",&ncores);
    COreSIze=tree->Branch("CoreSize","std::vector<int>",&coresize);
    COreBegin=tree->Branch("CoreBegin","std::vector<int>",&corebegin);
    COreEnd=tree->Branch("CoreEnd","std::vector<int>",&coreend);
    COrePosK=tree->Branch("CorePos","std::vector<double>",&coreposK);
    COreRadius=tree->Branch("CoreRadius","std::vector<double>",&coreradius);
    TRackBegin = tree->Branch("TrackBegin","std::vector<int>",&TrackBegin);
    TRackEnd = tree->Branch("TrackEnd","std::vector<int>",&TrackEnd);
    TRackLength = tree->Branch("TrackLength","std::vector<int>",&TrackLength);
    TRackMulti = tree->Branch("TrackMultiplicity",&TrackMultiplicity);
}

void ShowerProcessor::ClearVector()
{
  hitmap.clear();
  calohit.clear();
  vecdensity.clear();
  coresize.clear();
  corebegin.clear();
  coreend.clear();
  coreposK.clear();
  coreradius.clear();
  TrackBegin.clear(); 
  TrackEnd.clear();
  TrackLength.clear();
}

void ShowerProcessor::fillTree()
{
  file->cd();
  tree->Fill();   
}

int IJKToKey(int i, int j, int k){return 100*100*k+100*j+i;}



void ShowerProcessor::makeHitMap()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
  for (int j=0; j < numElements; ++j) {
    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
    int key=IJKToKey(IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]);
    hitmap[key]=hit;
  }
}

std::vector<int> ShowerProcessor::Nhit()
{
  std::vector<int> _nhit;
  int Nhit1 = 0; 
  int Nhit2 = 0; 
  int Nhit3 = 0;
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    float fThr = (*it)->getEnergy();
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

int ShowerProcessor::Nlayer()
{
  int nlayer = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int iK=0; iK<50; iK++){
    for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
      int layer=idDecoder(*it)["K-1"];
      if(iK==layer) { nlayer++; break; }
    }
  }
  return nlayer;
}

int ShowerProcessor::getShowerStartingLayer()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> decoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int Zbegin = -10000;
  float xmoy=CenterOfXGravity();
  float ymoy=CenterOfYGravity();
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    uint I=decoder(*it)["I"];
    uint J=decoder(*it)["J"];
    uint K=decoder(*it)["K-1"];
    if(abs(I-xmoy)>5||abs(J-ymoy)>5) continue;
    int count=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=calohit.begin(); jt!=calohit.end(); ++jt){
      if(jt==it) continue;
      if(K==decoder(*jt)["K-1"] &&
	 abs(decoder(*jt)["I"]-xmoy)<5 &&
	 abs(decoder(*jt)["J"]-ymoy)<5 ) 
	count++;
    }
    if(count<5) continue;
    int count2[3];
    for(int j=0; j<3; j++) count2[j]=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator jt=calohit.begin(); jt!=calohit.end(); ++jt){
      if(abs(decoder(*jt)["I"]-xmoy)>5 ||
	 abs(decoder(*jt)["J"]-ymoy)>5 ) continue;
      if(decoder(*jt)["K-1"]==K+1) count2[0]++;
      if(decoder(*jt)["K-1"]==K+2) count2[1]++;
      if(decoder(*jt)["K-1"]==K+3) count2[2]++;
    }
    if(count2[0]>3&&count2[1]>3&&count2[2]>3) {
      Zbegin=K;
      break;
    }
  }
  return Zbegin;
}
//int ShowerProcessor::Begin()
//{
//  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
//  int Zbegin = -1;
//  for(Int_t iK=0; iK<=50; iK++){
//    int count = 0;
//    for (int j=0; j<numElements; ++j){
//      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
//      int Layer=idDecoder(hit)["K-1"];
//      if(iK==Layer) count++;
//    }
//    if(count<5) continue;
//    int count2 = 0;
//    for (int j=0; j<numElements; ++j){
//      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
//      int Layer=idDecoder(hit)["K-1"];
//      if(iK+1==Layer || iK+2==Layer || iK+3==Layer) count2++;
//    }
//    if(count2>15) {Zbegin=iK+1; break;}
//  }
//  return Zbegin;
//}

int ShowerProcessor::End()
{
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  int end = 100;
  for(Int_t iK=0; iK<=48; iK++){
    int count = 0;
    for (int j=0; j<numElements; ++j){
      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
      int Layer=idDecoder(hit)["K-1"]+1;
      if((48-iK)==Layer) count++;
    }
    if(count>15) {end=48-iK; break;}
    if(count<5) continue;
    int count2 = 0;
    for (int j=0; j<numElements; ++j){
      EVENT::CalorimeterHit *hit = dynamic_cast<EVENT::CalorimeterHit*>( col->getElementAt( j ) );
      int Layer=idDecoder(hit)["K-1"]+1;
      if((48-iK-1)==Layer || (48-iK-2)==Layer || (48-iK-3)==Layer) count2++;
    }
    if(count2>15) {end=48-iK; break;}
  }
  return end;
}

float ShowerProcessor::CenterOfXGravity()
{
  float cog = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    cog +=idDecoder(*it)["I"];
  }
  return cog/numElements;
}


float ShowerProcessor::CenterOfYGravity()
{
  float cog = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    cog +=idDecoder(*it)["J"];
  }
  return cog/numElements;
}

float ShowerProcessor::CenterOfZGravity(int begin)
{
  float weight = 0;
  float cog = 0;
  float sumweight=0;
  if(begin<0) return 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(idDecoder(*it)["K-1"]<begin) continue;
    float fThr = (*it)->getEnergy();
    if(fThr>=2.5)      weight = 3;
    else if(fThr>=1.5) weight = 1;
    else if(fThr>=0.5) weight = 2;
    cog += weight*(idDecoder(*it)["K-1"]-begin);
    sumweight += weight;
  }
  return cog/sumweight;
}

float ShowerProcessor::Radius(int Zbegin, float* cog)
{
  float radius = 0;
  int compteur = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(idDecoder(*it)["K-1"]>=Zbegin){
      int I = idDecoder(*it)["I"];
      int J = idDecoder(*it)["J"];
      radius=radius+((cog[0]-I)*(cog[0]-I)+
		     (cog[1]-J)*(cog[1]-J));
      compteur++;
    }
  }
  if(compteur==0) return 0;
  return TMath::Sqrt(radius/compteur);
}

int ShowerProcessor::NhitLastPlan()
{
  int nlastplan = 0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(idDecoder(*it)["K-1"]>=41) nlastplan++;
  }
  return nlastplan;
}

float ShowerProcessor::HitinFirstPlates()
{
  float nhit=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(idDecoder(*it)["K-1"]<15)
      nhit++;
  }
  return nhit;
}


float ShowerProcessor::HitinCentralCells(float *cog)
{
  float nhit=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=calohit.begin(); it!=calohit.end(); ++it){
    if(abs(idDecoder(*it)["I"]-cog[0])<7&&
       abs(idDecoder(*it)["J"]-cog[1])<7)
      nhit++;
  }
  return nhit;
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
    if((*jt).cluster==clus.cluster) continue;
    float x=(*jt).pos[0];
    float y=(*jt).pos[1];
    float z=(*jt).pos[2];
    if(abs(X-x)<10 &&
       abs(Y-y)<10 &&
       Z==z) 
      neighbour++;
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

std::vector<Cluster> ShowerProcessor::GetClusters()
{
	std::vector<Cluster> vec;
	std::vector<EVENT::CalorimeterHit*> _temp;
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");

	for(std::map<int,EVENT::CalorimeterHit*>::iterator hit=hitmap.begin(); hit!=hitmap.end(); ++hit){

		std::vector<EVENT::CalorimeterHit*> _hit;
		if(std::find(_temp.begin(),_temp.end(), hit->second)!=_temp.end()) continue;
		_temp.push_back(hit->second);
		_hit.push_back(hit->second);

		for(std::map<int,EVENT::CalorimeterHit*>::iterator hit2=hitmap.begin(); hit2!=hitmap.end(); ++hit2){

			if(std::find(_temp.begin(),_temp.end(), hit2->second)!=_temp.end()) continue;
			if(idDecoder(hit2->second)["K-1"]!=idDecoder(hit->second)["K-1"]) continue;

			for(std::vector<EVENT::CalorimeterHit*>::iterator it=_hit.begin(); it!=_hit.end(); ++it){
				if(abs(idDecoder(*it)["I"]-idDecoder(hit2->second)["I"])<2 &&
				   abs(idDecoder(*it)["J"]-idDecoder(hit2->second)["J"])<2){
					_hit.push_back(hit2->second);
					_temp.push_back(hit2->second);
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

std::vector<bool> ShowerProcessor::TagTracks(std::vector<Cluster> clVec)
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

  for(int i=0; i<ClusterVecSize; i++) {

	  tagx.push_back(false);
	  tagy.push_back(false);
  }

  std::vector<IsATrack_t> isTrack;

  for(std::vector<Cluster>::iterator it=clVec.begin(); it!=clVec.end(); ++it) {

	if( (*it).cluster.size() > 8 ) continue;


	bool isol=SimNeighbourhood((*it),clVec);

	// si le cluster n'est pas isolé -> on prends pas
	if(isol==true) {
		IsATrack_t track;

		// rempli les r et theta pour chaque cluster
		for(size_t ith=0; ith<tetamax; ith++) {

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


  for(size_t ith=0; ith<tetamax; ith++) {
	  // itere sur les cluster candidats
	  for(std::vector<IsATrack_t>::iterator it=isTrack.begin(); it!=isTrack.end(); ++it) {

		  // si le binning de la matrice de hough dans la direction x < 6 -> on prends pas
		  if( houghx[ith][(*it).rhox[ith]] < 6 ) continue;
		  // sinon on tag comme etant bon dans la direction x
		  tagx[(*it).index]=true;


		  int nunu=0;

		  for(std::vector<IsATrack_t>::iterator jt=isTrack.begin(); jt!=isTrack.end(); ++jt) {

			  // si les deux clusters sont dans le meme plan ou si
			  if(abs((*jt).rhox[ith]-(*it).rhox[ith])>1 ||
					  (*jt).pos[2]==(*it).pos[2] ) continue;

			  if(abs((*it).pos[0]-(*jt).pos[0])<15 &&
				 abs((*it).pos[1]-(*jt).pos[1])<15 &&
				 abs((*it).pos[2]-(*jt).pos[2])<15) {
				  nunu++;
			  }
		  }
		  if(nunu<2) continue;

		  for(size_t jth=0; jth<tetamax; jth++) houghy[jth][(*it).rhoy[jth]]+=1;

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
         tagy[std::distance(clVec.begin(),it)]==true) {
      tag.push_back(true);
	  }
	  else tag.push_back(false);
  }


  return tag;
}

std::vector<Track*> ShowerProcessor::getRealTrack(std::vector<bool> Tag,std::vector<Cluster> clVec)
{
	std::vector<Track*> TrackVec;
	std::vector<Cluster> _temp;

	for(std::vector<Cluster>::iterator it=clVec.begin(); it!=clVec.end(); ++it) {

		if(Tag[std::distance(clVec.begin(),it)]==false) continue;
		bool append=false;
		for(std::vector<Cluster>::iterator kt=_temp.begin(); kt!=_temp.end(); ++kt) {
			if((*kt).cluster==(*it).cluster) {
				append=true;
				break;
			}
		}
		if(append==true) continue;

		std::vector<Cluster> aTrack;
		_temp.push_back(*it);
		aTrack.push_back(*it);

		for(std::vector<Cluster>::iterator jt=clVec.begin(); jt!=clVec.end(); ++jt) {

			if(Tag[std::distance(clVec.begin(),jt)]==false) continue;
			append=false;

			for(std::vector<Cluster>::iterator kt=_temp.begin(); kt!=_temp.end(); ++kt) {
				if((*jt).cluster==(*kt).cluster){
					append=true;
					break;
				}
			}
			if(append==true) continue;

			for(std::vector<Cluster>::iterator kt=aTrack.begin(); kt!=aTrack.end(); ++kt) {
				if(abs((*kt).pos[2]-(*jt).pos[2])<=2 &&
				abs((*kt).pos[1]-(*jt).pos[1])<=4 &&
				abs((*kt).pos[0]-(*jt).pos[0])<=4){
				aTrack.push_back(*jt);
				_temp.push_back(*jt);
				break;
				}
			}
		}
		Track *track=new Track();
		track->TrackClus=aTrack;
		track->getTrackStartingPoint();
		track->getLastTrackPoint();
		int length=track->end - track->begin;
		if(int(track->TrackClus.size())<4) continue;
		if(length<=2) continue;
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
  float deltaX=I[std::distance(I.begin(),std::max_element(K.begin(),K.end()))]-I[std::distance(I.begin(),std::min_element(K.begin(),K.end()))];
  float deltaY=J[std::distance(J.begin(),std::max_element(K.begin(),K.end()))]-J[std::distance(J.begin(),std::min_element(K.begin(),K.end()))];
  float deltaZ=*std::max_element(K.begin(),K.end())-*std::min_element(K.begin(),K.end());
  float angle=sqrt(deltaX*deltaX+deltaY*deltaY)/deltaZ; 
 return angle;
}

std::vector<Density_t> ShowerProcessor::Density()
{
  std::vector<Density_t> vec;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator hit=calohit.begin(); hit!=calohit.end(); ++hit){
    Density_t dens;
    int Icell=idDecoder(*hit)["I"];
    int Jcell=idDecoder(*hit)["J"];
    int Kcell=idDecoder(*hit)["K-1"];
    dens.pos[0]=Icell;
    dens.pos[1]=Jcell;
    dens.pos[2]=Kcell;
    dens.dens2D=0;
    dens.dens3D=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator hit2=calohit.begin(); hit2!=calohit.end(); ++hit2){
      int Icell2=idDecoder(*hit2)["I"];
      int Jcell2=idDecoder(*hit2)["J"];
      int Kcell2=idDecoder(*hit2)["K-1"];
      int weight=0;
      float fThr = (*hit2)->getEnergy();
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

int ShowerProcessor::Edge()
{
  int edge=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator hit=calohit.begin(); hit!=calohit.end(); ++hit){
    if( idDecoder(*hit)["I"]<10 || 
	idDecoder(*hit)["I"]>90 ||
	idDecoder(*hit)["J"]<10 || 
	idDecoder(*hit)["J"]>90  ){
      edge+=1;
    }
  }
  return edge;
}

int ShowerProcessor::hole1Finder(int begin, float *cog)
{
  if(begin<=0) return 0;
  int hole=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int k=begin-1; k>=0; k--){
    int nhit=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator hit=calohit.begin(); hit!=calohit.end(); ++hit){
      if(idDecoder(*hit)["K-1"]==k &&
	 abs(idDecoder(*hit)["I"]-cog[0])<10&&
	 abs(idDecoder(*hit)["J"]-cog[1])<10 ) nhit++;
    }
    if(nhit==0) hole++;
  }
  return hole;
}

int ShowerProcessor::holeFinder(int begin, float *cog)
{
  if(begin<0) return 0;
  int hole=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(int k=begin+1; k<begin+8; k++){
    int nhit=0;
    for(std::vector<EVENT::CalorimeterHit*>::iterator hit=calohit.begin(); hit!=calohit.end(); ++hit){
      if(idDecoder(*hit)["K-1"]==k &&
	 abs(idDecoder(*hit)["I"]-cog[0])<10 &&
	 abs(idDecoder(*hit)["J"]-cog[1])<10 ) nhit++;
    }
    if(nhit==0) hole++;
  }
return hole;
}

std::vector<MyCore_t> ShowerProcessor::CoreFinder(std::vector<Density_t> density)
{
  std::vector<MyCore_t> vec;
  std::vector<EVENT::CalorimeterHit*> _temp;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator hit=calohit.begin(); hit!=calohit.end(); ++hit){
    if(density[std::distance(calohit.begin(),hit)].dens3D<30) continue;
    //    std::cout << density[std::distance(hitmap.begin(),hit)] << std::endl;
    std::vector<EVENT::CalorimeterHit*> _hit;
    if(std::find(_temp.begin(),_temp.end(), *hit)!=_temp.end()) continue;
    _temp.push_back(*hit);
    _hit.push_back(*hit);
    for(std::vector<EVENT::CalorimeterHit*>::iterator hit2=calohit.begin(); hit2!=calohit.end(); ++hit2){
      if(density[std::distance(calohit.begin(),hit2)].dens3D<30) continue;
      if(std::find(_temp.begin(),_temp.end(), *hit2)!=_temp.end()) continue;
      for(std::vector<EVENT::CalorimeterHit*>::iterator it=_hit.begin(); it!=_hit.end(); ++it){
	if(abs(idDecoder(*it)["K-1"]-idDecoder(*hit2)["K-1"])<=1){ 
	  _hit.push_back(*hit2);
	  _temp.push_back(*hit2);
	  break;
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
  return vec;
}

int ShowerProcessor::NumberOfCores(std::vector<MyCore_t> vec)
{
  return vec.size();
}

int ShowerProcessor::CoreSize(MyCore_t core)
{
  return core.hits.size();
}

int ShowerProcessor::CoreBegin(MyCore_t core)
{
  int min=50;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    if(idDecoder(*it)["K-1"]<min) min=idDecoder(*it)["K-1"];
  }
  return min+1;
}

int ShowerProcessor::CoreEnd(MyCore_t core)
{
  int max=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    if(idDecoder(*it)["K-1"]>max) max=idDecoder(*it)["K-1"];
  }
  return max+1;
}

std::vector<float> ShowerProcessor::getCorePosition(MyCore_t core)
{
  float pos[3];
  for(size_t i=0; i<3; i++) pos[i]=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    pos[0]+=idDecoder(*it)["I"];
    pos[1]+=idDecoder(*it)["J"];
    pos[2]+=idDecoder(*it)["K-1"]+1;
  }
  std::vector<float> vec;
  for(size_t i=0; i<3; i++) vec.push_back(pos[i]/core.hits.size());
  return vec;
}

float ShowerProcessor::getCoreRadius(MyCore_t core)
{
  float radius=0;
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  for(std::vector<EVENT::CalorimeterHit*>::iterator it=core.hits.begin(); it!=core.hits.end(); ++it){
    radius+=(idDecoder(*it)["I"]-core.pos[0])*(idDecoder(*it)["I"]-core.pos[0])+
            (idDecoder(*it)["J"]-core.pos[1])*(idDecoder(*it)["J"]-core.pos[1]);
  }
  return sqrt(radius/CoreSize(core));
}

void ShowerProcessor::processRunHeader( LCRunHeader* run)
{
    _nRun++ ;
    _nEvt = 0;
} 



void ShowerProcessor::processEvent( LCEvent * evt )
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
      UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");    
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	calohit.push_back(hit);
      }
      makeHitMap();
      nlayer = Nlayer();
      begin = getShowerStartingLayer();
      zend = End();
      cog[0]=CenterOfXGravity();
      cog[1]=CenterOfYGravity();
      cog[2]=CenterOfZGravity(begin);
      radius = Radius(begin, cog);
      nlastplan=NhitLastPlan();
      hole=holeFinder(begin,cog);
      edge=Edge();
      if(//zend<100
	 begin>=0&&
	 //&& float(edge)/numElements<0.1 
	 //&& 
	 hole==0
	 /* && float(nlastplan)/numElements<15 
	    &&zend-begin>2 &&
	    && radius/cog[2]>0.20 &&
	 */){
	std::vector<int> nhitvec=Nhit();
	nhit=nhitvec[0];
	nhit1=nhitvec[1];
	nhit2=nhitvec[2];
	nhit3=nhitvec[3];
	longitudinal=HitinFirstPlates()/numElements;
	transversal=HitinCentralCells(cog)/numElements;
	hole1=hole1Finder(begin,cog);
        std::vector<Density_t> _density=Density();
	int size=_density.size();
	vecdensity.reserve(size);
	float ratio=0;
	for(int i=0; i<size; i++){
	  vecdensity.push_back(_density[i].dens3D);
	  if(_density[i].dens3D>30) ratio++;
	}
	EHratio=ratio/numElements;
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



void ShowerProcessor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ShowerProcessor::end(){ 
  file->Write();
  file->Close();
  std::cout << "ShowerProcessor::end()  " << name() 
     	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;  
}
