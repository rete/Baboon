#ifndef _CHAMBERANALYZER_H
#define _CHAMBERANALYZER_H
#include <limits.h>
#include "DHCALAnalyzer.h"
#include "DHCalEventReader.h"
#include "DCHistogramHandler.h"
#include <iostream>
#include <sys/timeb.h>
#include "sqlite3.h"
#include <mysql/mysql.h>



#include "IMPL/LCTOOLS.h"
#include "EVENT/RawCalorimeterHit.h" 
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/RawCalorimeterHitImpl.h"
#include "UTIL/CellIDDecoder.h"
#include "DifGeom.h"
#include "ChamberGeom.h"

#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"

#include "StructTree.h"
track_t theTrack;
point_t thePoint;
event_t theEvent;
cluster_t theCluster;
hit_t theHit;
shower_t theShower;
static double pad2cm=1.04125;
//static double pad2cm=1.0;

template <class T>
class array3D
{
	T* ptr_;
	uint32_t xs_,ys_,zs_;
public:
	array3D(){ptr_=0;}
	inline void initialise(T* p,uint32_t x_size1=60,uint32_t y_size1=96,uint32_t z_size1=96)
	{
		ptr_=p;
		xs_=x_size1;
		ys_=y_size1;
		zs_=z_size1;
	}

	inline void Fill(T* p){ memcpy(ptr_,p,xs_*ys_*zs_*sizeof(T));}
	inline void clear(){memset(ptr_,0,xs_*ys_*zs_*sizeof(T));}
	inline T* getPtr(){return ptr_;}
	inline T  getValue(uint32_t i,uint32_t j,uint32_t k){return ptr_[(i*ys_+j)*zs_+k];}
	inline void setValue(uint32_t i,uint32_t j,uint32_t k,T val){ptr_[(i*ys_+j)*zs_+k]=val;}

	inline uint32_t getXSize(){return xs_;}
	inline uint32_t getYSize(){return ys_;}
	inline uint32_t getZSize(){return zs_;}
};

class HTImage
{
public:
	HTImage(uint32_t nbinx,float xmin,float xmax,uint32_t nbiny,float ymin,float ymax);
	~HTImage();
	void Clear();
	void addPixel(float x,float y,float w=1.,uint32_t ch=INT_MAX,uint32_t pad=INT_MAX);
	void findMaximum(uint32_t& maxval,float& theta,float& r);
	void findMaxima(std::vector<uint32_t>& maxval,std::vector<float>& theta,std::vector<float>& r);
	void Draw(DCHistogramHandler* h);
private:
	uint16_t* theImage_;
	uint16_t* theOriginalImage_;
	int32_t theNbinx_;
	float theXmin_,theXmax_;
	int32_t theNbiny_;
	float theYmin_,theYmax_,theBinxSize_,theBinySize_;
};



class RecoHit;



class Shower
{
public:
	Shower(RecoHit&h);
	void clear();
	bool append(RecoHit& h,float dist_cut);
	void Add(RecoHit& h);
	double Distance(RecoHit& h);
	uint32_t getNumberOfHits(uint32_t plan,uint32_t threshold);
	double getCorrectedNumberOfHits(uint32_t plan,uint32_t threshold,std::map<uint32_t,double*> corr);   
	uint32_t getReduceNumberOfHits(uint32_t threshold,uint32_t fp=0,uint32_t lp=100);
	uint32_t getNumberOfHits(uint32_t threshold);
	uint32_t getNumberOfMips(uint32_t plan);
	uint32_t getFDHits(uint32_t* v,uint32_t thr);
	uint32_t getFDHitsN(uint32_t* v,uint32_t thr);
	std::map<uint32_t,std::vector<RecoHit> >& getPlans(){return thePlans_;}
	void PlayMatrix(uint32_t fp=1,uint32_t lp=60);
	void transverseProfile(uint32_t plan,uint32_t &nh,double &xb,double &yb, double &l0, double &l1,double* v0,double *v1,double &n9,double &n25);
	static void computePrincipalComponents(std::vector<RecoHit*> &v, double result[21]);
	static void culaPrincipalComponents(std::vector<RecoHit*> &v, double result[21]);

	double getl1(){return l1_;}
	double getl2(){return l2_;}
	double getl3(){return l3_;}
	double* getv1(){return v1_;}
	double* getv2(){return v2_;}
	double* getv3(){return v3_;}
	double* getxm(){return xm_;}
	void setSelected(bool t){selected_=t;}
	bool isSelected(){return selected_;}
	double closestDistance(Shower& sh);
	uint32_t getFirstPlan(){return firstPlan_;}
	uint32_t getLastPlan(){return lastPlan_;}
private:
	std::map<uint32_t,std::vector<RecoHit> > thePlans_;
	double l1_,l2_,l3_;
	double v1_[3],v2_[3],v3_[3];
	double xm_[3];
	bool selected_;
	uint32_t firstPlan_,lastPlan_;

};

class RecoHit
{
public:
	RecoHit(){;}
	RecoHit(DifGeom& d, ChamberGeom& c,IMPL::RawCalorimeterHitImpl* h,uint32_t hrtype=2);
	void initialise(DifGeom& d, ChamberGeom& c,IMPL::RawCalorimeterHitImpl* h,uint32_t hrtype=2);
	double X(){return x_;}
	double Y(){return y_;}
	double Z(){ return cg_.getZ();}
	uint8_t I(){return chamberLocalI_;}
	uint8_t  J(){return chamberLocalJ_;}
	unsigned short chamber(){return cg_.getId();}
	unsigned short dif(){return dg_.getId();}
	uint8_t getDifI(){ return difLocalI_;}
	uint8_t getDifJ(){ return difLocalJ_;}
	uint8_t getAsic(){ return 0xFF & (raw_->getCellID0()&0xFF00)>>8;}
	uint8_t getChannel(){ return 0xFF & (raw_->getCellID0()&0x3F0000)>>16;}
	uint16_t getAmplitude(){return raw_->getAmplitude();}
	//  inline void setIndices(uint32_t i,uint32_t n) { index_=i;next_=n;}
	//inline uint32_t getIndex(){return index_;}
	//inline uint32_t getNext(){return next_;}
	Shower* getShower(){return shower_;}
	void setShower(Shower* s){shower_=s;}
private:
	DifGeom dg_;
	ChamberGeom cg_;
	IMPL::RawCalorimeterHitImpl* raw_;
	uint8_t difLocalI_,difLocalJ_;
	uint8_t chamberLocalI_,chamberLocalJ_;
	// uint32_t index_,next_;
	double x_,y_;
	Shower* shower_;
};

/**
\class RECOCluster
\author L.Mirabito
\date May 2010
\version 1.0
\brief Vector of RecoHit. The position is the mean of X and Y. No usage of threshold
*/
class RECOCluster
{
public:
	RECOCluster(){valid_=true;}
	RECOCluster(RecoHit h);
	~RECOCluster();
	double dist(RecoHit h1,RecoHit h2);
	bool Append(RecoHit h);
	std::vector<RecoHit>* getHits();
	bool isAdjacent(RECOCluster &c);
	void setValidity(bool t){valid_=t;}
	bool isValid() {return valid_;}
	double Pos(int p);
	void Print();
	double X();
	double Y();
	double dX();
	double dY();
private:
	void calcPos();
	std::vector<RecoHit> hits_;
	double x_,y_,dx_,dy_;
	bool valid_;
};


/**
\class RecoPoint
\author L.Mirabito
\date May 2010
\version 1.0
\brief Vector of RecoHit. The position is the mean of X and Y. No usage of threshold
*/
class RecoPoint
{
public:
	RecoPoint(RECOCluster& h,unsigned int ch,double x , double y, double z,double dx=0.5,double dy=0.5);
	void Print();
	double X(){ return x_;}
	double Y(){return y_;}
	double dX(){ return dx_;}
	double dY(){return dy_;}
	double Z(){return z_;}
	double Charge(){return weight_;}
	void calcPos();
	unsigned int getChamberId(){ return chId_;}
	RECOCluster& getCluster(){return h_;}
	void setUsed(bool t){inTrack_=t;}
	bool isUsed(){return inTrack_;}
	void setPointId(uint32_t i){ptid_=i;}
	uint32_t getPointId(){return ptid_;}
private:
	unsigned int chId_;
	RECOCluster h_;
	double x_,y_,z_,x2_,y2_;
	double dx_,dy_;
	bool inTrack_;
	uint32_t ptid_;
	double weight_,weight2_;
};


class RecoCandTk
{
public:
	RecoCandTk();
	~RecoCandTk();
	void clear();
	void addPoint(RecoPoint& p);
	bool addNearestPoint(RecoPoint& p);
	bool addPoint(RecoPoint& p,double xcut, double ycut);
	bool addPoints(std::vector<RecoPoint> v, double zref,double xcut, double ycut);
	bool addPoints(std::vector<RecoPoint> v, double dcut);
	bool addChi2Points(std::vector<RecoPoint> v, double dcut,std::vector<std::vector<RecoPoint>::iterator>* used=NULL);
	void removeDistantPoint(float zcut);

	void Refit(RecoCandTk &t,float c);
	void regression();
	void regression1D(std::vector<double> vx,std::vector<double> weight,std::vector<double> y,double &chi2, double &alpha,double &beta);
	void calculateChi2();
	void clean();
	void Print();
	std::vector<RecoPoint*>& getList(){return list_;}
	std::vector<double>& getChi2(){return dist_;}
	double getXext(double z) { return ax_*z+bx_;}
	double getYext(double z) { return ay_*z+by_;}  
	double calculateDistance(RecoPoint& p);
	double ax_,bx_,ay_,by_,chi2_;
	double prChi2_;
	uint32_t firstChamber_,lastChamber_;
	inline bool isValid(){return valid_;}
	inline void setValid(bool t){valid_=t;}
	float zmin_,zmax_;
private:
	bool valid_;
	std::vector<RecoPoint*> list_;
	std::vector<double> dist_;
	float dmin_[61];

};

class Amas
{

public:
	Amas(RecoHit* h)
	{
		theHits_.push_back(h);
	}
	Amas(std::vector<RecoHit*>& vh)
	{
		theHits_.clear();
		this->copyFrom(vh);
		this->compute();
		//theHits_.push_back(h);
	}
	void add(RecoHit* h) {theHits_.push_back(h);}
	bool append(RecoHit* h,uint32_t del=1)
	{
		bool appended=false;
		for (std::vector<RecoHit*>::iterator it=theHits_.begin();it!=theHits_.end();it++)
		{
			if (abs((*it)->chamber()-h->chamber())>del) continue;
			if (abs((*it)->I()-h->I())>del) continue;
			if (abs((*it)->J()-h->J())>del) continue;
			appended=true;
		}
		if (appended) theHits_.push_back(h);
		return appended;
	}
	void compute()
	{
		Shower::computePrincipalComponents(theHits_,theComponents_);
		memset(nh,0,3*sizeof(uint32_t));
		for (std::vector<RecoHit*>::iterator it=theHits_.begin();it!=theHits_.end();it++)
		{
			int ithr= (*it)->getAmplitude()&0x3;
			if (ithr==1) nh[1]++;
			if (ithr==2) nh[0]++;
			if (ithr==3) nh[2]++;
		}
	}
	double getComponents(uint32_t i) {return theComponents_[i];}
	double* Components() {return &theComponents_[0];}
	uint32_t* Hits(){return &nh[0];}
	double X(){return theComponents_[0];}
	uint32_t size(){return theHits_.size();}
	std::vector<RecoHit*> &getHits(){return theHits_;}


	bool operator> (const Amas& other) const
	{
		return this->theHits_.size() > other.theHits_.size();
	}
	bool operator>=(const Amas& other) const
	{
		return this->theHits_.size() >= other.theHits_.size();
	}
	bool operator< (const Amas& other) const
	{
		return this->theHits_.size() < other.theHits_.size();
	}
	bool operator<=(const Amas& other) const
	{
		return this->theHits_.size() <= other.theHits_.size();
	}

	void copyTo(std::vector<RecoHit*>& vh)
	{
		for (std::vector<RecoHit*>::iterator it=theHits_.begin();it!=theHits_.end();it++)
		{
			vh.push_back((*it));
		}
	}
	void copyFrom(std::vector<RecoHit*>& vh)
	{
		for (std::vector<RecoHit*>::iterator it=vh.begin();it!=vh.end();it++)
		{
			theHits_.push_back((*it));
		}
	}
private:
	std::vector<RecoHit*> theHits_;
	double theComponents_[21];
	uint32_t nh[3];
};
class HC
{
public:
	HC(uint32_t m,float th,float r) : m_(m),th_(th),r_(r)
	{
		a_= -1./tan(th_);
		b_= r_/sin(th_)-50.*a_;
		points_.clear();
		for (int i=0;i<61;i++) dmin_[i]=1E9;
	}
	inline float Pos(float z) {return a_*z+b_;}
	inline void setDmin(float d,uint32_t ch){if (d<dmin_[ch]) dmin_[ch]=d;}
	inline void add(std::vector<RecoPoint>::iterator i){points_.push_back(i);}
	inline void Dump(){printf("=====> %d %f %f %d  \n",m_,a_,b_,points_.size());}
	inline uint32_t common(HC& other,std::vector<std::vector<RecoPoint>::iterator > &v)
	{
		uint32_t nc=0;
		for (std::vector< std::vector<RecoPoint>::iterator >::iterator ip=other.points_.begin();ip!=other.points_.end();ip++)
		if (std::find(points_.begin(),points_.end(),(*ip))!=points_.end()) 
		{nc++; v.push_back(*ip);}
		return nc;
	}
	uint32_t m_;
	float th_,r_,a_,b_,dmin_[61];

	std::vector<std::vector<RecoPoint>::iterator> points_;
};


class ChamberAnalyzer : public DHCALAnalyzer
{
public:
	ChamberAnalyzer(DHCalEventReader* r,DCHistogramHandler* h);
	virtual ~ChamberAnalyzer(){;}
	virtual void processEvent();
	virtual void initHistograms();
	virtual void processRunHeader()
	{
		if (writing_)
		reader_->writeRunHeader();
	}
	void presetParameters();
	void setWriting(bool t){writing_=t;}
	virtual void initJob();
	virtual void endJob();
	virtual void initRun(){;}
	virtual void endRun(){;}

	bool decodeTrigger(LCCollection* rhcol, double tcut);

	void trackHistos();
	void findTracks();
	void findTracks1();
	void CosmicFinder();
	void appendHits(RecoCandTk& t);
	void setCollectionName(std::string s){ collectionName_=s;}
	unsigned long long getExternalTriggerTime() { return (unsigned long long) long(externalTriggerTime_);}
	unsigned long long getAbsoluteFrameTime(int bc) { return (unsigned long long) long(externalTriggerTime_-bc);}


	void setDropFirstSpillEvent(bool t){dropFirstSpillEvent_=t;}

	void setClockSynchCut(unsigned int t){clockSynchCut_=t;}
	void setSpillSize(double t){spillSize_=t;}

	uint32_t NoiseStudy(std::map<uint32_t,std::bitset<255> > timeDif,std::map<uint32_t,std::bitset<61> > timeChamber);
	void FillTimeAsic(IMPL::LCCollectionVec* rhcol);
	void DIFStudy(IMPL::RawCalorimeterHitImpl* hit);

	void ShowerBuilder(std::vector<RecoHit*> vreco);
	void ImageBuilder(std::vector<RecoHit*> &vreco);
	void sobel_filtering(unsigned char i[60][96], float j[60][96]);

	void sobel_volume(unsigned char image1[60][96][96],float image2[60][96][96] );

	void drawDisplay();
	void HT();
	void HT2D();
	void HT3D();
	void HTOld();
	double checkTime();
	void findHoughCandidates3D(uint32_t stop,std::vector<RecoCandTk> &tkSeed);

	void findHoughCandidates(uint32_t type,uint32_t stop,std::vector<HC> &vX);
	void mergeHoughCandidate(std::vector<HC> &vX,std::vector<HC> &vY,std::vector<RecoCandTk> &tkSeed);
	inline void setuseSynchronised(bool t){useSynchronised_=t; }
	inline void setoldAlgo(bool t){oldAlgo_=t;}
	inline bool getoldAlgo(){return oldAlgo_;}
	inline void settkMinPoint(int t){tkMinPoint_=t;}
	inline void settkExtMinPoint(int t){tkExtMinPoint_=t;}
	inline void settkBigClusterSize(int t){tkBigClusterSize_=t;}
	inline void setspillSize(float t){spillSize_=t;}
	inline void settkChi2Cut(float t){tkChi2Cut_=t;}
	inline void settkDistCut(float t){tkDistCut_=t;}
	inline void settkExtChi2Cut(float t){tkExtChi2Cut_=t;}
	inline void settkExtDistCut(float t){tkExtDistCut_=t;}
	inline void settkAngularCut(float t){tkAngularCut_=t;}
	inline void setclockSynchCut(int t){clockSynchCut_=t;}
	inline void setminChambersInTime(int t){minChambersInTime_=t;}
	inline void setmaxHitCount(int t){maxHitCount_=t;}
	inline void setchamberEdge(float t){chamberEdge_=t;}
	inline void setrebuild(bool t) {rebuild_=t;}
	inline void setcollectionName(std::string t){collectionName_=t;}


	void createTrees(std::string s);
	void closeTrees();

	void PointsBuilder(std::vector<RecoHit*> &vrh);

	void TracksBuilder(Shower &ish,std::vector<RecoHit*> &vrh);
	void FillTrackTree();
	void readCalibration(uint32_t run) throw (std::string); 

	void openSqlite(std::string fileName);
	int32_t fillEventTable(uint32_t run, uint32_t event,uint64_t bcid,uint32_t time,uint32_t e3,uint32_t c3,int32_t w,uint32_t np);
	int32_t fillComponentsTable(uint32_t *nh,double *components);
	int32_t fillAmasTable(Amas* a);
	int32_t fillAmasSummaryTable(uint32_t core[3],uint32_t edge[3],uint32_t namas,uint32_t nafter,uint32_t tag,double zlast);
	int32_t fillShowerTable(shower_t &s);
	int32_t executeQuery(std::string stmt);
	uint32_t getLastInsertId();
	void decodeAccount(std::string account);
	void connect(std::string account);

	void findTimeSeeds( IMPL::LCCollectionVec* rhcol, int32_t nhit_min,std::vector<uint32_t>& candidate);
	void buildVolume(IMPL::LCCollectionVec* rhcol,uint32_t seed);
	//void EdgeDetection(unsigned char imagev[60][96][96],unsigned char core[60][96][96],unsigned char edge[60][96][96]);
	void EdgeDetection(array3D<unsigned char> &i,array3D<unsigned char> &c,array3D<unsigned char> &e);
	//void sobel_volume(unsigned char *image1,float *image2,uint32_t x_size1,uint32_t y_size1,uint32_t z_size1 );
	void sobel_volume(array3D<unsigned char> &im1,array3D<float> &im2);
	void draw(array3D<unsigned char> &all,array3D<unsigned char> &core,array3D<unsigned char> &edge);

	uint32_t mergeAmas(array3D<unsigned char> &core,array3D<unsigned char> &edge);

	void newHT(array3D<unsigned char> &core);

private:


	int nAnalyzed_;
	int nInSynch_,run_;
	double externalTriggerTime_,lastSpill_,lastPowerPulsedTime_;


	std::vector<RecoPoint> allpoints_;
	std::vector<RecoCandTk> tklist_;
	std::vector<RecoCandTk> tkgood_;
	std::map<uint32_t,std::vector<RecoPoint*> > chamberPoints_;
	std::map<uint32_t,uint32_t> asicCount_;
	std::map<uint32_t,double*> theCorreff_;
	double integratedTime_;
	RecoHit hitVolume_[60][96][96];

	DCHistogramHandler* rootHandler_;
	int hrtype_;
	int lastrunnb_;
	unsigned int currentTime_,currentEvent_,lastSpyEvent_;
	// Control
	bool rebuild_;
	bool dropFirstSpillEvent_;
	bool findTracks_;
	bool useSynchronised_;
	bool oldAlgo_;
	int tkMinPoint_;
	int tkExtMinPoint_;
	int tkBigClusterSize_;
	float spillSize_;
	float tkChi2Cut_;
	float tkDistCut_;
	float tkExtChi2Cut_;
	float tkExtDistCut_;
	float tkAngularCut_;
	int clockSynchCut_;
	int minChambersInTime_;
	int maxHitCount_;
	int tkFirstChamber_;
	int tkLastChamber_;
	float chamberEdge_;
	uint32_t trackIndex_;
	uint32_t houghIndex_;
	bool useTk4_;
	std::string collectionName_;
	uint32_t offTimePrescale_;


	// Reader
	DHCalEventReader* reader_;

	bool writing_;
	bool headerWritten_;
	bool draw_;

	IMPL::LCEventImpl* evt_;

	HTImage* theHTx_;
	HTImage* theHTy_;

	struct timeb theTime_;
	struct timeb theCurrentTime_;

	double theRhcolTime_;
	double theTimeSortTime_;
	double theTrackingTime_;
	double theHistoTime_;
	double zLastAmas_;
	int theSeuil_;
	int32_t theSkip_,npi_;
	unsigned long long theBCID_;
	uint32_t theDTC_,theGTC_;
	TTree* tEvents_;
	TTree* tShowers_;
	TTree* tTracks_;
	TFile* treeFile_;
	uint32_t theNall_,theNedge_,theTag_,theEdge3_,theCore3_,theBigCore_,theCore1_,theCore2_,theNplans_;
	int32_t theWeights_;
	float theZMax_,theZFirst_,theZLast_,thePlanAfter_;
	TNtuple* theNtuple_;
	TFile* theNtupleFile_;
	bool useSqlite_,useMysql_;
	sqlite3* theDb_;
	uint64_t theEventRowId_;
	MYSQL theMysql_;
	std::string myName_,myPwd_,myHost_,myDatabase_;
	array3D<unsigned char> theImage_;
	array3D<unsigned char> theImageWeight_;
	unsigned char theImageBuffer_[60*96*96];
	unsigned char theImageWeightBuffer_[60*96*96];
	unsigned char theImageCoreBuffer_[60*96*96];
	unsigned char theImageEdgeBuffer_[60*96*96];

	std::vector<Amas> theAmas_;
};
#endif
