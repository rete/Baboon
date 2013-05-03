#define Analyze_cxx
#include "Analyze.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <TMath.h>
#include <TF1.h>
#include <vector>
#include <algorithm>

std::vector<Variable_t> Analyze::Loop(int E)
{
  std::vector<Variable_t> variable;
  if (fChain == 0) {
    std::cout << "BUG:The chain is empty" << std::endl;
    return variable;
  }
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    Variable_t var;
    if(E>=40&&Nhit<250) continue;
    if(End-Begin>5
       && End<49
       && Begin>=0
       //&& (Begin>5||Nlayer>30)
       && sqrt(TMath::Power(TransversalCut-0.94,2)+TMath::Power(LongitudinalCut-0.94,2))>0.2
       && float(Lasthit)/Nhit<0.15
       && Radius/CoG[2]>0.30
       && float(Nhit)/Nlayer>3
       && Hole==0 ){
      var.nhit.push_back(Nhit);
      var.nhit.push_back(Nhit1);
      var.nhit.push_back(Nhit2);
      var.nhit.push_back(Nhit3);
      var.NhitPerLay=float(Nhit)/Nlayer;
      var.radius=Radius;
      var.end=End;
      var.begin=Begin*2;
      var.cog=CoG[2]*2;
      var.ebeam=E;
      var.trackMultiplicity=TrackMultiplicity;
      for(std::vector<int>::iterator it=TrackLength->begin(); it!=TrackLength->end(); ++it){
	var.trackLength.push_back(*it);
      }
      if(int(var.trackLength.size())!=var.trackMultiplicity)
	std::cout << "Problem  with vector size about track" << std::endl;
      if(CoreSize->size()!=0){
	var.coresize=*std::max_element(CoreSize->begin(), CoreSize->end());
	var.corebegin=CoreBegin->at(std::distance(CoreSize->begin(),std::max_element(CoreSize->begin(),CoreSize->end())));
	var.coreend=CoreEnd->at(std::distance(CoreSize->begin(),std::max_element(CoreSize->begin(),CoreSize->end())));
      }
      var.density=*Density;
      variable.push_back(var);
    }
  }
  return variable;
}

void bookHisto()
{
  hist0=new TH1D("nhit","",100,0,1800);
  hist1=new TH1D("nhit1","",100,0,1500);
  hist2=new TH1D("nhit2","",100,0,600);
  hist3=new TH1D("nhit3","",50,0,150);
  hist4=new TH1D("radius","",60,0,30);
  hist5=new TH1D("begin","",50,0,100);
  hist6=new TH1D("cog","",50,0,100);
  hist7=new TH1D("trackMultiplicity","",20,0,20);
  hist8=new TH1D("Erec","",120,0,120);
  hist9=new TH1D("CoreSize","",100,0,1000);
  hist10=new TH1D("CoreLength","",30,0,30);
  hist11=new TH1D("TrackLength","",48,0,48);
}

void getHistoInfo(int E)
{
  fstream out;
  
  TF1 *func=new TF1("func","gaus",0,2000);
  hist0->Fit(func);
  hist0->Fit(func,"","",func->GetParameter(1)-2*func->GetParameter(2),func->GetParameter(1)+2*func->GetParameter(2));
  out.open("../TXTFile/nhit_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << func->GetParameter(1) << " " 
      << func->GetParError(1) << " " 
      << func->GetParameter(2) << " " 
      << func->GetParError(2) << std::endl;
  out.close();
  
  hist1->Fit(func,"","",hist1->GetMean()-2*hist1->GetRMS(),hist1->GetMean()+2*hist1->GetRMS());
  out.open("../TXTFile/nhit1_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << func->GetParameter(1) << " " 
      << func->GetParError(1) << " " 
      << func->GetParameter(2) << " " 
      << func->GetParError(2) << std::endl;
  out.close();
  
  hist2->Fit(func,"","",hist2->GetMean()-2*hist2->GetRMS(),hist2->GetMean()+2*hist2->GetRMS());
  out.open("../TXTFile/nhit2_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << func->GetParameter(1) << " " 
      << func->GetParError(1) << " " 
      << func->GetParameter(2) << " " 
      << func->GetParError(2) << std::endl;
  out.close();

  hist3->Fit(func);
  out.open("../TXTFile/nhit3_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << func->GetParameter(1) << " " 
      << func->GetParError(1) << " " 
      << func->GetParameter(2) << " " 
      << func->GetParError(2) << std::endl;
  out.close();

  hist4->Fit(func);
  out.open("../TXTFile/radius_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << func->GetParameter(1) << " " 
      << func->GetParError(1) << " " 
      << func->GetParameter(2) << " " 
      << func->GetParError(2) << std::endl;
  out.close();
  
  out.open("../TXTFile/begin_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << hist5->GetMean() << " " 
      << hist5->GetMeanError() << " " 
      << hist5->GetRMS() << " " 
      << hist5->GetRMSError() << std::endl;
  out.close();

  out.open("../TXTFile/cogz_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << hist6->GetMean() << " " 
      << hist6->GetMeanError() << " " 
      << hist6->GetRMS() << " " 
      << hist6->GetRMSError() << std::endl;
  out.close();

  out.open("../TXTFile/trackMulti_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << hist7->GetMean() << " " 
      << hist7->GetMeanError() << " " 
      << hist7->GetRMS() << " " 
      << hist7->GetRMSError() << std::endl;
  out.close();

  out.open("../TXTFile/TrackLength_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << hist11->GetMean() << " " 
      << hist11->GetMeanError() << " " 
      << hist11->GetRMS() << " " 
      << hist11->GetRMSError() << std::endl;
  out.close();  

  out.open("../TXTFile/coreSize_pi-_qgsp_bert.txt",ios::out | ios::app);
  out << E << " " 
      << hist9->GetMean() << " " 
      << hist9->GetMeanError() << " " 
      << hist9->GetRMS() << " " 
      << hist9->GetRMSError() << std::endl;
  out.close();
}

void writeHisto(int E,const char* output)
{
  TF1 *func=new TF1("func","gaus",0,120);
  hist8->Fit(func,"","",hist8->GetMean()-2*hist8->GetRMS(),hist8->GetMean()+2*hist8->GetRMS());
  //  hist8->Fit(func,"","",func->GetParameter(1)-2*func->GetParameter(2),func->GetParameter(1)+2*func->GetParameter(2));
  std::cout << func->GetParameter(2)/func->GetParameter(1) << std::endl;
  fstream out;
  out.open("../TXTFile/erec_pi-qgsp_bert.txt",ios::out|ios::app);
  out << E << " " 
      << func->GetParameter(1) << " " 
      << func->GetParError(1) << " "
      << func->GetParameter(2) << " " 
      << func->GetParError(2) << std::endl;
  TFile *myfile = new TFile(output,"recreate");
  hist0->Write();
  hist1->Write();
  hist2->Write();
  hist3->Write();
  hist4->Write();
  hist5->Write();
  hist6->Write();
  hist7->Write();
  hist8->Write();
  hist9->Write();
  hist10->Write();
  hist11->Write();
  myfile->Write();
  myfile->Close();
}

void deleteHisto()
{
  delete hist0;
  delete hist1;
  delete hist2;
  delete hist3;
  delete hist4;
  delete hist5;
  delete hist6;
  delete hist7;
  delete hist8;
  delete hist9;
  delete hist10;
  delete hist11;
}


void writeResult(std::vector<Variable_t> _var)
{
  fstream out;
  out.open("../TXTFile/simnhit_pi-_qgsp_bert.txt",ios::out | ios::app);
  int size=_var.size();
  for(int i=0; i<size; i++){
    if(_var[i].nhit[0]>hist0->GetMean()-2*hist0->GetRMS()&& 
       _var[i].nhit[0]<hist0->GetMean()+2*hist0->GetRMS())
      out << _var[i].ebeam << " "
	  << _var[i].nhit[0] << " "
	  << _var[i].nhit[1] << " "  
	  << _var[i].nhit[2] << " "  
	  << _var[i].nhit[3] << std::endl; 
  }
  out.close();

  out.open("../TXTFile/simGeomVariable_pi-_qgsp_bert.txt",ios::out | ios::app);
  for(int i=0; i<size; i++){
    out << _var[i].ebeam << " "
	<< _var[i].begin << " "
	<< _var[i].cog << " "
	<< _var[i].radius << std::endl; 
  }
  out.close();

  out.open("../TXTFile/simTrackMulti_pi-_qgsp_bert.txt",ios::out | ios::app);
  for(int i=0; i<size; i++){
    out << _var[i].ebeam << " "
	<< _var[i].trackMultiplicity << std::endl;
  }
  out.close();

  out.open("../TXTFile/simuTrackLength_pi-_qgsp_bert.txt",ios::out | ios::app);
  for(int i=0; i<size; i++){
    for(int j=0; j<_var[i].trackMultiplicity; j++)
      out << _var[i].ebeam << " "
	  << _var[i].trackLength[j] << std::endl;
  }
  out.close();
}


std::vector<float> Erec(std::vector<Variable_t> _var)
{
  std::vector<float> erec;
  double bestW[9];
//  bestW[0] = 0.0323015; bestW[1] = 2.07791e-05; bestW[2] = -2.7591e-08; 
//  bestW[3] = 0.108808; bestW[4] = 4.86495e-05; bestW[5] = 2.91007e-08; 
//  bestW[6] = 0.119161; bestW[7] = 0.000132741; bestW[8] = 3.4553e-08; 
//bestW[0] = 0.0349562; bestW[1] = 8.8134e-06; bestW[2] = -1.38467e-08; 
//bestW[3] = 0.108319; bestW[4] = 8.14387e-05; bestW[5] = -2.05943e-08; 
//bestW[6] = 0.0929921; bestW[7] = 0.000146943; bestW[8] = 3.3868e-08; 
  bestW[0] = 0.0310504; bestW[1] = 2.61963e-05; bestW[2] = -2.69564e-08; 
  bestW[3] = 0.117587; bestW[4] = 2.04374e-05; bestW[5] = 3.17688e-08; 
  bestW[6] = 0.114468; bestW[7] = 0.000111675; bestW[8] = 4.05933e-08;
//  bestW[0] = 0.0362364; bestW[1] = 5.84913e-05; bestW[2] = -2.50951e-08; 
//  bestW[3] = 0.087075; bestW[4] = -8.6531e-05; bestW[5] = 1.52749e-08; 
//  bestW[6] = 0.177564; bestW[7] = -7.97984e-05; bestW[8] = 5.05509e-08; 
//bestW[0] = 0.0305359; bestW[1] = 3.79054e-05; bestW[2] = -3.88033e-08; 
//bestW[3] = 0.0980147; bestW[4] = 8.61153e-06; bestW[5] = 6.60923e-08; 
//bestW[6] = 0.206969; bestW[7] = -6.70553e-05; bestW[8] = 1.25789e-07; 
//Calculation of bestW from migrad minimisation
  int size=_var.size();
  fstream out; out.open("../TXTFile/simerec_pi-_qgsp_bert.txt",ios::out|ios::app);
  for(int i=0; i<size; i++){
    erec.push_back((bestW[0]+_var[i].nhit[0]*bestW[1]+_var[i].nhit[0]*_var[i].nhit[0]*bestW[2])*_var[i].nhit[1]+
		   (bestW[3]+_var[i].nhit[0]*bestW[4]+_var[i].nhit[0]*_var[i].nhit[0]*bestW[5])*_var[i].nhit[2]+
		   (bestW[6]+_var[i].nhit[0]*bestW[7]+_var[i].nhit[0]*_var[i].nhit[0]*bestW[8])*_var[i].nhit[3]);
//    erec.push_back(_var[i].nhit[1]*bestX1+
//		   _var[i].nhit[2]*bestY1+
//		   _var[i].nhit[3]*bestZ1);
    out << _var[i].ebeam << " " << erec.back() << std::endl;
  } 
  return erec;
}

void Process(int E)
{
  bookHisto();

  
  for(int i=0; i<10; i++){
    char input[200];
    sprintf(input,"%s%d%s%d%s","./QGSP_BERT/single_pi-_",E,"GeV_I",i,".root");
    printf(input);
    std::cout << " " << std::endl;
    Analyze *a = new Analyze(input);
    std::vector<Variable_t> _var = a->Loop(E);
    int size=_var.size();
    std::vector<float> _erec=Erec(_var);
    for(int j=0; j<size; j++){
      hist0->Fill(_var[j].nhit[0]);
      hist1->Fill(_var[j].nhit[1]);
      hist2->Fill(_var[j].nhit[2]);
      hist3->Fill(_var[j].nhit[3]);
      hist4->Fill(_var[j].radius);
      hist5->Fill(_var[j].begin);
      hist6->Fill(_var[j].cog);
      hist7->Fill(_var[j].trackMultiplicity);
      hist8->Fill(_erec[j]);
      hist9->Fill(_var[j].coresize);
      hist10->Fill(_var[j].coreend-_var[j].corebegin);
      for(int k=0; k<_var[j].trackMultiplicity; k++){
	hist11->Fill(_var[j].trackLength[k]);
      }
    }
    writeResult(_var);
    delete a;
  }
  getHistoInfo(E);
  char output[200];
  sprintf(output,"%s%d%s","../ResultFile/QGSP_BERT/result_pi-_",E,".root");
  writeHisto(E,output);
  deleteHisto();
}

void ProcessAll()
{
  int Energy[]={5,10,15,20,25,30,40,50,60,70,80,90};
  int size=sizeof(Energy)/sizeof(int);
  fstream out;
  out.open("../TXTFile/simnhit_pi-_qgsp_bert.txt",ios::out);
  out << "" ;
  out.close();
  out.open("../TXTFile/simGeomVariable_pi-_qgsp_bert.txt",ios::out);
  out << "" ;
  out.close();
  for(int i=0; i<size; i++){
    Process(Energy[i]);
  }
}
