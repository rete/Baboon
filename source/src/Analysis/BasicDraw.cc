  /// \file BasicDraw.cc
/*
 *
 * BasicDraw.cc source template generated by fclass
 * Creation date : ven. mai 3 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Analysis/BasicDraw.hh"

using namespace std;

namespace baboon {


	void CaliceStyle() {

	  gROOT->SetStyle("Plain");
	  gStyle->SetPalette(1);
	  gStyle->SetPadTickX(1);
	  gStyle->SetPadTickY(1);
	  gStyle->SetLabelFont(42,"xyz");
	  gStyle->SetTitleFont(42);
	  gStyle->SetTitleFont(42,"xyz");
	  gStyle->SetStatFont(42);

	  gStyle->SetFrameFillColor(kWhite);
	  gStyle->SetCanvasColor(kWhite);
	  gStyle->SetOptStat(0); /*don't show statistics box*/
	  gStyle->SetTitleSize(0.05, "xyz");
	  gStyle->SetLegendBorderSize(1);

	  gStyle->SetPadTopMargin(0.05);
	  gStyle->SetPadBottomMargin(0.14);
	  gStyle->SetPadLeftMargin(0.15);
	  gStyle->SetPadRightMargin(0.05);

	  gROOT->ForceStyle();
	}


	TPaveText *DrawText( TString string ) {

	  TPaveText *pt = new TPaveText(0.10, 0.91, 0.43, 0.96, "tbNDC");
	  pt->SetTextSize(0.04);
	  pt->SetTextColor(kGray);
	  pt->SetFillColor(0);
	  pt->SetLineWidth(0);
	  pt->SetBorderSize(0);
	  pt->SetTextAlign(22); //center
	  pt->AddText(string);
	  pt->Draw();
	  pt->SetTextFont(62);
	  return pt;
	}


	TGraph *NewTGraph( unsigned int nbOfPts , unsigned int colorLineID ) {

		TGraph *graph = new TGraph(nbOfPts);
		graph->SetLineColor(colorLineID);
		graph->SetLineWidth(2);
		graph->SetMarkerStyle(21);
		graph->SetMarkerColor(1);
		graph->SetMarkerSize(0.5);
		return graph;
	}


	TGraph *NewTGraph( unsigned int nbOfPts , double *x , double *y , unsigned int colorLineID ) {

		TGraph *graph = new TGraph(nbOfPts,x,y);
		graph->SetLineColor(colorLineID);
		graph->SetLineWidth(2);
		graph->SetMarkerStyle(21);
		graph->SetMarkerColor(1);
		graph->SetMarkerSize(0.5);
		return graph;
	}


	TLegend *NewTLegend( double x1 , double y1 , double x2 , double y2 ) {

		TLegend *legend = new TLegend(x1,y1,x2,y2);
		legend->SetMargin(0.3);
		legend->SetLineColor(0);
		legend->SetLineWidth(0);
		legend->SetFillColor(0);
		legend->SetFillStyle(0);
		return legend;
	}

	TH2D *NewTH2D( const string &histName
				   ,const string &histTitle
				   ,unsigned int nbBinX
				   ,int xMin
				   ,int xMax
				   ,unsigned int nbBinY
				   ,int yMin
				   ,int yMax
				   ,int color ) {

		TH2D *h2 = new TH2D( histName.c_str() , histTitle.c_str() , nbBinX , xMin , xMax ,  nbBinY , yMin , yMax );
		h2->SetMarkerStyle(20);
		h2->SetMarkerSize(.4);
		h2->SetMarkerColor( color );
		return h2;
	}

	TH1D *NewTH1D( const string &histName
				   ,const string &histTitle
				   ,unsigned int nbBinX
				   ,int xMin
				   ,int xMax
				   ,int color ) {

		TH1D *h1 = new TH1D( histName.c_str() , histTitle.c_str() , nbBinX , xMin , xMax );
//		h1->SetLineStyle(0);
		h1->SetLineWidth(2);
		h1->SetLineColor( color );
		return h1;
	}


}  // namespace 

