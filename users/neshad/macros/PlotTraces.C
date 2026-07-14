#include <iostream>
#include <sstream>
#include <fstream>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TMultiGraph.h>


#include "LendaEvent.hh"
#include "LendaFilter.hh"

void PlotTraces(){
	bool wait = true;
  // user input //
  //Using entry as gamma, entry2 as neutron
  Int_t entry = 46000;   //46000
  TString runnum = "5248";
  TString barname = "RMD1";
  Int_t entry2 = 25000;   //25000
  
  
  // end user input // 

  // Get tree
  //TFile* f = new TFile("./rootfiles/run-"+runnum+"-00.root");
    TFile* f = new TFile("/mnt/daqtesting/lenda/rootfiles/run-"+runnum+"-00.root");

  TTree* tree = (TTree*)f->Get("caltree");
  LendaEvent *lendaevent = new LendaEvent();
  tree->SetBranchAddress("lendaevent",&lendaevent);

  // Get entry 1

  tree->GetEntry(entry);

  // Various checks

  Int_t numBarsInEvent = 1;
    
  if (lendaevent->NumBars != numBarsInEvent) {
    cout << "Wrong number of bars in event: " << numBarsInEvent << " bars." << endl;
    return;
  }
	cout<<lendaevent->Bars[0].Name<<endl;
  Int_t barID = -1;
  for(Int_t i=0; i<numBarsInEvent; i++) {
    if(lendaevent->Bars[i].Name==barname) { barID=i; continue; }
  }
  if(barID<0){
    cout << "Can't find bar " << barname << endl;
    return;
  }
  

  // Get the traces for entry 1
  vector <UShort_t> tracetop = lendaevent->Bars[barID].Tops[0].GetTrace();
  double ph_gamma = lendaevent->Bars[barID].Tops[0].GetPulseHeight();
	int tracesize = tracetop.size();

	cout << "Trace is of length: "<< tracesize <<endl;
	
	cout<< "gamma PH: "<<ph_gamma<<endl;
	
	//Get entry 2
	
	
	tree->GetEntry(entry2);

  // Various checks

  
    
  if (lendaevent->NumBars != numBarsInEvent) {
    cout << "Wrong number of bars in event: " << numBarsInEvent << " bars." << endl;
    return;
  }

  
  for(Int_t i=0; i<numBarsInEvent; i++) {
    if(lendaevent->Bars[i].Name==barname) { barID=i; continue; }
  }
  if(barID<0){
    cout << "Can't find bar " << barname << endl;
    return;
  }
  
cout<<"got here" <<endl;
  // Get the traces
  vector <UShort_t> trace2 = lendaevent->Bars[barID].Tops[0].GetTrace();
  double ph_neutron = lendaevent->Bars[barID].Tops[0].GetPulseHeight();
  
  cout <<"got traces"<< endl;
  
  LendaFilter aFilter;
  Int_t start = 1;
  Double_t rise = -1.;
  Double_t background = -1.;
  Double_t noise = -1.;
  Int_t ph_neutron_max = aFilter.GetMaxPulseHeight(trace2,start);
  Int_t ph_neutron_net = aFilter.GetPulseComplete(trace2,start,rise,background,noise);
	
	
  cout<< "neutron PH: "<<ph_neutron<<endl;
  cout<< "max PH: " << ph_neutron_max <<endl;
  cout<< "Net neutron PH: "<<ph_neutron_net <<endl;
	int trace2size = trace2.size();
	
	
  
  // Fill the graphs
  Double_t time[tracesize];
  Double_t gamma[tracesize];
  Double_t neutron[trace2size];
  Double_t maxgam = 0, mingam = 999999;
  Double_t maxneut = 0, minneut = 999999;
  for (int i=0; i < tracesize; i++){
    time[i] = i;
    gamma[i] = tracetop[i];
	neutron[i] = trace2[i];
    if(gamma[i]>maxgam) maxgam=gamma[i];
    if(gamma[i]<mingam) mingam=gamma[i];
	if(neutron[i]>maxneut) maxneut=gamma[i];
    if(neutron[i]<minneut) minneut=gamma[i];

  }

  
  
  
  TCanvas* c = new TCanvas("c");
  TMultiGraph *mg = new TMultiGraph();
  TGraph * GammaTraceGraph = new TGraph(tracesize,time,gamma);
  GammaTraceGraph->SetTitle(barname+" Traces ("+runnum+")");
  GammaTraceGraph->SetMarkerColorAlpha(kBlue,1);
  GammaTraceGraph->SetLineColorAlpha(kBlue,1);
  GammaTraceGraph->GetHistogram()->GetXaxis()->SetRangeUser(0,130);
  GammaTraceGraph->Draw("AL");

	
  TGraph *NeutronTraceGraph = new TGraph(trace2size,time,neutron);
  NeutronTraceGraph->SetMarkerColorAlpha(kRed,1);
  NeutronTraceGraph->SetLineColorAlpha(kRed,1);
  NeutronTraceGraph->Draw("L");

  
  TLegend* leg =new TLegend(0.5, 0.5, 0.9, 0.9);
  leg->SetHeader(barname+" Traces ("+runnum+")");
  leg->AddEntry(GammaTraceGraph, "Gamma trace", "l");
  leg->AddEntry(NeutronTraceGraph, "Neutron trace", "l");
  leg->Draw("P");
  
  mg->Add(GammaTraceGraph);
  mg->Add(NeutronTraceGraph);
  mg->Draw("L*");


  if(wait) { gPad->WaitPrimitive(); }
  
  //c->SaveAs(barname+" Traces_("+runnum+")" + "_g_"+ Form("%d",entry) + "_n_" + Form("%d",entry2) + ".jpg");

  delete c;

 f->Close();

}