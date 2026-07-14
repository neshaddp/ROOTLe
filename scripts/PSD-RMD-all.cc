#include <stdio.h>
#include <stdlib.h>
//test
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <TRandom1.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "TRandom3.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TChain.h"
#include "TRandom2.h"

//Local Headers
#include "LendaEvent.hh"
#include "ddaschannel.hh"
#include "DDASEvent.hh"


#include <ctime>
#include <iomanip>

#include "LendaPacker.hh"

#include "TFitResult.h"
#include "TF1.h"
#define BAD_NUM -10008

#include "R00TLeSettings.hh"

using namespace std;


int main(int argc, char **argv){

  R00TLeSettings aSettings;
  //Make a packer object
  LendaPacker *thePacker = new LendaPacker(&aSettings);
  

  

  //prepare files 
  ////////////////////////////////////////////////////////////////////////////////////
  TChain *inT=new TChain("caltree");
  //inT->Add("./../users/jorge/rootfiles/run-4207-00-rmd.root");
  inT->Add("./../users/jorge/rootfiles/run-4601-00-rmd.root");
  //inT->Add("./../users/jorge/rootfiles/run-4607-00-rmd.root");


  Long64_t nentry=(Long64_t) (inT->GetEntries());
  cout <<"The number of entries is : "<< nentry << endl ;
  



  // Openning output Tree and output file
  ////////////2////////////////////////////////////////////////////////////////////////
  TFile *outFile = new TFile("PSD-2017test.root","recreate");
 

 
  int shortgate_low = 15;
  int shortgate_high = 20;
  int longgate_low = 180;
  int longgate_high = 200;
 
  //int shortgate_low = 50;
  //int shortgate_high = 130;
  //int longgate_low = 500;
  //int longgate_high = 900;

  int NumberOfParameterSets = ((shortgate_high -shortgate_low)/5.)*((longgate_high - longgate_low)/20.);
  int NumberOfHistogramSets = NumberOfParameterSets * 9; // Nine RMD detectors
  int NumberRMD = 9; //Number of RMD detectors
  //int NumberOfHistogramSets = NumberOfParameterSets * 1; // One RMD detectors
  //int NumberRMD = 1; //Number of RMD detectors

  cout<<"YOU ASKED FOR "<<NumberOfParameterSets <<" parameter sets...which gives "<<NumberOfHistogramSets<< " histogram sets"<<endl;
 


  //vector <TH2F*> PSD_vs_E(NumberOfParameterSets);
  //vector <TH1F*> PSD1(NumberOfParameterSets);
  //vector <TH1F*> PSD2(NumberOfParameterSets);
  //vector <TH1F*> PSD3(NumberOfParameterSets);
 
  vector < vector <TH2F*> > RMD_PSD_vs_E(NumberRMD);
  vector < vector <TH1F*> > RMD_PSD1(NumberRMD);
  vector < vector <TH1F*> > RMD_PSD2(NumberRMD);
  vector < vector <TH1F*> > RMD_PSD3(NumberRMD);


  stringstream nameStream;
  float PSDlow=0.;
  float PSDhigh=1.0;
  int PSDBins=500;
  int Elow=0;
  int Ehigh=2000;
  int EBins=1000;
  
  for (int rmd_index=0; rmd_index<NumberRMD; rmd_index++) { // Loop over all RMDs
    int RMD = rmd_index + 1; 
    //int RMD = rmd_index + 4; // We are only anlyzing RMD104 

    //cout << "test 1" << endl;
    RMD_PSD1[rmd_index].resize(NumberOfParameterSets);
    RMD_PSD2[rmd_index].resize(NumberOfParameterSets);
    RMD_PSD3[rmd_index].resize(NumberOfParameterSets);
    RMD_PSD_vs_E[rmd_index].resize(NumberOfParameterSets);

    int count =0;
    //for (int shgate=shortgate_low; shgate < shortgate_high; shgate=shgate+20) {
    // for (int lgate=longgate_low; lgate < longgate_high; lgate=lgate+100) {
    for (int shgate=shortgate_low; shgate < shortgate_high; shgate=shgate+5) {
      for (int lgate=longgate_low; lgate < longgate_high; lgate=lgate+20) {

	nameStream.str("");
	nameStream<<"RMD"<<RMD<<"-SG"<<shgate<<"LG"<<lgate<<" PSD (100 keVee)";
	RMD_PSD1[rmd_index][count]=new TH1F(nameStream.str().c_str(),"PSD1",PSDBins,PSDlow,PSDhigh);

	nameStream.str("");
	nameStream<<"RMD"<<RMD<<"-SG"<<shgate<<"LG"<<lgate<<" PSD (200 keVee)";
	RMD_PSD2[rmd_index][count]=new TH1F(nameStream.str().c_str(),"PSD2",PSDBins,PSDlow,PSDhigh);

	nameStream.str("");
	nameStream<<"RMD"<<RMD<<"-SG"<<shgate<<"LG"<<lgate<<" PSD (300 keVee)";
	RMD_PSD3[rmd_index][count]=new TH1F(nameStream.str().c_str(),"PSD3",PSDBins,PSDlow,PSDhigh);

	nameStream.str("");
	nameStream<<"RMD"<<RMD<<"-SG"<<shgate<<"LG"<<lgate<<" PSD vs E";
	RMD_PSD_vs_E[rmd_index][count]=new TH2F(nameStream.str().c_str(),"PSD_vs_E",EBins,Elow,Ehigh,PSDBins,PSDlow,PSDhigh);


	count++;
	
      }
    }

    //cout << count << " " << NumberOfParameterSets << endl;

    //RMD_PSD1.push_back(PSD1);
    //PSD1.clear();
  }


  cout<<"DONE Initializing all the histograms"<<endl;
  // set input tree branvh variables and addresses
  ////////////////////////////////////////////////////////////////////////////////////
  
  //Specify the  branch
  LendaEvent* inEvent = new LendaEvent();
  inT->SetBranchAddress("lendaevent",&inEvent);
  //  outT->BranchRef();
  
  LendaEvent *outEvent = new LendaEvent();

 

  clock_t startTime;
  clock_t otherTime;
  double timeRate=0;
  bool timeFlag=true;
  startTime = clock();

  TString previousFile="";
  for (Long64_t jentry=0; jentry<nentry;jentry++) { // Main analysis loop. Go through all entries in the root file
    
    if (inT->GetCurrentFile()->GetName() !=previousFile){ // Check file name to extract run number
      cout<<"NEWFILE"<<endl;
      cout<<"Prev: "<<previousFile<<endl;
      cout<<"Curr: "<<inT->GetCurrentFile()->GetName()<<endl;
      
      previousFile=inT->GetCurrentFile()->GetName();


      string name = previousFile.Data();
      //run-0###-00.root
      cout<<"\n\n\n";
      int index = name.find("run-");
      name=  name.substr(index+4,4);
      int RunNum=atoi(name.c_str());
      thePacker->FindAndSetMapAndCorrectionsFileNames(RunNum); // Use run number to read map and correction files
      // it also builds the array of MapInfo containers GlobalIDToMapInfo with important information for each channel, 
      // like e.g. CFD fit parameters, IDs, trace analysis flag, etc...

      cout << "I just read map file" << endl;
    }



    inT->GetEntry(jentry); // Get the event from the input tree 

    //cout << "Number of bars: " << inEvent->NumBars << endl;


    //if (inEvent->NumBars == 1) { 

    int NumRMD = inEvent->NumBars;
    bool goodRMD = false;
    for (int ibar=0; ibar < NumRMD; ibar++) { // Loop over bars in event

      //int RMDId = inEvent->Bars[0].BarId;
      //string RMDName = inEvent->Bars[0].Name; 
      int RMDId = inEvent->Bars[ibar].BarId;
      string RMDName = inEvent->Bars[ibar].Name; 






      //Double_t E = inEvent->Bars[ibar].Tops[0].GetEnergy();
      Double_t E = inEvent->Bars[ibar].Tops[0].GetPulseHeight();


      goodRMD = true;
      // Energy calibration 
      if (RMDName == "RMD101") { 
	E = 0.6660 * E + 15.624;
      } else if (RMDName == "RMD102") { 
	E = 0.6624 * E + 12.433;
      } else if (RMDName == "RMD103") { 
	E = 0.5977 * E + 14.307;
      } else if (RMDName == "RMD104") { 
	E = 0.7305 * E + 10.015;
      } else if (RMDName == "RMD105") { 
	E = 0.6590 * E + 14.154;
      } else if (RMDName == "RMD106") { 
	E = 0.8348 * E + 14.529;
      } else if (RMDName == "RMD107") { 
	E = 0.6852 * E + 14.295;
      } else if (RMDName == "RMD108") { 
	E = 0.6948 * E + 15.200;
      } else if (RMDName == "RMD109") { 
	E = 0.6842 * E + 11.573;
      } else {
	E = -9990;
	goodRMD = false;
      }

      //cout << "Energy: " << E << endl; 

      
      if (goodRMD) { // Fill histograms only with RMD stuff

		
	LendaFilter aFilter;
	//LendaChannel rmd_ch = inEvent->Bars[0].Tops[0];
	LendaChannel rmd_ch = inEvent->Bars[ibar].Tops[0];
	int size = rmd_ch.GetTrace().size();
	vector <UShort_t> trace = rmd_ch.GetTrace();
	
	int count =0;
	for (int shgate=shortgate_low; shgate < shortgate_high; shgate=shgate+5) {
	  for (int lgate=longgate_low; lgate < longgate_high; lgate=lgate+20) {
	    
	    int start = -1;
	    Double_t sG = aFilter.GetGateRMD(trace,start,shgate);
	    Double_t lG = aFilter.GetGateRMD(trace,start,lgate);
	    
	    if ((lG+sG) > 0) {
	      //double psd = lG/(lG+sG); // PSD defined as RMD 
	      double psd = (lG-sG)/lG; // PSD defined as RMD 
	      if (E > 100) RMD_PSD1[RMDId][count]->Fill(psd);
	      if (E > 200) RMD_PSD2[RMDId][count]->Fill(psd);
	      if (E > 300) RMD_PSD3[RMDId][count]->Fill(psd);
	      
	      RMD_PSD_vs_E[RMDId][count]->Fill(E,psd);
	    }
	    
	    count++;
	    
	  }
	}
	
      }
      
      
    } // End of loop over bars in event
      











    //Periodic printing
    if (jentry % 10000 <10 && jentry >=20000 && timeFlag){
      timeFlag=false;
      otherTime=clock();
      timeRate = TMath::Abs( double((startTime-otherTime))/double(CLOCKS_PER_SEC));
      timeRate = timeRate/jentry;
    }
    //Periodic printing
    if (jentry % 100 ==0 ){
      cout<<flush<<"\r"<<"                                                                                          "<<"\r";
      cout<<"On Event "<<jentry<<" "<<((double)jentry)/(nentry)*100<<"% minutes remaining "<<
	(1.0/60)*timeRate*(nentry-jentry)<<" hours remaining "<<(1.0/3600)*timeRate*(nentry-jentry);
    }



    
    
  }//End main analysis loop
  


  Double_t par[6];
  TF1 *aGauss1 = new TF1("aGauss1","gaus",0.25,0.37);
  TF1 *aGauss2 = new TF1("aGauss2","gaus",0.37,0.50);
  TF1 *aGauss = new TF1("aGauss","gauss(0)+gauss(3)",0.20,0.60);
  aGauss->SetLineColor(4);

  vector <vector<double> > theFOMs(NumberRMD);
  vector <vector<double> > mean1(NumberRMD);
  vector <vector<double> > sigma1(NumberRMD);
  vector <vector<double> > mean2(NumberRMD);
  vector <vector<double> > sigma2(NumberRMD);
  for (auto & a : theFOMs){
    a.resize(NumberOfParameterSets);
  }
  for (auto & a : mean1){
    a.resize(NumberOfParameterSets);
  }
  for (auto & a : sigma1){
    a.resize(NumberOfParameterSets);
  }
  for (auto & a : mean2){
    a.resize(NumberOfParameterSets);
  }
  for (auto & a : sigma2){
    a.resize(NumberOfParameterSets);
  }

  for (int i=0;i<NumberRMD;i++) {
    for (int j=0; j<NumberOfParameterSets;j++) {
      //RMD_PSD1[i][j]->Fit("aGauss1","QSR");
      //RMD_PSD1[i][j]->Fit("aGauss2","QSR+");
      //RMD_PSD2[i][j]->Fit("aGauss1","QSR");
      //RMD_PSD2[i][j]->Fit("aGauss2","QSR+");
      RMD_PSD3[i][j]->Fit("aGauss1","QSR");
      RMD_PSD3[i][j]->Fit("aGauss2","QSR+");
      aGauss1->GetParameters(&par[0]);
      aGauss2->GetParameters(&par[3]);
      aGauss->SetParameters(par);
      //RMD_PSD1[i][j]->Fit("aGauss","QSR+");
      //RMD_PSD2[i][j]->Fit("aGauss","QSR+");
      RMD_PSD3[i][j]->Fit("aGauss","QSR+");
      aGauss->GetParameters(&par[0]);

      mean1[i][j] = par[1]; 
      sigma1[i][j] = par[2]; 
      mean2[i][j] = par[4]; 
      sigma2[i][j] = par[5]; 
      theFOMs[i][j]= (mean2[i][j] - mean1[i][j]) / (2.35*(sigma1[i][j] + sigma2[i][j]));

      cout << "Detector RMD" << i+1 << endl; 
      cout << "Gauss 1: " << mean1[i][j] << " " << sigma1[i][j] << "  Gauss 2:" << mean2[i][j] << " " << sigma2[i][j] << endl; 
      cout << "FOM=" << theFOMs[i][j] << endl;
    }
  }


   
  ofstream out("./TheFOMs-2017test.txt");
  for (int i=0;i<NumberRMD;i++) {
    for (int j=0; j<NumberOfParameterSets;j++) {
      out<<RMD_PSD3[i][j]->GetName()<<"  "<<theFOMs[i][j]<<" "<<mean1[i][j]<<" "<<sigma1[i][j]<<" "<<mean2[i][j]<<" "<<sigma2[i][j]<<endl;
      //out<<RMD_PSD2[i][j]->GetName()<<"  "<<theFOMs[i][j]<<" "<<mean1[i][j]<<" "<<sigma1[i][j]<<" "<<mean2[i][j]<<" "<<sigma2[i][j]<<endl;
    }
  }
  out<<endl;
  
 
  
  //Close the file
  outFile->Write();
  outFile->Close();
  
  cout<<"\n\n**Finished**\n\n";
  
  return  0;
  
}


