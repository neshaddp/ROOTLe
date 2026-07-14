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

  //inT->Add("./rootfiles/run-0406-??.root");
  // inT->Add("./rootfiles/run-0398-00.root");
  // inT->Add("./rootfiles/run-0397-00.root");
  // inT->Add("./rootfiles/run-0396-00.root");
  inT->Add("./../users/jorge/rootfiles/run-4204-00-rmd.root");


  Long64_t nentry=(Long64_t) (inT->GetEntries());
  cout <<"The number of entires is : "<< nentry << endl ;
  



  // Openning output Tree and output file
  ////////////2////////////////////////////////////////////////////////////////////////
  //TFile *outFile = new TFile("Result-4204-L5-6-G0-1-W3-6-D4-7.root","recreate");
  TFile *outFile = new TFile("Result-temp.root","recreate");
 
  //int FL_low=8;
  //int FL_high=9;
  //int FG_low=0;
  //int FG_high=1;
  //int w_low=3;
  //int w_high=4;
  //int d_low=5;
  //int d_high=6;

  int FL_low=5;
  int FL_high=6;
  int FG_low=0;
  int FG_high=1;
  int w_low=5;
  int w_high=6;
  int d_low=5;
  int d_high=6;

 

  int NumberOfFilterSets = (FL_high-FL_low)*(FG_high-FG_low)*(w_high-w_low)*(d_high-d_low);
  cout<<"YOU ASKED FOR "<<NumberOfFilterSets <<" filter sets"<<endl;
  vector <TH1F*> TheHistograms(NumberOfFilterSets);
  vector <TH1F*> TheCubicHistograms(NumberOfFilterSets);
  vector <TH2F*> TheCubic_vs_Energy(NumberOfFilterSets);

  vector <TH2F*> TheCubic_vs_Noise(NumberOfFilterSets);
  vector <TH2F*> TheCubic_vs_SigNoise(NumberOfFilterSets);
  vector <TH2F*> TheNoise_vs_ERMDHistograms(NumberOfFilterSets);

  vector <TH1F*> TheSignalToNoiseRMDHistograms(NumberOfFilterSets);
  vector <TH1F*> TheSignalToNoiseLENDAHistograms(NumberOfFilterSets);
  vector <TH1F*> TheNoiseRMDHistograms(NumberOfFilterSets);
  vector <TH1F*> TheNoiseLENDAHistograms(NumberOfFilterSets);


  map < string, int> MapOfRejectedEvents;

  stringstream nameStream;
  int count =0;
  int xlow=600;
  int xhigh=1000;
  int nBins=8000;
  
  // Loop over all possible timing parameter combinations. Create histograms for each combination
  for (int FL=FL_low;FL<FL_high;FL++){ 
    for (int FG=FG_low;FG<FG_high;FG++){
      for (int w=w_low;w<w_high;w++){
	for (int d=d_low;d<d_high;d++){
	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d;
	  TheHistograms[count]=new TH1F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh);
	  MapOfRejectedEvents[nameStream.str()]=0; // This keeps track of "failed" parameter combinations

	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Cubic";
	  TheCubicHistograms[count]=new TH1F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh);
	  TheCubicHistograms[count]->GetXaxis()->SetTitle("Tdiff (clock ticks)");
	  //	  MapOfRejectedEvents[nameStream.str()]=0;

	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Cubic_vs_Energy";
	  TheCubic_vs_Energy[count]=new TH2F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh,4000,0,TMath::Power(2.0,14.0));
	  TheCubic_vs_Energy[count]->GetXaxis()->SetTitle("Tdiff (clock ticks)");
	  TheCubic_vs_Energy[count]->GetYaxis()->SetTitle("Energy (keVee)");

	  //TheCubic_vs_Energy[count]=new TH2F(nameStream.str().c_str(),"Title",nBins,675,682,4000,0,700);
	  //	  MapOfRejectedEvents[nameStream.str()]=0;
	  
	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Signal-to-Noise";
	  TheSignalToNoiseRMDHistograms[count]=new TH1F(nameStream.str().c_str(),"Title",3000,0,300);
	  TheSignalToNoiseRMDHistograms[count]->GetXaxis()->SetTitle("Signal-to-Noise ratio");
	  
	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Signal-to-Noise Lenda Top";
	  TheSignalToNoiseLENDAHistograms[count]=new TH1F(nameStream.str().c_str(),"Title",30000,0,30000);
	  TheSignalToNoiseLENDAHistograms[count]->GetXaxis()->SetTitle("Signal-to-Noise ratio");
	  
	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Noise";
	  TheNoiseRMDHistograms[count]=new TH1F(nameStream.str().c_str(),"Title",1000,0,100);
	  TheNoiseRMDHistograms[count]->GetXaxis()->SetTitle("Noise (a.u.)");

	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Noise_vs_E";
	  TheNoise_vs_ERMDHistograms[count]=new TH2F(nameStream.str().c_str(),"Title",5000,0,1000,500,0,20);
	  TheNoise_vs_ERMDHistograms[count]->GetXaxis()->SetTitle("Energy (keVee)");
	  TheNoise_vs_ERMDHistograms[count]->GetYaxis()->SetTitle("Noise (a.u.)");
	  
	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Noise Lenda Top";
	  TheNoiseLENDAHistograms[count]=new TH1F(nameStream.str().c_str(),"Title",1000,0,100);
	  TheNoiseLENDAHistograms[count]->GetXaxis()->SetTitle("Noise (a.u.)");

	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Cubic_vs_Noise";
	  //TheCubic_vs_Noise[count]=new TH2F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh,1000,0,20);
	  TheCubic_vs_Noise[count]=new TH2F(nameStream.str().c_str(),"Title",100,677,681,1000,0,20);
	  TheCubic_vs_Noise[count]->GetXaxis()->SetTitle("Tdiff (clock ticks)");
	  TheCubic_vs_Noise[count]->GetYaxis()->SetTitle("Noise (a.u.)");

	  nameStream.str("");
	  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Cubic_vs_sigNoise";
	  //TheCubic_vs_Noise[count]=new TH2F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh,1000,0,300);
	  TheCubic_vs_SigNoise[count]=new TH2F(nameStream.str().c_str(),"Title",100,677,681,500,0,200);
	  TheCubic_vs_SigNoise[count]->GetXaxis()->SetTitle("Tdiff (clock ticks)");
	  TheCubic_vs_SigNoise[count]->GetYaxis()->SetTitle("Signal-to-Noise (a.u.)");

	  count++;
	}
      }
    }
  }

  cout<<"DONE Initializing all the histograms"<<endl;
  // set input tree branvh variables and addresses
  ////////////////////////////////////////////////////////////////////////////////////
  
  //Specify the  branch
  LendaEvent* inEvent = new LendaEvent();
  inT->SetBranchAddress("lendaevent",&inEvent);
  //  outT->BranchRef();
  
  LendaEvent *outEvent = new LendaEvent();

  //non branch timing variables 
     ////////////////////////////////////////////////////////////////////////////////////


  //  Filter theFilter; // Filter object
  ////////////////////////////////////////////////////////////////////////////////////


  //  vector <Sl_Event*> CorrelatedEvents(4,NULL);
  //  map <Long64_t,bool> mapOfUsedEntries;//Used to prevent double counting


  // This are corrections for the time (or cubictime) dependence on Dt (or CubicDt) (where Dt is time diff top - bottom)
  Double_t cubicCor=-0.154953;
  Double_t linCor=-0.153091;

  Double_t cubicCor2=2.74170e-02;
  Double_t linCor2=2.49526e-02;

  

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

    }

    inT->GetEntry(jentry); // Get the event from the input tree 



    // We require to have TWO bars with single-events (for ToF calculation). The bars are "NL01" and "SL12"
    // It is also assumed that the first bar Bar[0] will always be NL01, and the second one will be SL12 (?)
    // Moreover, the PulseHeight in SL12 must be > 2500 (?)
    if (inEvent->NumBars == 2 &&  //inEvent->Bars[0].SimpleEventBit &&inEvent->Bars[1].SimpleEventBit
	inEvent->Bars[0].Name=="RMD104" && inEvent->Bars[1].Name=="SL01") { // && inEvent->Bars[1].Bottoms[0].GetPulseHeight()>2500){



      //Loop over the all the filters in the same way as above
      int count =0;
      for (int FL=FL_low;FL<FL_high;FL++){
	for (int FG=FG_low;FG<FG_high;FG++){
	  for (int w=w_low;w<w_high;w++){
	    for (int d=d_low;d<d_high;d++){

	      // Record LENDA time
	      Double_t lenda_time = 0.5 * (inEvent->Bars[1].Tops[0].GetTime() + inEvent->Bars[1].Bottoms[0].GetTime());
	      Double_t lenda_cubictime = 0.5 * (inEvent->Bars[1].Tops[0].GetCubicTime() + inEvent->Bars[1].Bottoms[0].GetCubicTime());


	      Double_t iniE = inEvent->Bars[0].Tops[0].GetEnergy();
	      Double_t inizerocubic = inEvent->Bars[0].Tops[0].GetCubicCFD();


	      // And now re-define filter paramters for RMDs
	      thePacker->ForceAllFilters(FL,FG,d,w);
	      thePacker->ReMakeLendaEvent(inEvent,outEvent); // Re-make a new LendaEvent using the loop CFD parameters
	                                                     // instead of the original ones defined in the MapInfo container
	      
	      Double_t time = outEvent->Bars[0].Tops[0].GetTime() - lenda_time;
	      Double_t CubicTime = outEvent->Bars[0].Tops[0].GetCubicTime() - lenda_cubictime;
	      Double_t E = outEvent->Bars[0].Tops[0].GetEnergy();
	      if (outEvent->Bars[0].Name=="RMD104") {
		E = 0.7104 * E - 2.6486;
	      }

	      TheHistograms[count]->Fill(time);
	      
	      TheCubicHistograms[count]->Fill(CubicTime);
	      
	      //TheCubic_vs_Energy[count]->Fill(CubicTime,outEvent->Bars[0].Tops[0].GetPulseHeight());
	      TheCubic_vs_Energy[count]->Fill(CubicTime,E);





	      // Calculate now trace parameters for RMD
	      LendaFilter aFilter;
	      LendaChannel rmdch = outEvent->Bars[0].Tops[0];
	      int size = rmdch.GetTrace().size();
	      vector <UShort_t> trace = rmdch.GetTrace();

	      int MaxSpotInTrace_temp = -1;
	      Double_t rising = -1;
	      Double_t background = -1;
	      Double_t Noise = -1;
	      Double_t amplitude = aFilter.GetPulseComplete(trace,MaxSpotInTrace_temp, rising, background, Noise);	      
	      Double_t sigNoise = amplitude / Noise;

	      //if (sigNoise > 60) TheCubicHistograms[count]->Fill(CubicTime);


	      TheSignalToNoiseRMDHistograms[count]->Fill(sigNoise);
	      TheNoiseRMDHistograms[count]->Fill(Noise);
	      TheCubic_vs_Noise[count]->Fill(CubicTime,Noise);
	      TheCubic_vs_SigNoise[count]->Fill(CubicTime,sigNoise);

	      TheNoise_vs_ERMDHistograms[count]->Fill(E,Noise);

	      // Calculate now trace parameters for LENDA top
	      LendaChannel lendach = outEvent->Bars[1].Tops[0];
	      int sizelenda = lendach.GetTrace().size();
	      vector <UShort_t> lendatrace = lendach.GetTrace();

	      MaxSpotInTrace_temp = -1;
	      rising = -1;
	      background = -1;
	      Noise = -1;
	      amplitude = aFilter.GetPulseComplete(lendatrace,MaxSpotInTrace_temp, rising, background, Noise);	      
	      sigNoise = amplitude / Noise;
	      TheSignalToNoiseLENDAHistograms[count]->Fill(sigNoise);
	      TheNoiseLENDAHistograms[count]->Fill(Noise);
	      






	      //cout << "Amplitude: " << amplitude << " Maximum: " << MaxSpotInTrace_temp << " Rasing time: " << rising << 
	      //" Baseline: " << background << " Noise: " << Noise <<  " Signal-to-noise: " << sigNoise << endl;

	      outEvent->Finalize();


	      outEvent->Clear();
	      count++;
	    }
	  }
	}
      }

      
      
    }
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
      cout<<"On Event "<<jentry<<" "<<((double)jentry)/(nentry)*100<<"% minutes remaining "<<(1.0/60)*timeRate*(nentry-jentry)<<" hours remaining "<<(1.0/3600)*timeRate*(nentry-jentry);
    }
    
    
  }//End main analysis loop
  

  vector < vector <TH1F*> > theVecs;
  //theVecs.push_back(TheHistograms);
  theVecs.push_back(TheCubicHistograms);
  
  int NumOfHistVectors=theVecs.size();

  cout << "Num of HistVectors: " <<  NumOfHistVectors << endl;

  TF1 * aGauss = new TF1("aGauss","gaus",660,690);
  TFitResultPtr result;
  Int_t status;
  vector <vector<double> > theResolutions(NumOfHistVectors);
  vector <vector<double> > theChi2s(NumOfHistVectors);

  for (auto & a : theResolutions){
    a.resize(NumberOfFilterSets);
  }

  for (auto & a : theChi2s){
    a.resize(NumberOfFilterSets);
  }
 
  double BestRes=1000;
  string BestName="NONE";
  double BestFit=1000;
  string BestFitName="NONE";

  for (int i=0;i<NumOfHistVectors;i++){
    for (int j=0;j<NumberOfFilterSets;j++){
      result = theVecs[i][j]->Fit("aGauss","QSR");
      status=result;
      if (status==0) {
  	//theResolutions[i][j]=result->Value(2)*2.35*4/TMath::Sqrt(2);
	Double_t totfwhm = result->Value(2)*2.35*4; // RMD FWHM in ns
 	theResolutions[i][j]=TMath::Sqrt(pow(totfwhm,2) - 0.25); //RMD FWHM in ns after subtracting 500-ps (squared) resultion from LENDA

  	cout<<"Res is "<<theResolutions[i][j]<<endl;
  	if (theResolutions[i][j]<BestRes){
  	  BestRes=theResolutions[i][j];
  	  BestName=theVecs[i][j]->GetName();
  	}
  	theChi2s[i][j]=result->Chi2()/result->Ndf();
  	if (result->Chi2()/result->Ndf() < BestFit){
  	  BestFit=result->Chi2()/result->Ndf();
  	  BestFitName=theVecs[i][j]->GetName();
  	}
      } else {
  	cout<<"***Bad Fit For "<<theVecs[i][j]->GetName()<<endl;
      }
    }
  }
  
  ofstream out("./TheResoultions-temp.txt");
  for (int i=0;i<theResolutions.size();i++){
    for (int j=0;j<theResolutions[i].size();j++){
      out<<theVecs[i][j]->GetName()<<"  "<<theResolutions[i][j]<<MapOfRejectedEvents[theVecs[i][j]->GetName()]<<"  "<<theChi2s[i][j]<<endl;
    }
  }
  out<<endl;
  
  
  //Close the file
  outFile->Write();
  outFile->Close();
  cout<<"Best Res is "<<BestName<<" "<<BestRes<<endl;
  cout<<"Best Fit is "<<BestFitName<<" "<<BestFit<<endl;
  
  cout<<"\n\n**Finished**\n\n";
  
  return  0;
  
}


