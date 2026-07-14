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

  inT->Add("./../users/jorge/rootfiles/run-4605-00-rmd.root");


  Long64_t nentry=(Long64_t) (inT->GetEntries());
  cout <<"The number of entires is : "<< nentry << endl ;
  



  // Openning output Tree and output file
  ////////////2////////////////////////////////////////////////////////////////////////
  //TFile *outFile = new TFile("Result-test2017-RMD101-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD102-Efunction.root","recreate");
  TFile *outFile = new TFile("Result-test2017-RMD103-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD104-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD105-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD106-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD107-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD108-Efunction.root","recreate");
  //TFile *outFile = new TFile("Result-test2017-RMD109-Efunction.root","recreate");



  int NumberOfEnergies = 10;
  vector <TH1F*> TheCubicHistogramsE(NumberOfEnergies);
  TH1F* TheCubicHistograms;
  TH2F* TheCubic_vs_Energy;
  TH1F* TheNoiseRMDHistograms;


  map < string, int> MapOfRejectedEvents;

  stringstream nameStream;
  int count =0;
  int xlow=-25;
  int xhigh=-5;
  int nBins=1000;
  //int nBins=200;
  int ERMDBins=1000; 
  int ERMDMax=2000;
 
  
  int FL=2;
  int FG=4;
  int w=6;
  int d=3;

  
  for (int E=60; E < 360; E=E+30) {
    nameStream.str("");
    nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<"_Cubic-E-"<<E;
    TheCubicHistogramsE[count]=new TH1F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh);
    TheCubicHistogramsE[count]->GetXaxis()->SetTitle("Tdiff (clock ticks)");
    
    count++;
  }

  nameStream.str("");
  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Cubic";
  TheCubicHistograms=new TH1F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh);
  TheCubicHistograms->GetXaxis()->SetTitle("Tdiff (clock ticks)");


  nameStream.str("");
  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Cubic_vs_Energy";
  TheCubic_vs_Energy=new TH2F(nameStream.str().c_str(),"Title",nBins,xlow,xhigh,ERMDBins,0,ERMDMax);
  TheCubic_vs_Energy->GetXaxis()->SetTitle("Tdiff (clock ticks)");
  TheCubic_vs_Energy->GetYaxis()->SetTitle("Energy (keVee)");


  nameStream.str("");
  nameStream<<"FL"<<FL<<"FG"<<FG<<"w"<<w<<"d"<<d<<" Noise";
  TheNoiseRMDHistograms=new TH1F(nameStream.str().c_str(),"Title",1000,0,100);
  TheNoiseRMDHistograms->GetXaxis()->SetTitle("Noise (a.u.)");



  cout<<"DONE Initializing all the histograms"<<endl;
  // set input tree branvh variables and addresses
  ////////////////////////////////////////////////////////////////////////////////////
  
  //Specify the  branch
  LendaEvent* inEvent = new LendaEvent();
  inT->SetBranchAddress("lendaevent",&inEvent);
  
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

    }

    inT->GetEntry(jentry); // Get the event from the input tree 



    // We require to have TWO bars with single-events (for ToF calculation): one from RMD104 and another from NL10
    // It is also assumed that the first bar Bar[0] will always be NL10
    if (inEvent->NumBars == 2 &&  inEvent->Bars[0].SimpleEventBit &&
	//inEvent->Bars[1].Name=="RMD101" && inEvent->Bars[0].Name=="NL10") { 
	//inEvent->Bars[1].Name=="RMD102" && inEvent->Bars[0].Name=="NL10") { 
	inEvent->Bars[1].Name=="RMD103" && inEvent->Bars[0].Name=="NL10") { 
        //inEvent->Bars[1].Name=="RMD104" && inEvent->Bars[0].Name=="NL10") { 
	//inEvent->Bars[1].Name=="RMD105" && inEvent->Bars[0].Name=="NL10") { 
	//inEvent->Bars[1].Name=="RMD106" && inEvent->Bars[0].Name=="NL10") { 
	//inEvent->Bars[1].Name=="RMD107" && inEvent->Bars[0].Name=="NL10") { 
	//inEvent->Bars[1].Name=="RMD108" && inEvent->Bars[0].Name=="NL10") { 
	//inEvent->Bars[1].Name=="RMD109" && inEvent->Bars[0].Name=="NL10") { 

      int count;
      
      int FL=2;
      int FG=4;
      int w=6;
      int d=3;

      // Record LENDA time
      Double_t lenda_time = 0.5 * (inEvent->Bars[0].Tops[0].GetTime() + inEvent->Bars[0].Bottoms[0].GetTime());
      Double_t lenda_cubictime = 0.5 * (inEvent->Bars[0].Tops[0].GetCubicTime() + inEvent->Bars[0].Bottoms[0].GetCubicTime());
      Double_t lenda_E = 0.5 * (inEvent->Bars[0].Tops[0].GetEnergy() + inEvent->Bars[0].Bottoms[0].GetEnergy());
      
      
      // And now re-define filter paramters for RMDs
      thePacker->ForceAllFilters(FL,FG,d,w);
      thePacker->ReMakeLendaEvent(inEvent,outEvent); // Re-make a new LendaEvent using the loop CFD parameters
      // instead of the original ones defined in the MapInfo container
      
      Double_t CubicTime = outEvent->Bars[1].Tops[0].GetCubicTime() - lenda_cubictime;
      Double_t E = outEvent->Bars[1].Tops[0].GetPulseHeight();
      
      if (outEvent->Bars[1].Name == "RMD101") { 
	E = 0.6660 * E + 15.624;
      } else if (outEvent->Bars[1].Name == "RMD102") { 
	E = 0.6624 * E + 12.433;
      } else if (outEvent->Bars[1].Name == "RMD103") { 
	E = 0.5977 * E + 14.307;
      } else if (outEvent->Bars[1].Name == "RMD104") { 
	E = 0.7305 * E + 10.015;
      } else if (outEvent->Bars[1].Name == "RMD105") { 
	E = 0.6590 * E + 14.154;
      } else if (outEvent->Bars[1].Name == "RMD106") { 
	E = 0.8348 * E + 14.529;
      } else if (outEvent->Bars[1].Name == "RMD107") { 
	E = 0.6852 * E + 14.295;
      } else if (outEvent->Bars[1].Name == "RMD108") { 
	E = 0.6948 * E + 15.200;
      } else if (outEvent->Bars[1].Name == "RMD109") { 
	E = 0.6842 * E + 11.573;
      }
      
      
      
      if (lenda_E > 5000) { // Plot events with LENDA energies in the "good" region (i.e., where the 500-ps resolution makes sense
	// (using lower energies would include LENDA resolutions > 500 ps)

	
	
	TheCubicHistograms->Fill(CubicTime);
	TheCubic_vs_Energy->Fill(CubicTime,E);
	
	
	//Calculate now trace parameters for RMD
	LendaFilter aFilter;
	LendaChannel rmdch = outEvent->Bars[1].Tops[0];
	int size = rmdch.GetTrace().size();
	vector <UShort_t> trace = rmdch.GetTrace();
	
	int MaxSpotInTrace_temp = -1;
	Double_t rising = -1;
	Double_t background = -1;
	Double_t Noise = -1;
	Double_t amplitude = aFilter.GetPulseComplete(trace,MaxSpotInTrace_temp, rising, background, Noise);	      
	Double_t sigNoise = amplitude / Noise;
	
	TheNoiseRMDHistograms->Fill(Noise);

	

	//cout <<  "RMD104 energy: " << E << " count: " << count << endl;
 
	//Double_t Ebin = 5; // Define slize of E-projetion on Cubic
	Double_t Ebin = 15;
	
	if (E >=50-Ebin && E<=50+Ebin) {
	  count = 0;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=90-Ebin && E<=90+Ebin) {
	  count = 1;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=120-Ebin && E<=120+Ebin) {
	  count = 2;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=150-Ebin && E<=150+Ebin) {
	  count = 3;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=180-Ebin && E<=180+Ebin) {
	  count = 4;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=210-Ebin && E<=210+Ebin) {
	  count = 5;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=240-Ebin && E<=240+Ebin) {
	  count = 6;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=270-Ebin && E<=270+Ebin) {
	  count = 7;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=300-Ebin && E<=300+Ebin) {
	  count = 8;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	} else if (E >=330-Ebin && E<=330+Ebin) {
	  count = 9;
	  TheCubicHistogramsE[count]->Fill(CubicTime);
	}
	

	//cout << " I did it" << endl;

	
      }
      
      
      outEvent->Finalize();
           
      outEvent->Clear();
      
      
      
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
  
  

  cout << "Num of Histograms: " <<  NumberOfEnergies << endl;

  TF1 * aGauss = new TF1("aGauss","gaus",-20,-6);
  TFitResultPtr result;
  Int_t status;
  vector <double> theResolutions(NumberOfEnergies);
  vector <double> theChi2s(NumberOfEnergies);

  
  double BestRes=1000;
  string BestName="NONE";
  double BestFit=1000;
  string BestFitName="NONE";

  for (int i=0;i<NumberOfEnergies;i++){
    result = TheCubicHistogramsE[i]->Fit("aGauss","QSR");
    status=result;
    if (status==0) {
      Double_t totfwhm = result->Value(2)*2.35*4; // RMD FWHM in ns
      theResolutions[i]=TMath::Sqrt(pow(totfwhm,2) - 0.25); //RMD FWHM in ns after subtracting 500-ps (squared) resultion from LENDA

      cout<<"Res is "<<theResolutions[i]<<endl;
      if (theResolutions[i]<BestRes){
	BestRes=theResolutions[i];
	BestName=TheCubicHistogramsE[i]->GetName();
      }

      theChi2s[i]=result->Chi2()/result->Ndf();
      if (result->Chi2()/result->Ndf() < BestFit){
	BestFit=result->Chi2()/result->Ndf();
	BestFitName=TheCubicHistogramsE[i]->GetName();
      }
    } else {
      cout<<"***Bad Fit For "<<TheCubicHistogramsE[i]->GetName()<<endl;
    }
  }

  
  //ofstream out("./TheResoultions-energy-RMD101.txt");
  //ofstream out("./TheResoultions-energy-RMD102.txt");
  ofstream out("./TheResoultions-energy-RMD103.txt");
  //ofstream out("./TheResoultions-energy-RMD104.txt");
  //ofstream out("./TheResoultions-energy-RMD105.txt");
  //ofstream out("./TheResoultions-energy-RMD106.txt");
  //ofstream out("./TheResoultions-energy-RMD107.txt");
  //ofstream out("./TheResoultions-energy-RMD108.txt");
  //ofstream out("./TheResoultions-energy-RMD109.txt");
  for (int i=0;i<theResolutions.size();i++){
    out<< TheCubicHistogramsE[i]->GetName() << " Res: " << theResolutions[i] << " Chi2: " <<theChi2s[i] << endl;
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


