#define Analoop_cxx

#include "Analoop.h"


/////////////////////////////////////////////////
// This method is called once at the beginning //
// it is where histograms are allocated        //
// ie when new is called for them	       //
/////////////////////////////////////////////////
void Analoop::SlaveBegin(TTree * t) {
  TString option = GetOption();

  

  // Set signal handler
  signal_received = kFALSE;
  signal(SIGINT, signalhandler);


  
  //////////////////////////////////////////////////////////////////
  // Look for a settings object in the Input list to the selector //
  // if it is not there then all is lost			  //
  //////////////////////////////////////////////////////////////////

  TheSettings = (R00TLeSettings*) fInput->FindObject("Settings0");

  //  TheSettings=(R00TLeSettings*)gDirectory->Get("TheSettings");
  if (TheSettings != NULL){
    cout<<"Found A R00TLeSettings Object!"<<endl;
  }else{
    cout<<"Could Not Find a R00TLeSettings Object.  Hard Exit"<<endl;
    exit(1);
  }

  PreviousTreeNumber=-9999;//Set tree number to an intial bad value

  TString user =gSystem->Getenv("R00TLe_User");
  TString install =gSystem->Getenv("R00TLeInstall");
  TString CorrectionsFileName = install+"/users/"+user+"/AllTheCorrections.txt";


  MyCorrections = new TEnv(CorrectionsFileName);
  
  CurrentRunNumber =-1;
  
  TString S800MapFile= install+"/users/"+user+"/invmap.inv";
  loadInverseMap(S800MapFile,&MapA,&MapY,&MapB,&MapD);

  ////////////////////////////////////////////////////////////////////
  // make histograms for standard quantities (TOFs and PulsHeights) //
  ////////////////////////////////////////////////////////////////////
  
  Int_t NumBars =-1;
  NumBars = TheSettings->GetNumBars();//The number of bars in this file from TheSettings
  
  //Resize the vectors for the standard quanties for LendaBars
  AvgEnergies.resize(NumBars);
  AvgTOFs.resize(NumBars);
  TopTOFs.resize(NumBars);
  BottomTOFs.resize(NumBars);
  BottomEnergies.resize(NumBars);
  TopEnergies.resize(NumBars);
  TopETimeDiff.resize(NumBars);
  BotETimeDiff.resize(NumBars);
  GeoAvg.resize(NumBars);
  TopPHTimeDiff.resize(NumBars);
  BotPHTimeDiff.resize(NumBars);
  GeoAvgPH.resize(NumBars);
  TopPulseHeight.resize(NumBars);
  BotPulseHeight.resize(NumBars);
  AvgPulseHeight.resize(NumBars);
  Top_PSD_vs_E.resize(NumBars);
  Bot_PSD_vs_E.resize(NumBars);

  //Binning of the TOF and energy histograms 
  int EBins=5000; 
  int EMax=100000;
  int TOFBins=3000;
  int TOFMinMax=100;
  double PSDmin = -0.2;  
  double PSDmax = 0.6; //default 0.6
  int PSDbins = 400;
  int MaxE = 100000;
  int MinE = 0;
  int Ebins = 10000;




  //For each Bar in the file make these standar histograms
  stringstream nameStream;
  for (int i=0;i<NumBars;i++){
     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_AvgE";
     AvgEnergies[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);

     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_TopE";
     TopEnergies[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);

     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_BottomE";
     BottomEnergies[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);

     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_TopTOF";
     TopTOFs[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax);

     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_BottomTOF";
     BottomTOFs[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax);

     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_AvgTOF";
     AvgTOFs[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_TopE_vs_TimeDiff";
     TopETimeDiff[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax,EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_BottomE_vs_TimeDiff";
     BotETimeDiff[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax,EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_GeoAvg";
     GeoAvg[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax,EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_AvgPH";
     AvgPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_TopPH";
     TopPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);

     nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_BottomPH";
     BotPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_TopE_vs_TimeDiff_using_PH";
     TopPHTimeDiff[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax,EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_BottomE_vs_TimeDiff_using_PH";
     BotPHTimeDiff[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax,EBins,0,EMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_GeoAvg_using_PH";
     GeoAvgPH[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-TOFMinMax,TOFMinMax,EBins,0,EMax);
     
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_Top_PSD_vs_E";
     Top_PSD_vs_E[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),Ebins,MinE,MaxE,PSDbins,PSDmin,PSDmax);
     
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_Bot_PSD_vs_E";
     Bot_PSD_vs_E[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),Ebins,MinE,MaxE,PSDbins,PSDmin,PSDmax);
     
  }
  
 /* 
  // ///////////////////////////////
  // // new more histograms here  //
  ///////////////////////////////
  double PSDMin = -0.2;  
  double PSDMax = 0.5; //default 0.6
  int PSDBins = 400;
  int MaxCalibE = 100000;
  int MinCalibE = 0;
  int CalibEBins = 10000;
  double PSDMin2 = -.2;  
  double PSDMax2 = .6;
  int PSDBins2 = 200;
  int NumLongGates;
  int NumShortGates;
  TString LongGates[] = {"Max+20","Max+30","Max+40","Max+50","Max+60","Max+70"};
  TString ShortGates[] = {"Max-4","Max",};
  
  NumLongGates = sizeof(LongGates)/sizeof(LongGates[0]);
  NumShortGates = sizeof(ShortGates)/sizeof(ShortGates[0]);
  
  RMD_PSD1_vs_E.resize(NumBars);
  RMDPSD1.resize(NumBars);
  RMD_PSD2_vs_E.resize(NumBars);
  RMDPSD2.resize(NumBars);
  PSD1Cut1.resize(NumBars);
  PSD1Cut2.resize(NumBars);
  PSD1Cut3.resize(NumBars);
 
 
 for (int i=0;i<NumBars;i++){

  RMD_PSD1_vs_E[i].resize(NumShortGates);
  RMDPSD1[i].resize(NumShortGates);
  RMD_PSD2_vs_E[i].resize(NumShortGates);
  RMDPSD2[i].resize(NumShortGates);
  PSD1Cut1[i].resize(NumShortGates);
  PSD1Cut2[i].resize(NumShortGates);
  PSD1Cut3[i].resize(NumShortGates);
	 
	for (int shortgate=0;shortgate<NumShortGates;shortgate++)
	{
		RMD_PSD1_vs_E[i][shortgate].resize(NumLongGates);
		RMDPSD1[i][shortgate].resize(NumLongGates);
		RMD_PSD2_vs_E[i][shortgate].resize(NumLongGates);
		RMDPSD2[i][shortgate].resize(NumLongGates);
		PSD1Cut1[i][shortgate].resize(NumLongGates);
		PSD1Cut2[i][shortgate].resize(NumLongGates);
		PSD1Cut3[i][shortgate].resize(NumLongGates);
		
		for(int longate = 0; longate<NumLongGates;longate++)
		{
			nameStream.str("");
			nameStream<<TheSettings->GetBarName(i)<<"E vs PSD1-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate];
			RMD_PSD1_vs_E[i][shortgate][longate] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),CalibEBins,MinCalibE,MaxCalibE,PSDBins,PSDMin,PSDMax);
			
			
			nameStream.str("");
			nameStream<<TheSettings->GetBarName(i)<<"PSD1"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] <<"";
			RMDPSD1[i][shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
			
			//nameStream.str("");
			//nameStream<<"PSD2"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] <<"";
			//RMDPSD2[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
						
			//nameStream.str("");
			//nameStream<<"E vs PSD2: Short gate: "<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate];
			//RMD_PSD2_vs_E[shortgate][longate] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),CalibEBins,MinCalibE,MaxCalibE,PSDBins2,PSDMin2,PSDMax2);
						
			//nameStream.str("");
			//nameStream<<TheSettings->GetBarName(i)<<"PSD1Cut1-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] <<"(";
			//PSD1Cut1[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
					
			//nameStream.str("");
			//nameStream<<"PSD1Cut2-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] << "";
			//PSD1Cut2[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
			
			//nameStream.str("");
			//nameStream<<"PSD1Cut3-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate]<<"";
			//PSD1Cut3[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
			
			
		}
	}
 }

*/
  // Add all the histograms to Output list of the Selector
  // fOutput is a TSelectorList which will take ownership of
  // an object when it is added
  fOutput->AddAll(gDirectory->GetList());
}




//////////////////////////////////////////////////////////////////
// This is where you fill histograms and make cuts on the data. //
// Process is called for every entry in the ROOT tree	        //
// 							        //
// The tree contains two branches s800calc and lendaevent       //
//////////////////////////////////////////////////////////////////
Bool_t Analoop::Process(Long64_t entry) {
  //  cout<<"IN PROCESS "<<entry<<endl;
  fChain->GetTree()->GetEntry(entry);
  
  if (signal_received) {
    Info(__FUNCTION__, "At the entry No. %lld", entry);
    this->Abort("Signal received.", kAbortProcess);
  }

  ///////////////////////////////////////////////////////////////////////
  // The following catch is to figure out if the file has been changed //
  // since the last time time process was called.  PROOF will process  //
  // the files in an arbitrary order so one needs to check	       //
  ///////////////////////////////////////////////////////////////////////
 if (fChain->GetTreeNumber() != PreviousTreeNumber){
    PreviousTreeNumber= fChain->GetTreeNumber();

    TString temp=fChain->GetCurrentFile()->GetName();
    string name = temp.Data(); 

    int index = name.find("run-");
    name=  name.substr(index+4,4);
    int RunNum=atoi(name.c_str());

    CurrentRunNumber=RunNum;
  } 


  //Nice load Bar feature from S.Noji
  loadBar(entry, nentries, 1000, 50);
  
                                 /////////////////////////
                                 // Begin s800 analysis //
                                 /////////////////////////

  double IC = s800calc->GetIC()->GetSum();
  double TOF=0;
  if (s800calc->GetMultiHitTOF()->fObj.size() !=0){
    TOF=s800calc->GetMultiHitTOF()->GetFirstObjHit();
  }

  TOF+=1000;//Make the TOF a more reasonable Positive number
  Hist("PID",TOF,IC,8000,0,3000,1000,0,2000);


  //Extract CRDC information
  Float_t crdc1Xcog = s800calc->GetCRDC(0)->GetXcog();
  Float_t crdc2Xcog = s800calc->GetCRDC(1)->GetXcog();
  
  Float_t crdc1TAC =  s800calc->GetCRDC(0)->GetTAC();
  Float_t crdc2TAC =  s800calc->GetCRDC(1)->GetTAC();

  //Get the mask calibrations for the current run
  Double_t Crdc1XOffset = MyCorrections->GetValue(Form("run%04d.crdc1.xOffset",CurrentRunNumber),0.0);
  Double_t Crdc1YOffset = MyCorrections->GetValue(Form("run%04d.crdc1.yOffset",CurrentRunNumber),0.0);
  Double_t Crdc1YSlope = MyCorrections->GetValue(Form("run%04d.crdc1.ySlope",CurrentRunNumber),0.0);
  
  Double_t Crdc2XOffset = MyCorrections->GetValue(Form("run%04d.crdc2.xOffset",CurrentRunNumber),0.0);
  Double_t Crdc2YOffset = MyCorrections->GetValue(Form("run%04d.crdc2.yOffset",CurrentRunNumber),0.0);
  Double_t Crdc2YSlope = MyCorrections->GetValue(Form("run%04d.crdc2.ySlope",CurrentRunNumber),0.0);

  //Form the calibrated CRDC information
  Float_t xfp =crdc1Xcog*2.54 + Crdc1XOffset; 
  Float_t yfp =crdc1TAC *Crdc1YSlope + Crdc1YOffset;

  Float_t x_crdc2 = crdc2Xcog*2.54 +Crdc2XOffset;
  Float_t y_crdc2 = crdc2TAC*Crdc2YSlope +Crdc2YOffset;

  Float_t afp = TMath::ATan((x_crdc2 - xfp)/1073.0);
  Float_t bfp = TMath::ATan((y_crdc2-yfp)/1073.0);

  //Prefom the ray tracing
  //The ray tracing routine needs input in meters and radian
  Float_t ata,yta,bta,dta;
  ata=yta=bta=dta=0;
  Raytracing(MapA,MapY,MapB,MapD, xfp/1000. ,afp ,yfp/1000, bfp ,ata,yta,bta,dta);
  //Returns in Meters and raidian
  //Put the angles in degrees
  //Put yta in mm
  //put Di in percent
  ata*=(180./3.14159265);
  bta*=(180./3.14159265);
  yta*=1000;
  dta*=100;

  //Put in Degrees
  afp*=(180./3.14159265);
  bfp*=(180./3.14159265);

  Hist("afp",afp,1000,-100,100);
  Hist("bfp",bfp,1000,-100,100); 
  Hist("xfp",xfp,1000,-400,400);
  Hist("yfp",yfp,1000,-400,400);
    
  Hist("dta",dta,1000,-100,100);
  Hist("ata",ata,1000,-100,100);
  Hist("bta",bta,1000,-100,100);
  Hist("yta",yta,1000,-100,100);
  


  
  //////////////////////////////////////////////////////
  // Each lendaevent has variable number of bars      //
  // Each Bar has varibale number of tops and bottoms //
  // Loop over all of them below                      //
  //////////////////////////////////////////////////////
  
  
  
  int NumBarsInEvent= lendaevent->NumBars; // Get them number of bars in the event
  
  for (int i =0;i<NumBarsInEvent;i++){ //Loop over all Bars in the Lenda event
    
    int NumTopsInBar = lendaevent->Bars[i].NumTops; //Number of tops in ith bar
    int NumBottomsInBar = lendaevent->Bars[i].NumBottoms;//Number of bottoms in ith bar
    int BarId = lendaevent->Bars[i].BarId; //The ith Bars Unique Bar Id

    for (int t=0;t<NumTopsInBar;t++){ // Loop over all the TOPS in the Bar
      double E=lendaevent->Bars[i].Tops[t].GetEnergy();
	  double PH=lendaevent->Bars[i].Tops[t].GetPulseHeight();
      if (E>0){
	TopEnergies[BarId]->Fill(E); //Fill Top Energy if >0
      }
	  if (PH>0){
	TopPulseHeight[BarId]->Fill(PH); //Fill Top Energy if >0
      }
    }//End for over tops
    
    for (int b=0;b<NumBottomsInBar;b++){//Loop over all the Bottoms in the Bar
      double E=lendaevent->Bars[i].Bottoms[b].GetEnergy();
	  double PH=lendaevent->Bars[i].Bottoms[b].GetPulseHeight();
      if (E>0){
	BottomEnergies[BarId]->Fill(E); //Fill Bottom Energy if >0
      }
	  if (PH>0){
	BotPulseHeight[BarId]->Fill(PH); //Fill Top Energy if >0
      }
    }//End for over bottoms
	
	double TopE;
	double BotE;
	double TopPH;
	double BotPH;
	double PSD_top = -999;
	double PSD_bot = -999;

	
	 for (int k=0; k<NumTopsInBar; k++){ // Loop over all the TOPS in the Bar
      TopE=lendaevent->Bars[i].Tops[k].GetEnergy();
	  TopPH=lendaevent->Bars[i].Tops[k].GetPulseHeight();
	  PSD_top = (1-lendaevent->Bars[i].Tops[k].GetShortGate()/lendaevent->Bars[i].Tops[k].GetLongGate());
	  //calibrated_E = lendaevent->Bars[i].Tops[k].GetEnergy() * rmd_slope + rmd_intercept;
	  //rmd_ch = lendaevent->Bars[i].Tops[k];
	  //refPSD = 1.-lendaevent->Bars[i].Tops[k].GetShortGate()/lendaevent->Bars[i].Tops[k].GetLongGate();
		
	 }
	 for (int m=0; m<NumBottomsInBar; m++){ // Loop over all the TOPS in the Bar
	  BotE=lendaevent->Bars[i].Bottoms[m].GetEnergy();
	  BotPH=lendaevent->Bars[i].Bottoms[m].GetPulseHeight();
	  PSD_bot = (1-lendaevent->Bars[i].Bottoms[m].GetShortGate()/lendaevent->Bars[i].Bottoms[m].GetLongGate());
	 }
	 
	 for (int n=0; n< NumTopsInBar; n++){ // Loop over all the TOPS in the Bar
	 
	  double TopTOF=lendaevent->Bars[i].GetTopTOF();
	  double BottomTOF=lendaevent->Bars[i].GetBottomTOF();
	  
	  double lendaTDiff = TopTOF-BottomTOF;
	  double Geoavg = sqrt(TopE*BotE);
	  double GeoavgPH = sqrt(TopPH*BotPH);
	  
      if (TopE>0){		  
	   TopETimeDiff[BarId]->Fill(lendaTDiff,TopE); //Fill Top Energy if > 0
      }
	  if (BotE>0){		  
	   BotETimeDiff[BarId]->Fill(lendaTDiff,BotE); //Fill Top Energy if > 0
      }
	  if (Geoavg>0){		  
	  GeoAvg[BarId]->Fill(lendaTDiff,Geoavg); //Fill Top Energy if > 0
      }
	  
	  if (TopPH>0){		  
	   TopPHTimeDiff[BarId]->Fill(lendaTDiff,TopPH); //Fill Top Energy if > 0
      }
	  if (BotPH>0){		  
	   BotPHTimeDiff[BarId]->Fill(lendaTDiff,BotPH); //Fill Top Energy if > 0
      }
	  if (GeoavgPH>0){		  
	  GeoAvgPH[BarId]->Fill(lendaTDiff,GeoavgPH); //Fill Top Energy if > 0
      }
	  
	  if (TopE>0){		  
	   Top_PSD_vs_E[BarId]->Fill(TopE,PSD_top); //Fill Top Energy if > 0
      }
	  if (BotE>0){		  
	   Bot_PSD_vs_E[BarId]->Fill(BotE,PSD_bot); //Fill Top Energy if > 0
      }
    }//End for over tops
	

    //Fill histograms that are only valid when there was 1 top and 1 bottom
    if (lendaevent->Bars[i].SimpleEventBit){//has 1 top and 1 bottom
      if (lendaevent->Bars[i].Tops[0].GetEnergy() > 0 && lendaevent->Bars[i].Bottoms[0].GetEnergy()>0){
	AvgEnergies[BarId]->Fill(lendaevent->Bars[i].GetAvgEnergy());
      }
	  if (lendaevent->Bars[i].Tops[0].GetPulseHeight() > 0 && lendaevent->Bars[i].Bottoms[0].GetPulseHeight()>0){
	AvgPulseHeight[BarId]->Fill(lendaevent->Bars[i].GetAvgPulseHeight());
      }
	  
      AvgTOFs[BarId]->Fill(lendaevent->Bars[i].GetAvgTOF());
    }
    
    if ( NumTopsInBar == 1 ){
      TopTOFs[BarId]->Fill(lendaevent->Bars[i].GetTopTOF());
    }
    if (NumBottomsInBar == 1 ){
      BottomTOFs[BarId]->Fill(lendaevent->Bars[i].GetBottomTOF());
    }
	
  }//End for over bars in event
  
 
 

  for (int i = 0; i<NumBarsInEvent;i++) //this is set up to do timing for
  {                                    //the LSD/RMD detectors. Will now play with energy cuts      

		  		   
    int NumTopsInBar = lendaevent->Bars[i].NumTops; //Number of tops in ith bar
    int NumBottomsInBar = lendaevent->Bars[i].NumBottoms;//Number of bottoms in ith bar
    int BarId = lendaevent->Bars[i].BarId; //The ith Bars Unique Bar Id
    string temprefname = lendaevent->Bars[i].GetBarName();
	double E = lendaevent->Bars[i].Tops[0].GetEnergy();
	double E_rmd = -999;
	double E_lsc = -999;
	double PH = -999; 
	double PSD = -999;
	double PSD_top1 = -999;
	double PSD_bot1 = -999;
	double PSD_top2 = -999;
	double PSD_bot2 = -999;
	double PH_calibrated_RMD = -999;
	double PH_lsc = -999;
	double RMD_PSD = -999;
	double PH_calibrated_LSC = -999;
	double lenda_tof = -999;
	string detector = "RMD4";
	bool run = true;
	bool calibrated = true;
	double rmd_softtime = -999.;
	double rmd_cubictime = -999.;
	//cout << " " << endl;
	double rmd_slope = 1;
	double rmd_intercept = 0;
	TString LENDAname = "NL01";
	double cubicToF = -999.;
	long double lenda_top = -999; ///lenda top time cubic
    long double lenda_bottom = -999; ///////lenda bottom timme cubic
	Double_t sG = -999;
	Double_t lG = -999;
	
/*	
	if(temprefname=="NL01") //This for PSD analysis with default gates changing
	{
		double e_rmd_top = lendaevent->Bars[i].Tops[0].GetEnergy();
		double e_rmd_bot = lendaevent->Bars[i].Bottoms[0].GetEnergy();
		PSD_top1 = (1-lendaevent->Bars[i].Tops[0].GetShortGate()/lendaevent->Bars[i].Tops[0].GetLongGate());
		PSD_bot1 = (1-lendaevent->Bars[i].Bottoms[0].GetShortGate()/lendaevent->Bars[i].Bottoms[0].GetLongGate());
		AutoHisto("NL01 Top: PSD vs E", e_rmd_top, PSD_top1, 10000, 0, 100000, 400,-.2,.6);
		AutoHisto("NL01 Bottom: PSD vs E", e_rmd_bot, PSD_bot1, 10000, 0, 100000, 400,-.2,.6);
		AutoHisto("NL01 Top PSD All",  PSD_top1, 400, -.2,.6);
		AutoHisto("NL01 Bottom PSD All",  PSD_bot1, 400, -.2,.6);
	}
	
	if(temprefname=="NL03") //This for PSD analysis with default gates changing
	{
		double e_rmd_top = lendaevent->Bars[i].Tops[0].GetEnergy();
		double e_rmd_bot = lendaevent->Bars[i].Bottoms[0].GetEnergy();
		PSD_top2 = (1-lendaevent->Bars[i].Tops[0].GetShortGate()/lendaevent->Bars[i].Tops[0].GetLongGate());
		PSD_bot2 = (1-lendaevent->Bars[i].Bottoms[0].GetShortGate()/lendaevent->Bars[i].Bottoms[0].GetLongGate());
		AutoHisto("NL03 Top: PSD vs E", e_rmd_top, PSD_top2, 10000, 0, 100000, 400,-.2,.6);
		//AutoHisto("NL03 Bottom: PSD vs E", e_rmd_bot, PSD_bot2, 10000, 0, 100000, 400,-.2,.6);
		AutoHisto("NL03 Top PSD All",  PSD_top2, 400, -.2,.6);
		//AutoHisto("NL03 Bottom PSD All",  PSD_bot2, 400, -.2,.6);
	}
*/
/*
  for (int t=0;t<NumTopsInBar;t++){
	//if(temprefname=="NL01") {	 //This for PSD analysis with manual gates changing
	
		LendaFilter aFilter;
		//LendaChannel rmd_ch = inEvent->Bars[0].Tops[0];
		LendaChannel rmd_ch = lendaevent->Bars[i].Tops[t];
		double calibrated_E = lendaevent->Bars[i].Tops[t].GetEnergy() * rmd_slope + rmd_intercept;
		int size = rmd_ch.GetTrace().size();
		vector <UShort_t> trace = rmd_ch.GetTrace();
		int Max = -999;
		double Max_val = 0;
		double refPSD = 1.-lendaevent->Bars[i].Tops[t].GetShortGate()/lendaevent->Bars[i].Tops[t].GetLongGate();
	
					
		//if(calibrated_E >3875 && calibrated_E <5850){
					
			//	if(refPSD >.14 && refPSD < .17){
			//	cout << "this is a neutron "<< entry<<" and has pulse height "<<calibrated_E <<endl;
			//}
		
			//if(refPSD > -.02 && refPSD <.01){
			//	cout << "this is a gamma "<< entry<<" and has pulse height "<<calibrated_E <<endl; 
			//}
				
			for(int j = 0; j<size; j++)/////Finding the max of the trace to base integration on
			{
				if(trace[j] > Max_val) {
				
					Max = j;	
					Max_val = trace[j];
				
				}// end of if
			}//end of for
				
		double shortgate_high = Max;
		double shortgate_low = Max - 4;
		double longgate_high = Max -4 + 70;
		double longgate_low = Max -4 +20;
		
		int k = 0;
		int count =0;
		
		for (int shgate=shortgate_low; shgate < shortgate_high+1; shgate=shgate+4) 
		{
			int l = 0;
			int start = shgate;
			for (int lgate=longgate_low; lgate < longgate_high+1; lgate=lgate+10) 
			{
				double L = lgate - shgate;
				
				sG = aFilter.GetGate(trace,start,10); //Integrates trace from "raising" point (calculated from trace derivative) over range L. 
				lG = aFilter.GetGate(trace,start,L);  //Integrates trace from "raising" point (calculated from trace derivative) over range L. 				
				
				if((lG + sG) > 0)
				{
					double psd1 = (lG-sG)/lG;
										
					//if(calibrated_E >4000 && calibrated_E <6000){
						
						RMD_PSD1_vs_E[BarId][k][l]->Fill(calibrated_E,psd1);		
						RMDPSD1[BarId][k][l]->Fill(psd1);
						//PSD1Cut1[k][l]->Fill(psd1);
						//}
						
					//if(calibrated_E >3000 && calibrated_E <8000) {PSD1Cut2[k][l]->Fill(psd1);}
					//if(calibrated_E >3000 && calibrated_E <10000){PSD1Cut2[k][l]->Fill(psd1);}															
				}
				l++;
				//if(lgate == longgate_high) {cout << psd1 <<endl;	}
			}//end of for loop of long gate
			k++;
									
		}//end of for loop of short gate
		//} end of if loop of energy cut	
	 //} //end of if statement
   }
*/	
	/*
	if (lendaevent->NumBars == 8 && temprefname == LENDAname) //This is for plots of E vs LendaDiff
	{
		double topE = lendaevent->Bars[i].Tops[0].GetEnergy();
		double botE = lendaevent->Bars[i].Bottoms[0].GetEnergy();
		double E_calibrated = -999;
		
		
		//if (topE > 3000 && botE > 4000) /////making an energy cut, don't have top and bottom calibrated but this should be fine using 3k and 4k respectively {
			
			for(int j = 0; j < NumBarsInEvent;j++)
			{
				double E = lendaevent->Bars[j].Tops[j].GetEnergy();
				string tempname2 = lendaevent->Bars[j].GetBarName();
				if (tempname2 == LENDAname)
				{
				   topE = lendaevent->Bars[j].Tops[j].GetEnergy();
				   botE = lendaevent->Bars[j].Bottoms[j].GetEnergy();
				   
				   //lenda_tof = .5 * (lendaevent->Bars[j].Tops[0].GetCubicTime() + lendaevent->Bars[j].Bottoms[0].GetCubicTime());
				   lenda_top = (long double)lendaevent->Bars[j].Tops[j].GetCubicTime();
				   lenda_bottom = (long double)lendaevent->Bars[j].Bottoms[j].GetCubicTime();
				}
				else
				
			
				//if(tempname2 =="NL02")
				{
					PH = lendaevent->Bars[j].Tops[j].GetEnergy();
					PH_calibrated_RMD = PH * rmd_slope + rmd_intercept;
					E_rmd = E;
					rmd_softtime = lendaevent->Bars[j].Tops[j].GetSoftTime();
					rmd_cubictime = (long double)lendaevent->Bars[j].Tops[j].GetCubicTime();
					//PSD = (1-lendaevent->Bars[j].Tops[0].GetShortGate()/lendaevent->Bars[j].Tops[0].GetLongGate());
					//AutoHisto("PSD All",  PSD, 400, -.2,.6);
				}//end else
				
				
				
			}//end inner for

			//cubicToF = (rmd_cubictime - lenda_tof);
			long double lendaTDiff = lenda_top-lenda_bottom;
			double GeoAvg = sqrt(topE * botE);
			
			//Energy_vs_PSD[BarId]->Fill(E_calibrated, PSD);
			AutoHisto("Lenda Time Diff",lendaTDiff,2000, -20, 20);
			AutoHisto("Lenda Top time",lenda_top,10000, -20, 20);
			AutoHisto("Lenda Bottom time",lenda_bottom,10000, -20, 20);
			AutoHisto("RMD Top time",rmd_cubictime,2000, -20, 20);
			AutoHisto("RMD ToF vs E",lendaTDiff,PH_calibrated_RMD,2000, -20, 20, 10000,0,20000);
			AutoHisto("RMD Top: ToF vs Energy", lendaTDiff, topE, 2000, -100, 100, 10000,0,100000);
			AutoHisto("RMD Bottom: ToF vs Energy", lendaTDiff, botE, 2000, -100, 100, 10000,0,100000);
			AutoHisto("Geometrical Average", lendaTDiff,GeoAvg, 2000, -100, 100, 10000,0,100000);
			
			
				if (lendaTDiff >=1.5 && lendaTDiff <=1.9){
				
				
				AutoHisto("Lenda Time Diff cut1",lendaTDiff,2000, -20, 20);
				AutoHisto("RMD ToF vs E cut1",lendaTDiff,PH_calibrated_RMD,2000, -20, 20, 10000,0,16000);
				AutoHisto("NL01T ToF vs Energy cut1", lendaTDiff, topE, 2000, -100, 100, 10000,0,20000);
				AutoHisto("NL01B ToF vs Energy cut1", lendaTDiff, botE, 2000, -100, 100, 10000,0,20000);
				}
			
				else if (lendaTDiff >=1.9 && lendaTDiff <=2.3){
			
				AutoHisto("Lenda Time Diff cut2",lendaTDiff,2000, -20, 20);
				AutoHisto("RMD ToF vs E cut2",lendaTDiff,PH_calibrated_RMD,2000, -20, 20, 10000,0,16000);
				AutoHisto("NL01T ToF vs Energy cut2", lendaTDiff, topE, 2000, -100, 100, 10000,0,20000);
				AutoHisto("NL01B ToF vs Energy cut2", lendaTDiff, botE, 2000, -100, 100, 10000,0,20000);
				}
		
				else if (lendaTDiff >=2.3 && lendaTDiff <=2.7){
			
				AutoHisto("Lenda Time Diff cut3",lendaTDiff,2000, -20, 20);
				AutoHisto("RMD ToF vs E cut3",lendaTDiff,PH_calibrated_RMD,2000, -20, 20, 10000,0,16000);
				AutoHisto("NL01T ToF vs Energy cut3", lendaTDiff, topE, 2000, -100, 100, 10000,0,20000);
				AutoHisto("NL01B ToF vs Energy cut3", lendaTDiff, botE, 2000, -100, 100, 10000,0,20000);
			
				}; //end of inner else statement
		    }//end of inner if statement
         */
	//}//end of outer if statement 

	
  }//End of for loop 

  return kTRUE;
}








///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/////////////////The casual user can ignore the vodoo magic below//////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

Analoop::~Analoop() {
  // Destructor
  Info(__FUNCTION__,"Destructing.");  
  delete s800calc;
  delete lendaevent;
}

void Analoop::Init(TTree *tree) {
  s800calc    = new S800Calc();
  lendaevent  = new LendaEvent();

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(0);
  fChain->SetBranchAddress("s800calc",    &s800calc,    &b_s800calc);
  fChain->SetBranchAddress("lendaevent",  &lendaevent,  &b_lendaevent);

  
  accum_nentries = 0;
}

Bool_t Analoop::Notify() {
  if (fChain->GetCurrentFile() != NULL){
    
    filename = fChain->GetCurrentFile()->GetName();

    Info("Notify","File: %s (%7.2f MB)", filename.Data(), fChain->GetTree()->GetTotBytes()/1024./1024.);

    treenum = fChain->GetTreeNumber() + 1;
    entrynum = fChain->GetChainEntryNumber(0);

    nentries = fChain->GetTree()->GetEntries();
    accum_nentries += nentries;
    Info("Notify","%lld entries (%6.2f%%) in this tree",nentries, (double)nentries / fChain->GetEntries() * 100.);
    Info("Notify","%lld entries (%6.2f%%) up to this tree",accum_nentries, (double)accum_nentries / fChain->GetEntries() * 100.);
  }
  return kTRUE;
}

void Analoop::Begin(TTree *) {
  TString option = GetOption();

  OutFileName=option;
}


void Analoop::SlaveTerminate() {

}

void Analoop::Terminate() {

  TheSettings = (R00TLeSettings*) fInput->FindObject("Settings0");

  cout<<"\n Writing histograms to "<<OutFileName<<endl;
  TFile out(OutFileName,"recreate");
  fOutput->Write();//Write all the histograms to disk
  TheSettings->Write();

  out.Close();

  this->ResetAbort();
  signal_received = kFALSE;
}



/**Make a 1D histogram with a Name

 */
void Analoop::MakeHistogram(TString name,Int_t bins,Double_t xlow,Double_t xhigh){
  
  fOutput->AddLast( new TH1F(name,name,bins,xlow,xhigh));

}

/**Fill a 1D histogram with a name

 */
void Analoop::FillHistogram(TString name,Float_t value){

  TObject * object = fOutput->FindObject(name);

  if (object == NULL){
    Error("Analoop::FillHistogram",name+" not found");
    return;
  }
  TString className=object->ClassName();
  if (className !="TH1F"){
    Error("Analoop::FillHistogram",name+" not a histogram");
  }

  ((TH1F*)object)->Fill(value);

}



/**Makea 2D histogram with a Name
 */
void Analoop::MakeHistogram(TString name,Int_t binsX,Double_t xlow,Double_t xhigh,Int_t binsY,Double_t yLow,Double_t yHigh){
  fOutput->AddLast( new TH2F(name,name,binsX,xlow,xhigh,binsY,yLow,yHigh));
}


/**Fill a 2D histogram with a Name
 */
void Analoop::FillHistogram(TString name,Float_t Xvalue,Float_t Yvalue){

  TObject * object = fOutput->FindObject(name);
  
  if (object == NULL){
    Error("Analoop::FillHistogram",name+" not found");
    return;
  }
  TString className=object->ClassName();
  if (className !="TH2F"){
    Error("Analoop::FillHistogram",name+" not a histogram");
  }

  ((TH2F*)object)->Fill(Xvalue,Yvalue);

}




void Analoop::AutoHisto(TString name,Float_t value,Int_t bins, Double_t xlow, Double_t xhigh){
  
  TObject * object = fOutput->FindObject(name);
  if ( object == NULL) {//The histogram is not there 
    MakeHistogram(name,bins,xlow,xhigh);
  }
  FillHistogram(name,value);
}

void Analoop::AutoHisto(TString name,Float_t Xvalue,Float_t Yvalue,Int_t binsX, Double_t xlow, Double_t xhigh,Int_t binsY,Double_t yLow,Double_t yHigh){

  TObject * object = fOutput->FindObject(name);
  if ( object == NULL) {//The histogram is not there 
    MakeHistogram(name,binsX,xlow,xhigh,binsY,yLow,yHigh);
  }
  FillHistogram(name,Xvalue,Yvalue);
  
}


void Analoop::AutoHisto(Int_t HistNumber,Float_t value,Int_t bins, Double_t xlow, Double_t xhigh){
  stringstream s;
  s<<"h"<<HistNumber;
  AutoHisto(s.str().c_str(),value,bins,xlow,xhigh);

}
void Analoop::AutoHisto(Int_t HistNumber,Float_t Xvalue,Float_t Yvalue,Int_t binsX, Double_t xlow, Double_t xhigh,Int_t binsY,Double_t ylow,Double_t yhigh){
  stringstream s;
  s<<"h"<<HistNumber;
  AutoHisto(s.str().c_str(),Xvalue,Yvalue,binsX,xlow,xhigh,binsY,ylow,yhigh);

}


