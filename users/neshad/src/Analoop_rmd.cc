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
  Int_t NumRMDS = 1;
  
  //CHANGE THIS FOR OTHER ONES
  
  cout << "Number of bars in file: " << NumBars << endl;  
  
  //Resize the vectors for the standard quanties for LendaBars
  AvgEnergies.resize(NumBars);
  BottomEnergies.resize(NumBars);
  TopEnergies.resize(NumBars);
  TopPulseHeight.resize(NumBars);
  BottomPulseHeight.resize(NumBars);
  AvgPulseHeight.resize(NumBars);
  CalibratedEnergy.resize(NumBars);

  
  
  //Binning of the TOF and energy histograms 
  int EBins=30000; 
  int EMax=300000;
  int EMaxCal = 600;
  int EMinCal = 0;
  int EBinsCal = 600;
  int PHBins = 4250;
  int PHMax = 17000;
  int PHMin = 0;
  double PSDMin = -0.2;  
  double PSDMax = .6;
  int PSDBins = 200;
  
  //For each Bar in the file make these standard histograms
  stringstream nameStream;
  string chname;
  int barchid;

  
  
  for (int i=0;i<NumBars;i++){
	  
		 chname = TheSettings->GetBarName(i);
		 barchid = TheSettings->GetBarId(chname);
		 
		 
		cout << "Bar name: " << chname << " has ID number: " << barchid << endl;  
		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_AvgE";
		AvgEnergies[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);


		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_AvgPH";
		AvgPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PHBins,0,PHMax);

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Calibrated_E";
		CalibratedEnergy[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBinsCal,EMinCal,EMaxCal);

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_TopE";
		TopEnergies[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);


		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Top_PH";
		TopPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PHBins,PHMin,PHMax);

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Bottom_PH";
		BottomPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PHBins,PHMin,PHMax);


		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_BottomE";
		BottomEnergies[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),EBins,0,EMax);

  }

  // ///////////////////////////////
  // // new more histograms here  //
  ///////////////////////////////
  
  ////Attempting to use this idea for PSD analysis, want to make different figures w/ various long and short gates
  ////will define what these gates are, value-wise, at some point
  int MaxCalibE = 1400;
  int MinCalibE = 0;
  int CalibEBins = 700;
  double PSDMin2 = .4;  
  double PSDMax2 = 1;
  int PSDBins2 = 200;
  int NumLongGates;
  int NumShortGates;
  TString LongGates[] = {"Max+20","Max+30","Max+40","Max+50"};
  TString ShortGates[] = {"Max-4","Max",};
  
  NumLongGates = sizeof(LongGates)/sizeof(LongGates[0]);
  NumShortGates = sizeof(ShortGates)/sizeof(ShortGates[0]);
  
  RMD_PSD1_vs_E.resize(NumShortGates);
  RMDPSD1.resize(NumShortGates);
  RMD_PSD2_vs_E.resize(NumShortGates);
  RMDPSD2.resize(NumShortGates);
  PSD1Cut1.resize(NumShortGates);
  PSD1Cut2.resize(NumShortGates);
  PSD1Cut3.resize(NumShortGates);
  PSD2Cut1.resize(NumShortGates);
  PSD2Cut2.resize(NumShortGates);
  PSD2Cut3.resize(NumShortGates);
  
	for (int shortgate=0;shortgate<NumShortGates;shortgate++)
	{
		RMD_PSD1_vs_E[shortgate].resize(NumLongGates);
		RMDPSD1[shortgate].resize(NumLongGates);
		RMD_PSD2_vs_E[shortgate].resize(NumLongGates);
		RMDPSD2[shortgate].resize(NumLongGates);
		PSD1Cut1[shortgate].resize(NumLongGates);
		PSD1Cut2[shortgate].resize(NumLongGates);
		PSD1Cut3[shortgate].resize(NumLongGates);
		PSD2Cut1[shortgate].resize(NumLongGates);
		PSD2Cut2[shortgate].resize(NumLongGates);
		PSD2Cut3[shortgate].resize(NumLongGates);
		
		
		for(int longate = 0; longate<NumLongGates;longate++)
		{
			nameStream.str("");
			nameStream<<"E vs PSD1-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate];
			RMD_PSD1_vs_E[shortgate][longate] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),CalibEBins,MinCalibE,MaxCalibE,PSDBins,PSDMin,PSDMax);
			
			
			nameStream.str("");
			nameStream<<"PSD1-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] <<"(100keV)";
			PSD1Cut1[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
			
			/*
			nameStream.str("");
			nameStream<<"E vs PSD2: Short gate: "<<ShortGates[shortgate]<<"; Long gate: " <<LongGates[longate];
			RMD_PSD2_vs_E[shortgate][longate] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),CalibEBins,MinCalibE,MaxCalibE,PSDBins2,PSDMin2,PSDMax2);
			*/
			
			nameStream.str("");
			nameStream<<"PSD1-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] << "(200keV)";
			PSD1Cut2[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
			
			nameStream.str("");
			nameStream<<"PSD1-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate]<<"(500keV)";
			PSD1Cut3[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin,PSDMax);
			
			nameStream.str("");
			nameStream<<"PSD2-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] <<"(100keV)";
			PSD2Cut1[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin2,PSDMax2);
			
			nameStream.str("");
			nameStream<<"PSD2-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate] << "(200keV)";
			PSD2Cut2[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin2,PSDMax2);
			
			nameStream.str("");
			nameStream<<"PSD2-SG-"<<ShortGates[shortgate]<<"-LG-" <<LongGates[longate]<<"(500keV)";
			PSD2Cut3[shortgate][longate] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PSDBins,PSDMin2,PSDMax2);
		}
	}

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
  //cout<<"IN PROCESS "<<entry<<endl;
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

  /*
  Hist("afp",afp,1000,-100,100);
  Hist("bfp",bfp,1000,-100,100); 
  Hist("xfp",xfp,1000,-400,400);
  Hist("yfp",yfp,1000,-400,400);
    
  Hist("dta",dta,1000,-100,100);
  Hist("ata",ata,1000,-100,100);
  Hist("bta",bta,1000,-100,100);
  Hist("yta",yta,1000,-100,100);
  */


  
  //////////////////////////////////////////////////////
  // Each lendaevent has variable number of bars      //
  // Each Bar has variable number of tops and bottoms //
  // Loop over all of them below                      //
  //////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //  This is a basic loop, fills bottom and top energy and pulse height plus the averages.  ////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  

 
  int NumBarsInEvent= lendaevent->NumBars;
  double slope = 1;
  double intercept = 0;
  string LENDAname = "NL01";
  string RMD1name = "RMD1";
  string RMD2name = "RMD2";
  string RMD3name = "RMD3";
  string RMD4name = "RMD4";
  string RMD5name = "RMD5";
  string RMD6name = "RMD6";
  string RMD7name = "RMD7";
  string RMD8name = "RMD8";
  string RMD9name = "RMD9";
  string LSC1name = "LSC1";
  string LSC2name = "LSC2";
  string refdet = "LSC2";
  double calibrated_E = -999.;
  bool Cf_source = true;
  double PSD = -999.;
	//TString runnum = "5103";
  double shortgate_low = 0;
  double shortgate_high = 0;
  vector<int> short_gates;
  vector<int>long_gates;
  
  for (int i = 0; i<NumBarsInEvent; i++)
  {
	bool goodchannel = false;
	 	 //cout << "it is inside the loop"<< endl;
	string BarName = lendaevent->Bars[i].GetBarName(); //bar name
	int NumTopsInBar = lendaevent->Bars[i].NumTops; //Number of tops in ith bar
	int NumBottomsInBar = lendaevent->Bars[i].NumBottoms;//Number of bottoms in ith bar
	int BarId = lendaevent->Bars[i].BarId; //The ith Bars Unique Bar Id	  
	

	if (BarName == RMD1name)    {slope = .0230337; intercept = -71.69704;} //Energy calibrations for the bars/samples. Mostly doing this here to plot calibrated energies
	if (BarName == RMD4name)	{slope = .023552; intercept = -61.9719;}
	if (BarName == RMD5name) 	{slope = .0553628; intercept = -154.29709;}
	if (BarName == RMD6name) 	{slope = .0535289; intercept = -156.9073;}
	if (BarName == RMD7name)	{slope = .1829949; intercept = -125.93706;}		////    if (BarName == RMD7name)	{slope = .1738190; intercept = -46.930903;}
	if (BarName == RMD8name)	{slope = .1211519; intercept = -79.40985;}
	if (BarName == RMD9name)	{slope = .1633322; intercept = -38.418701; goodchannel = true;}
//	if (BarName == LSC1name)	{slope = .021405; intercept = -56.53591;}		//////slope = .02474467; intercept = -103.72709;
//	if (BarName == LSC2name)	{slope = .072599; intercept = -38.97435;}		//////slope = .088313; intercept = -99.8487;
	 	
		
	/////This is where I will go through each of the traces and adjust the PSD manually by calling the LendaFilter class
	if(goodchannel == true)
	{		
		LendaFilter aFilter;
		//LendaChannel rmd_ch = inEvent->Bars[0].Tops[0];
		LendaChannel rmd_ch = lendaevent->Bars[i].Tops[0];
		calibrated_E = lendaevent->Bars[i].Tops[0].GetPulseHeight() * slope + intercept;
		int size = rmd_ch.GetTrace().size();
		vector <UShort_t> trace = rmd_ch.GetTrace();
		int Max = -999;
		double Max_val = 0;
		
		for(int j = 0; j<size; j++)/////Finding the max of the trace to base integration on
		{
			if(trace[j] > Max_val) {Max = j;	Max_val = trace[j];}
		}
		shortgate_high = Max;
		double longgate_high = Max -4 + 50;
		double longgate_low = Max -4 +20;
		
		shortgate_low = Max - 4;
		int k = 0;
		int count =0;
		
		for (int shgate=shortgate_low; shgate < shortgate_high+1; shgate=shgate+4) 
		{
			int l = 0;
			int start = shgate;
			for (int lgate=longgate_low; lgate < longgate_high+1; lgate=lgate+10) 
			{
				double L = lgate - shgate;
				Double_t sG = aFilter.GetGate(trace,start,10);
				Double_t lG = aFilter.GetGate(trace,start,L);
				if((lG + sG) > 0)
				{
					double psd1 = (lG-sG)/lG;
					double psd2 = lG/(lG+sG);
					if(calibrated_E > 100 && calibrated_E< 1400) {RMD_PSD1_vs_E[k][l]->Fill(calibrated_E,psd1); PSD1Cut1[k][l]->Fill(psd1); PSD2Cut1[k][l]->Fill(psd2);}
					if(calibrated_E > 200 && calibrated_E< 1400) {PSD1Cut2[k][l]->Fill(psd1); PSD2Cut2[k][l]->Fill(psd2);}
					if(calibrated_E > 500 && calibrated_E< 1400) {PSD1Cut3[k][l]->Fill(psd1); PSD2Cut3[k][l]->Fill(psd2);}
					
					
				}
				l++;
				//if(lgate == longgate_high) {cout << psd1 <<endl;	}
			}
			k++;
		}
	
	
	}
	
	 for (int b=0;b<NumBottomsInBar;b++)
	{//Loop over all the Bottoms in the Bar
		double E=lendaevent->Bars[i].Bottoms[b].GetEnergy();
		double PH=lendaevent->Bars[i].Bottoms[b].GetPulseHeight();
		if (E>0){
			BottomEnergies[BarId]->Fill(E); //Fill Bottom Energy if >0
			BottomPulseHeight[BarId]->Fill(PH);
			}
	}//End for over bottoms
	  
	for (int t=0;t<NumTopsInBar;t++)
	{
		double E=lendaevent->Bars[i].Tops[t].GetEnergy();
		double PH=lendaevent->Bars[i].Tops[t].GetPulseHeight();
		//CalibratedEnergy[BarId]->Fill(PH*slope + intercept);
		if (E>0){
			TopEnergies[BarId]->Fill(E); //Fill top Energy if >0
			TopPulseHeight[BarId]->Fill(PH);
			CalibratedEnergy[BarId]->Fill(PH*slope + intercept);
			}
	}

	 /*
    //Fill histograms that are only valid when there was 1 top and 1 bottom
    if (lendaevent->Bars[i].SimpleEventBit){//has 1 top and 1 bottom
      if (lendaevent->Bars[i].Tops[0].GetEnergy() > 0 && lendaevent->Bars[i].Bottoms[0].GetEnergy()>0){
	AvgEnergies[BarId]->Fill(lendaevent->Bars[i].GetAvgEnergy());
	AvgPulseHeight[BarId]->Fill(lendaevent->Bars[i].GetAvgPulseHeight());
      }
      //AvgTOFs[BarId]->Fill(lendaevent->Bars[i].GetAvgTOF());
    }
    */
    if ( NumTopsInBar == 1 ){
     // TopTOFs[BarId]->Fill(lendaevent->Bars[i].GetTopTOF());
    }
    if (NumBottomsInBar == 1 ){
      //BottomTOFs[BarId]->Fill(lendaevent->Bars[i].GetBottomTOF());
    }
	  
	  
	  
	  	  
	 if (BarName.find("RMD") == 0 && Cf_source == true)
	{
		//cout << BarName << endl;
		PSD = (1-lendaevent->Bars[i].Tops[0].GetShortGate()/lendaevent->Bars[i].Tops[0].GetLongGate());
		AutoHisto(BarName + " PSD All Events", PSD, 200, -.1, .3);
		double tempE = lendaevent->Bars[i].Tops[0].GetPulseHeight();
		calibrated_E = lendaevent->Bars[i].Tops[0].GetPulseHeight() * slope + intercept;
		if(calibrated_E > 100){
			AutoHisto(BarName + "PSD-E<100keV", PSD, 200,-.2,.6);
			AutoHisto(BarName + " Energy vs PSD All Events, GetGate default", calibrated_E, PSD, 700,0,1400,200,-.2,.6);
		}
		double psd_threshold_200 = .146;
		double psd_threshold_500 = .14;
		if (calibrated_E > 500 && calibrated_E < 505)
		{
			/*
			if (PSD > 0.1091 && PSD < 0.1197)
				
			{
				cout << "" <<entry << " is a gamma and has pulse height: " << tempE<< endl;
			}
			*/
			/*
			if (PSD > 0.1534 && PSD < 0.1684)
			{
				cout << "" << endl;
				cout << "" <<entry << " is a neutron and has energy: " <<tempE << endl;
				cout << "" << endl;
			}*/
			/*
			AutoHisto(BarName + " All events over 200 keV PSD v E", calibrated_E, PSD,350,0,1400,200,-.1,.3);
			AutoHisto(BarName + " All events over 200 keV PSD", PSD, 200, -.1, .3);
			
			if (PSD > psd_threshold_200)
			{
				AutoHisto(BarName + "_over200_neutrons", PSD, 200, -.1, .3);
			}
			else
			{
				AutoHisto(BarName + "_over200_gammas", PSD, 200, -.1, .3);
			}
			
			
			
			if (calibrated_E > 500)
			{
				AutoHisto(BarName + " All events over 500 keV PSD v E", calibrated_E, PSD,350,0,1400,200,-.1,.3);
				AutoHisto(BarName +" All events over 500 keV PSD", PSD, 200, -.1, .3);
				if (PSD > psd_threshold_500)
				{
					AutoHisto(BarName + "_over500_neutrons", PSD, 200, -.1, .3);
					AutoHisto(BarName + "_over200_neutrons", PSD, 200, -.1, .3);
				}
				else
				{
					AutoHisto(BarName + "_over500_gammas", PSD, 200, -.1, .3);
					AutoHisto(BarName + "_over200_gammas", PSD, 200, -.1, .3);
				}
			}
			*/

		}
		
	}
	
	/*
	if (NumBarsInEvent == 2 && BarName.find("LSC") == 0)
	{
		double refE = lendaevent->Bars[i].Tops[0].GetPulseHeight();
		double refToF = lendaevent->Bars[i].Tops[0].GetCubicTime();
		double refPSD = 1.-lendaevent->Bars[i].Tops[0].GetShortGate()/lendaevent->Bars[i].Tops[0].GetLongGate();
		AutoHisto("LSC PSD", refPSD, 300, -.2,.6);
		AutoHisto("LSC E vs PSD", refE, refPSD, 1000, 0, 17000, 300, -.2,.6);

		if(refE > 4000 && refE <4500)
		{
			if(refPSD >.14 && refPSD < .17)
			{
			cout << "this is a neutron "<< entry<<" and has pulse height "<<refE <<endl;
			}
			if(refPSD > -.02 && refPSD <.01)
			{
			cout << "this is a gamma "<< entry<<" and has pulse height "<<refE <<endl;
			}
		}
		
		double RMDE = -9999;
		double RMDToF = -9999;
		double RMDPSD = -9999;
		double tof = -9999;
		
		for (int j = 0; j <NumBarsInEvent; j++)
		{
			string tempname = lendaevent->Bars[j].GetBarName();
			if (tempname.find("RMD") == 0)
			{
				RMDToF = lendaevent->Bars[j].Tops[0].GetCubicTime();
				RMDE = lendaevent->Bars[j].Tops[0].GetPulseHeight();
				RMDPSD = 1-lendaevent->Bars[j].Tops[0].GetShortGate()/lendaevent->Bars[j].Tops[0].GetLongGate();
				//cout <<"got here" <<endl;
			}
		}
		
		//cout << "refToF: " <<refToF<<endl;
		//cout << "RMDToF: " <<RMDToF<<endl;
		
		tof = refToF-RMDToF;
		//cout <<"Tdiff: "<<tof<<endl;
		AutoHisto("RMD9 PSD vs cubic ToF", tof, RMDPSD, 500,-10,10,300,-.1,.4);
		
	}
	*/
	
	
	if (lendaevent->NumBars == 2 && BarName == LENDAname && lendaevent->Bars[i].SimpleEventBit) //This is for the timing measurements, can also throw neutron efficiency measurements in this loop
	{
		
		
		double topE = topE = lendaevent->Bars[i].Tops[0].GetPulseHeight(); 
		double botE = lendaevent->Bars[i].Bottoms[0].GetPulseHeight(); 
		double lenda_tof = lenda_tof = .5 * (lendaevent->Bars[i].Tops[0].GetCubicTime() +
				   lendaevent->Bars[i].Bottoms[0].GetCubicTime());
		double PH_calibrated_RMD = -999.0;
		double rmd_cubictime = -999.0;
		double cubicToF = -999;
		string nameforgraph = "none";
		
		if (topE > 3710 && botE > 4620) /////making an energy cut, don't have top and bottom calibrated but this should be fine using 3k and 4k respectively
		{
			
			for(int j = 0; j < NumBarsInEvent;j++)
			{
				double E = lendaevent->Bars[j].Tops[0].GetEnergy();
				string tempname = lendaevent->Bars[j].GetBarName();
				if (tempname == RMD1name)    {slope = .0230337; intercept = -71.69704;} 
				if (tempname == RMD4name)	{slope = .023552; intercept = -61.9719;}
				if (tempname == RMD5name) 	{slope = .0553628; intercept = -154.29709;}
				if (tempname == RMD6name) 	{slope = .0535289; intercept = -156.9073;}
				if (tempname == RMD7name)	{slope = .1738190; intercept = -46.930903;}
				if (tempname == RMD8name)	{slope = .1211519; intercept = -79.40985;}
				if (tempname == RMD9name)	{slope = .1633322; intercept = -38.418701;}
				if (tempname == LSC1name)	{slope = .021405; intercept = -56.53591;}		
				if (tempname == LSC2name)	{slope = .072599; intercept = -38.97435;}
				if(tempname.find("RMD") == 0)
				{
					PH_calibrated_RMD = lendaevent->Bars[j].Tops[0].GetPulseHeight() * slope + intercept;
					double lG = lendaevent->Bars[j].Tops[0].GetLongGate();
					double sG = lendaevent->Bars[j].Tops[0].GetShortGate();
					PSD = 1.-sG/lG;
					//rmd_softtime = lendaevent->Bars[j].Tops[0].GetSoftTime();
					rmd_cubictime = lendaevent->Bars[j].Tops[0].GetCubicTime();
					nameforgraph = tempname;
				}//end if
				
				
				
			}//end inner for loop over bars

			cubicToF = (rmd_cubictime - lenda_tof);
			
			//cout << "Calibrated energy: " << PH_calibrated_RMD<< endl;
			
			double bin_begin = -4;
			double bin_end = 6;
			
			if (PH_calibrated_RMD > 0.00 && PH_calibrated_RMD < 50.0) AutoHisto("Bin_01", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 50.0 && PH_calibrated_RMD < 75.0) AutoHisto("Bin_02", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 75.0 && PH_calibrated_RMD < 100.0) AutoHisto("Bin_03", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 100.0 && PH_calibrated_RMD < 125.0) AutoHisto("Bin_04", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 125.0 && PH_calibrated_RMD < 150.0) AutoHisto("Bin_05", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 150.0 && PH_calibrated_RMD < 175.0) AutoHisto("Bin_06", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 175.0 && PH_calibrated_RMD < 200.0) AutoHisto("Bin_07", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 200.0 && PH_calibrated_RMD < 225.0) AutoHisto("Bin_08", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 225.0 && PH_calibrated_RMD < 250.0) AutoHisto("Bin_09", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 250.0 && PH_calibrated_RMD < 275.0) AutoHisto("Bin_10", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 275.0 && PH_calibrated_RMD < 300.0) AutoHisto("Bin_11", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 300.0 && PH_calibrated_RMD < 325.0) AutoHisto("Bin_12", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 325.0 && PH_calibrated_RMD < 350.0) AutoHisto("Bin_13", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 350.0 && PH_calibrated_RMD < 375.0) AutoHisto("Bin_14", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 375.0 && PH_calibrated_RMD < 400.0) AutoHisto("Bin_15", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 400.0 && PH_calibrated_RMD < 425.0) AutoHisto("Bin_16", cubicToF, 500,bin_begin,bin_end);
			if (PH_calibrated_RMD > 425.0 && PH_calibrated_RMD < 450.0) AutoHisto("Bin_17", cubicToF, 500,bin_begin,bin_end);

			
			//cout << "RMD Energy is: " << E_rmd<< endl;
			//cout << "LSC Energy is: " << E_lsc<< endl;

			//AutoHisto("RMD Calibrated PH vs RMD_PSD", PH_calibrated_RMD, RMD_PSD, 500,0,1000, 500,-.1,.4);
			
			//AutoHisto("Soft ToF", rmd_softtime - lsq_softtime, 35, -.4,.40);
			//AutoHisto("Cubic ToF", rmd_cubictime - lsq_cubictime,50,-2,2);
			//AutoHisto("RMD4 ToF vs PSD", cubicToF, RMD_PSD, 250,-10,10, 125,-.1,.4);
			AutoHisto(nameforgraph + " ToF Vs Energy", cubicToF, PH_calibrated_RMD, 200,-2.0,4.0,500,0,500);
			//AutoHisto("RMD5 ToF Vs Energy", cubicToF, PH_calibrated_RMD, 250,bin_begin,bin_end,500,0,500);
			
			
		}//end of inner if statement
			
		
	}//end of outer if statement
	
	
	
	
  }//end simple for loop
  
  return kTRUE;
  
  
}//end of main function








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


