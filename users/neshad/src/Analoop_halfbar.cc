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
  TopPulseHeight.resize(NumBars);
  BottomPulseHeight.resize(NumBars);
  AvgPulseHeight.resize(NumBars);

  //Binning of the TOF and energy histograms 
  int EBins=5000; 
  int EMax=100000;
  int TOFBins=2000;
  int TOFMinMax=100;
  int PHBins = 3000;
  int PHMax = 17000;
  



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
     nameStream<<TheSettings->GetBarName(i)<<"_TopPH";
     TopPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PHBins,0,PHMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_BottomPH";
     BottomPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PHBins,0,PHMax);
	 
	 nameStream.str("");
     nameStream<<TheSettings->GetBarName(i)<<"_AvgPH";
     AvgPulseHeight[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),PHBins,0,PHMax);
     
	   }

  // ///////////////////////////////
  // // new more histograms here  //
  ///////////////////////////////
  


  // Add all the histograms to Output list of the Selector
  // fOutput is a TSelectorList which will take ownership of
  // an object when it is added
  fOutput->AddAll(gDirectory->GetList());
}


#define BAD_NUM -10008

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
  int refID = -999;
  double refE = -999;
  long double refCubicTime = -999;
  long double refSoftTime = -999;
  long double TopCubicTime = -999;
  long double TopSoftTime = -999;
  long double TopE = -999;
  long double BotE = -999;
  long double BotCubicTime = -999;
  long double BotSoftTime = -999;
  
  for (int i =0;i<NumBarsInEvent;i++){ //Loop over all Bars in the Lenda event
	
    int NumTopsInBar = lendaevent->Bars[i].NumTops; //Number of tops in ith bar
    int NumBottomsInBar = lendaevent->Bars[i].NumBottoms;//Number of bottoms in ith bar
    int BarId = lendaevent->Bars[i].BarId; //The ith Bars Unique Bar Id
	
	
    for (int t=0;t<NumTopsInBar;t++){ // Loop over all the TOPS in the Bar
      double E=lendaevent->Bars[i].Tops[t].GetEnergy();
	  double ph = lendaevent->Bars[i].Tops[t].GetPulseHeight();
      if (E>0){
	TopPulseHeight[BarId]->Fill(ph);
	TopEnergies[BarId]->Fill(E); //Fill Top Energy if >0
      }
    }//End for over tops
    
    for (int b=0;b<NumBottomsInBar;b++){//Loop over all the Bottoms in the Bar
      double E=lendaevent->Bars[i].Bottoms[b].GetEnergy();
	  double ph = lendaevent->Bars[i].Bottoms[b].GetPulseHeight();
      if (E>0){
	BottomPulseHeight[BarId]->Fill(ph);
	BottomEnergies[BarId]->Fill(E); //Fill Bottom Energy if >0
      }
    }//End for over bottoms

	
    
    if ( NumTopsInBar == 1 ){
      TopTOFs[BarId]->Fill(lendaevent->Bars[i].GetTopTOF());
    }
    if (NumBottomsInBar == 1 ){
      BottomTOFs[BarId]->Fill(lendaevent->Bars[i].GetBottomTOF());
    }
  }//End for over bars in event
  

  bool goodchannel = false;
  long double t1_GE = -999; //variable for raw timing from top PMT of half bar (SL08B)
  long double t2_GE = -999; //variable for raw timing from bottom PMT of half bar (SL08T)
  long double lenda_top = -999; ///lenda top time cubic
  long double lenda_bottom = -999; ///////lenda bottom timme cubic
  double LTopE = -999;
  double LBotE = -999;
  double TopPH = -999;
  double BotPH = -999;
  double t1PH = -999;
  double t2PH = -999;
  double testE = -999;
  double ETopRMD = -999;
  double EBotRMD = -999;
  
  for(int i = 0; i<NumBarsInEvent;i++) /////using this loop for RMD half bar timing.
  {
	 
	 
	 
	  
	bool ref = false;
    string BarName = lendaevent->Bars[i].GetBarName();  
	int BarId = lendaevent->Bars[i].BarId;
	  
	

    //Fill histograms that are only valid when there was 1 top and 1 bottom
	//For LENDA bar coincidences, we want the simple event to be true; both PMTs need to see the event.
    if (lendaevent->Bars[i].SimpleEventBit && NumBarsInEvent ==2)
	{//has 1 top and 1 bottom
		goodchannel = true;
		if(BarName == "NL01") {ref = true; refID = BarId;} //Looking for the reference channel, which is the LENDA bar.
		
		
		if(ref == true) //getting all the values for the ref bar
		{
			 refCubicTime = .5 * ((long double)lendaevent->Bars[i].Tops[0].GetCubicTime() + (long double)lendaevent->Bars[i].Bottoms[0].GetCubicTime()); //ref time from LENDA bar, cubic fit to the timing signal
			 refSoftTime = .5 * ((long double)lendaevent->Bars[i].Tops[0].GetSoftTime() + (long double)lendaevent->Bars[i].Bottoms[0].GetSoftTime()); //ref time from LENDA bar, software time
			 lenda_top = (long double)lendaevent->Bars[i].Tops[0].GetCubicTime();
			 lenda_bottom = (long double)lendaevent->Bars[i].Bottoms[0].GetCubicTime();
			 LTopE = (long double)lendaevent->Bars[i].Tops[0].GetEnergy();
			 LBotE = (long double)lendaevent->Bars[i].Bottoms[0].GetEnergy();
		}
		else //This is the RMD half-bar (if it isn't the LENDA bar)
		  {
			  Int_t start = 1;
			  Double_t rise = -1.;
			  Double_t background = -1.;
			  Double_t noise = -1.;
			  LendaFilter aFilter;
		//LendaChannel rmd_ch = inEvent->Bars[0].Tops[0];
			  LendaChannel top_ch = lendaevent->Bars[i].Tops[0];
			  LendaChannel bot_ch = lendaevent->Bars[i].Bottoms[0];
			  vector <UShort_t> tracetop = top_ch.GetTrace(); //Need the traces to get the energy. Could also add a new function to R00TLe class, might want to do that in the future.
			  vector <UShort_t> tracebot = bot_ch.GetTrace();
			  TopCubicTime = (long double)lendaevent->Bars[i].Tops[0].GetCubicTime(); //Using long doubles because I've had issues previously with rounding
			  BotCubicTime = (long double)lendaevent->Bars[i].Bottoms[0].GetCubicTime(); //uncorrected RMD times
			  
			  
			  ETopRMD = GetEnergyRMD(tracetop, start); //Using a function defined below to get the energy. Sample has a longer tail and integrating more of the signal
			  EBotRMD = GetEnergyRMD(tracebot, start);	  //gives a much better value for the energy.
			 
			 if (ETopRMD>0  && EBotRMD > 0){  //filling the graphs with energies
				TopEnergies[BarId]->Fill(ETopRMD*1.05);  //1.05 is used to scale top energy to the bottom. This was used for center of gravity plot below.
				BottomEnergies[BarId]->Fill(EBotRMD);				//Fill Top Energy if >0
			  }
			
			
			 
			  //TopPH = aFilter.GetPulseComplete(tracetop,start,rise,background,noise);
			  //BotPH = aFilter.GetPulseComplete(tracebot,start,rise,background,noise);
			  TopE = (long double)lendaevent->Bars[i].Tops[0].GetEnergy();
			  BotE = (long double)lendaevent->Bars[i].Bottoms[0].GetEnergy();
			  //TopCubicTime = t1_GE-((-1.86094453E-24) * (TopE*TopE*TopE*TopE*TopE*TopE) + (8.75799325E-20) * (TopE*TopE*TopE*TopE*TopE) + (-1.54835643E-15) * TopE*TopE*TopE*TopE + (1.26676802E-11) *TopE*TopE*TopE + (-4.27633324E-08) *TopE*TopE + (-3.32243163E-05) * TopE + 1.40318873E+00);
			  
				
				
			//if (TopPH == -10008 || BotPH == -10008)  {cout << entry <<" Top: " << TopPH <<" Bot: " << BotPH << endl;}
			  TopSoftTime = lendaevent->Bars[i].Tops[0].GetSoftTime();
			  //BotCubicTime =  t1_GE-((9.848157E-27) * (BotE*BotE*BotE*BotE*BotE*BotE) + (-9.63503164E-22) * (BotE*BotE*BotE*BotE*BotE) + (5.12692129E-17) * BotE*BotE*BotE*BotE + (-1.19696704E-12) *BotE*BotE*BotE + (0.0000000172079657) *BotE*BotE + (-0.000157399065) * BotE + 0.942492944);
			  
			  BotSoftTime = lendaevent->Bars[i].Bottoms[0].GetSoftTime();
			  
			  //cout << "Top Time: " << lendaevent->Bars[i].Tops[0].GetTime() << endl;
			  /*
			  double TdiffCubic = TopCubicTime - BotCubicTime;
			  double TdiffSoft = TopSoftTime - BotSoftTime;
			  */
			  
			  //AutoHisto("Top PH RMD", TopPH, 3000, 0, 17000);
			  //AutoHisto("Bottom PH RMD", BotPH, 3000, 0, 17000);
		  } 
    }
  } //end for loop to fill RMD times
  	if (goodchannel == true)
	{
		long double TtopCub = TopCubicTime - refCubicTime; //timing stuffs, cubic refers to cubic fit of timing, softtime refers to software time (cubic is better)
		long double TbotCub = BotCubicTime - refCubicTime;
		long double TtopSoft = TopSoftTime - refSoftTime;
		long double TbotSoft = BotSoftTime - refSoftTime;
		
		
		double CoG = (1.05 * ETopRMD - EBotRMD)/(1.05 * ETopRMD + EBotRMD); //CoG plot: (aE1 = E2)/(aE1+E2) this isn't a good approximation here. See LENDA bar thesis for more info
		
		t1_GE = TtopCub-(5.34620190e-01 * TMath::Exp( - 6.47921682e-05 * ETopRMD) + 8.34459014e-01);
		t2_GE =  TbotCub-(5.40018088e-01 * TMath::Exp( - 6.54050986e-05 * EBotRMD) + 2.78576985e-01);
				
		//t1PH = TtopCub-((-1.86094453E-24) * (TopPH*TopPH*TopPH*TopPH*TopPH*TopPH) + (8.75799325E-20) * (TopPH*TopPH*TopPH*TopPH*TopPH) + (-1.54835643E-15) * TopPH*TopPH*TopPH*TopPH + (1.26676802E-11) *TopPH*TopPH*TopPH + (-4.27633324E-08) *TopPH*TopPH + (-3.32243163E-05) * TopPH + 1.40318873E+00);
		//t2PH =  TbotCub-((9.848157E-27) * (BotPH*BotPH*BotPH*BotPH*BotPH*BotPH) + (-9.63503164E-22) * (BotPH*BotPH*BotPH*BotPH*BotPH) + (5.12692129E-17) * BotPH*BotPH*BotPH*BotPH + (-1.19696704E-12) *BotPH*BotPH*BotPH + (0.0000000172079657) *BotPH*BotPH + (-0.000157399065) * BotPH + 0.942492944);
		
		
		long double Tcub_avg = (t1_GE+t2_GE)/2;
		long double lendaTDiff = lenda_top-lenda_bottom;
		long double Tavg_exp = (t1_GE + t2_GE)/2;
		double rmd_tdiff = t1_GE - t2_GE;
		double rmd_tdiff_uncorrected = TtopCub - TbotCub;
		
		
		AutoHisto("TopE (RMD) vs uncorrectedT",TtopCub,ETopRMD, 600,-3,3,10000,0,100000);
		AutoHisto("BotE (RMD) vs uncorrectedT",TbotCub,EBotRMD, 600,-3,3,10000,0,100000);
		AutoHisto("TopE vs RMD tdiff",t1_GE-t2_GE, ETopRMD, 600,-3,3,10000,0,100000);
		
		AutoHisto("TopE vs RMDTdiff (uncorrected)", TtopCub - TbotCub, ETopRMD,600,-3,3,10000,0,100000);
		AutoHisto("BotE vs RMDTdiff (uncorrected)", TtopCub - TbotCub,EBotRMD,600,-3,3,10000,0,100000);
		
		AutoHisto("Top E vs LENDATdiff", lendaTDiff, ETopRMD, 600,-3,3,10000,0,100000);
		AutoHisto("Bottom E vs LENDATdiff", lendaTDiff, EBotRMD, 600,-3,3,10000,0,100000);
		
		AutoHisto("BottomE vs RMD tdiff",t1_GE-t2_GE, EBotRMD, 600,-3,3,10000,0,100000);
		AutoHisto("RMD Avg (y) vs RMD Tdiff", t1_GE-t2_GE, Tavg_exp, 500, -2,2,500, -2,2);		
		//AutoHisto("RMD Avg vs NL01 Tdiff", lendaTDiff, Tcub_avg,500, -2,2,500, -2,2);
		AutoHisto("TopE vs RMDAvgt", Tavg_exp, ETopRMD, 600,-3,3,10000,0,100000);
		AutoHisto("BottomE vs RMDAvgt", Tavg_exp, EBotRMD, 600,-3,3,10000,0,100000);
		
		//AutoHisto("NL01T vs RMD Tdiff",t1_GE-t2_GE,LTopE, 600,-3,3,5000,0,100000);
		//AutoHisto("NL01B vs RMD Tdiff",t1_GE-t2_GE,LBotE, 600,-3,3,5000,0,100000);
		
		
		AutoHisto("CoG (y) vs avgT RMD",Tcub_avg, CoG,600,-2,2,1000,-1.5,1.5);
		AutoHisto("CoG (y) vs RMD Tdiff",t1_GE-t2_GE, CoG,600,-2,2,1000,-1.5,1.5);
		AutoHisto("TopE vs avgrefT", refCubicTime,ETopRMD,600,-30,30,10000,0,100000);
		AutoHisto("BottomE vs avgrefT", refCubicTime,EBotRMD,600,-30,30,10000,0,100000);
		
		
		//AutoHisto("Half Bar TopE vs uncorrectedT", TtopCub, TopPH, 1000, -30,30, 10000, 0, 100000);
		//AutoHisto("Half Bar BotE vs uncorrectedT", TbotSoft, BotPH, 600, -30,30, 10000, 0, 100000);
		//AutoHisto("Half Bar TopPH vs T Cubic",TtopCub, TopE, 600,-30,30,10000,0,17000);
		//AutoHisto("Top Time uncorrected",TtopCub, 500, -1,3);
		//AutoHisto("Top E vs T Cubic corrected",t1_GE, ETopRMD, 600,-3,3,10000,0,100000);
		AutoHisto("Bottom E vs T Cubic corrected",t2_GE, EBotRMD, 600,-3,3,10000,0,100000);
		AutoHisto("Top E vs T Corrected (exponential)", t1_GE,ETopRMD, 600,-3,3,10000,0,100000);
		
		//AutoHisto("Half Bar BotE vs T Cubic", TbotCub, BotE, 600, -30,30, 10000, 0, 100000);
		//AutoHisto("Half Bar BotE vs T Software", TbotSoft, BotE, 600, -30,30, 10000, 0, 100000);
		//AutoHisto("Half Bar BotPH vs T Cubic",TbotCub, BotE, 600,-30,30,10000,0,17000);
		//AutoHisto("Bottom Time uncorrected",TbotCub, 500, -1,3);
		//AutoHisto("Half Bar BotPH vs T Cubic",t2_GE, BotE, 600,-30,30,10000,0,17000);
		
		//AutoHisto("Top Corrected time", t1_GE, 500,-2,2);
		//AutoHisto("Top Energy vs Top Corrected Time", t1_GE,TopE, 500, -2,2, 10000,0,17000);
		//AutoHisto("Bottom Corrected time", t2_GE, 500, -2,2);
		//AutoHisto("Bottom Energy vs Bottom Corrected Time", t2_GE,BotE, 500, -3,3, 10000,0,17000);
		/*
		AutoHisto("corrected top time vs Top time - bottom time",t1_GE-t2_GE,t1_GE, 500, -2,2,500, -2,2);
		AutoHisto("corrected bottom time vs top time - bottom time ",t1_GE-t2_GE,t2_GE, 500, -2,2,500, -2,2);
		*/
		AutoHisto("Top Corrected (exp) vs Bottom Corrected", t1_GE,t2_GE, 500, -2,2, 500,-2,2);
		
		if(t1_GE-t2_GE < .1 && t1_GE-t2_GE > -.1) {AutoHisto("RMD AvgT (-.1<tdiff<.1) ",Tcub_avg, 600, -1,1);}
 
		double bin_begin = -.7;
		double bin_end = 2.3;
		/*
		if (ETopRMD > 0.00 && ETopRMD < 3000.00) AutoHisto("Bin_01", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 3000.0 && ETopRMD < 6000.0) AutoHisto("Bin_02", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 6000.0 && ETopRMD < 9000.0) AutoHisto("Bin_03", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 9000.0 && ETopRMD < 12000.0) AutoHisto("Bin_04", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 12000.0 && ETopRMD < 15000.0) AutoHisto("Bin_05", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 15000.0 && ETopRMD < 18000.0) AutoHisto("Bin_06", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 18000.0 && ETopRMD < 21000.0) AutoHisto("Bin_07", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 21000.0 && ETopRMD < 24000.0) AutoHisto("Bin_08", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 24000.0 && ETopRMD < 27000.0) AutoHisto("Bin_09", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 27000.0 && ETopRMD < 30000.0) AutoHisto("Bin_10", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 30000.0 && ETopRMD < 33000.0) AutoHisto("Bin_11", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 33000.0 && ETopRMD < 36000.0) AutoHisto("Bin_12", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 36000.0 && ETopRMD < 39000.0) AutoHisto("Bin_13", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 39000.0 && ETopRMD < 42000.0) AutoHisto("Bin_14", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 42000.0 && ETopRMD < 45000.0) AutoHisto("Bin_15", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 45000.0 && ETopRMD < 48000.0) AutoHisto("Bin_16", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 48000.0 && ETopRMD < 51000.0) AutoHisto("Bin_17", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 51000.0 && ETopRMD < 54000.0) AutoHisto("Bin_18", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 54000.0 && ETopRMD < 57000.0) AutoHisto("Bin_19", TtopCub, 500,bin_begin,bin_end);
		if (ETopRMD > 57000.0 && ETopRMD < 60000.0) AutoHisto("Bin_20", TtopCub, 500,bin_begin,bin_end);
		*/
		/*
		if (EBotRMD > 0.00 && EBotRMD < 3000.00) AutoHisto("Bin_01", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 3000.0 && EBotRMD < 6000.0) AutoHisto("Bin_02", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 6000.0 && EBotRMD < 9000.0) AutoHisto("Bin_03", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 9000.0 && EBotRMD < 12000.0) AutoHisto("Bin_04", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 12000.0 && EBotRMD < 15000.0) AutoHisto("Bin_05", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 15000.0 && EBotRMD < 18000.0) AutoHisto("Bin_06", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 18000.0 && EBotRMD < 21000.0) AutoHisto("Bin_07", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 21000.0 && EBotRMD < 24000.0) AutoHisto("Bin_08", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 24000.0 && EBotRMD < 27000.0) AutoHisto("Bin_09", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 27000.0 && EBotRMD < 30000.0) AutoHisto("Bin_10", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 30000.0 && EBotRMD < 33000.0) AutoHisto("Bin_11", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 33000.0 && EBotRMD < 36000.0) AutoHisto("Bin_12", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 36000.0 && EBotRMD < 39000.0) AutoHisto("Bin_13", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 39000.0 && EBotRMD < 42000.0) AutoHisto("Bin_14", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 42000.0 && EBotRMD < 45000.0) AutoHisto("Bin_15", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 45000.0 && EBotRMD < 48000.0) AutoHisto("Bin_16", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 48000.0 && EBotRMD < 51000.0) AutoHisto("Bin_17", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 51000.0 && EBotRMD < 54000.0) AutoHisto("Bin_18", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 54000.0 && EBotRMD < 57000.0) AutoHisto("Bin_19", TbotCub, 500,bin_begin,bin_end);
		if (EBotRMD > 57000.0 && EBotRMD < 60000.0) AutoHisto("Bin_20", TbotCub, 500,bin_begin,bin_end);
		*/
	}		
  
  
  
  
  return kTRUE;
}




//Below is the energy integration from the LendaFilter class, GetEnergyRMD, with different energy windows. Takes the trace and MaxSpot
Double_t Analoop::GetEnergyRMD(std::vector <UShort_t> &trace,Int_t &MaxSpot){
	
	/////first find the maximum of the pulse
	int _maxSpot=-1;
	Double_t max=0;
	for (int i=0;i<trace.size();i++){
	if (trace[i]>max){
	max=trace[i];
	_maxSpot=i;
	}

	}
	MaxSpot=_maxSpot;
	
	/////now this is the GetEnergyRMD function
	Double_t thisEventsIntegral;
	Double_t sumBegin=0;
	Double_t sumEnd=0;

	Double_t signalIntegral=0;

	int traceLength=trace.size();
	int LengthForBackGround=0.2*traceLength;

	Int_t windowForEnergyleft=4;
	Int_t windowForEnergyright=46;

	if ( MaxSpot - windowForEnergyleft < 0 ){
	return BAD_NUM;
	} else if (MaxSpot + windowForEnergyright > (traceLength -1) ){
	return BAD_NUM;
	}


	for ( int i=0 ;i<LengthForBackGround;i++){
	sumBegin = sumBegin + trace[i];
	sumEnd = sumEnd + trace[traceLength-1-i];
	}

	Double_t BackGround=BAD_NUM;

	if (MaxSpot > LengthForBackGround && MaxSpot < (traceLength-LengthForBackGround) ){
	//The peak is in the middle use the beginning of trace as background
	BackGround=sumBegin/LengthForBackGround;    
	} else if (MaxSpot < LengthForBackGround) {
	//The peak is in the beginning of trace use end as background
	BackGround=sumEnd/LengthForBackGround;
	} else if (MaxSpot > (traceLength-LengthForBackGround)){
	//Peak is at end use beginning for backaground
	BackGround=sumBegin/LengthForBackGround;
	} else{
	//Something makes no sense
	return BAD_NUM;
	}




	for (int i=MaxSpot-windowForEnergyleft;i< MaxSpot+windowForEnergyright;++i) {
	signalIntegral = trace[i]+ signalIntegral;
	}

	if (  signalIntegral - BackGround *(windowForEnergyleft + windowForEnergyright)>0){
	thisEventsIntegral = signalIntegral - BackGround *(windowForEnergyleft + windowForEnergyright);
	}  else{
	thisEventsIntegral = BAD_NUM;
	}
	
	return thisEventsIntegral;
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


