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
  AvgTOFs.resize(NumBars);
  TopTOFs.resize(NumBars);
  BottomTOFs.resize(NumBars);
  BottomEnergies.resize(NumBars);
  TopEnergies.resize(NumBars);
  TopPulseHeight.resize(NumBars);
  Energy_vs_PSD.resize(NumBars);
  ToF_vs_Energy.resize(NumBars);
  BottomPulseHeight.resize(NumBars);
  AvgPulseHeight.resize(NumBars);
  RMD_Calibrated_Energy.resize(NumRMDS);
  RMDPSD.resize(NumRMDS);
  CalibratedEnergy.resize(NumBars);
  PSDs.resize(NumBars);
  PSD_all.resize(NumBars);
  ToF_vs_refEnergy.resize(NumBars);
  ToF_vs_refEnergy_uncorrected.resize(NumBars);
  neutron_KE_30.resize(NumBars);
  neutron_KE_40.resize(NumBars);
  neutron_KE_50.resize(NumBars);
  neutron_KE_60.resize(NumBars);
  neutron_KE_70.resize(NumBars);
  neutron_KE_Am_peak.resize(NumBars);
  KE_vs_E.resize(NumBars);
  
  
  //Binning of the TOF and energy histograms 
  int EBins=30000; 
  int EMax=300000;
  int TOFBins=500;
  int TOFMinMax=50;
  /*
  int PHBins=10000; 
  int PHMax=65000;
  */
  int PHBins = 5000;
  int PHMax = 60000;
  int PHMin = 0;
  int TOFMax = 30;
  int TOFMin = -30;
  int E_cal_min = 0;
  string det_name = "RMD4";
  int EMaxCal = 300;
  int EMinCal = 0;
  int EBinsCal = 1200;
  double PSDMin = -.2;
  double PSDMax = .6;
  int PSDBins = 200;
  string source = "137Cs";
  bool calibrated = false;
  bool Cf_source = false;
  
  
  //For each Bar in the file make these standard histograms
  stringstream nameStream;
  string chname;
  int barchid;

  
  
  for (int i=0;i<NumBars;i++){
	  
		 chname = TheSettings->GetBarName(i);
		 barchid = TheSettings->GetBarId(chname);

		 if (chname == "LSC2")
		 {
			PHBins = 1000;
			PHMax = 20000;
		 }
		 
		 else
		 {
			 PHBins = 5000;
			 PHMax = 63000;
		 }
		 
		 
		 
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
		nameStream<<TheSettings->GetBarName(i)<<"_PSD";
		PSDs[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),200,-.2,.6);


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

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Calibrated_ToF";
		TopTOFs[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),TOFBins,-10,90);

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_PSD_all";
		PSD_all[i] = new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),200,-.2,.6);

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Tdiff_vs_E";
		ToF_vs_Energy[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),500,-10,90,500,0,1500);
		ToF_vs_Energy[i]->GetXaxis()->SetTitle("TDiff (ns)");
		ToF_vs_Energy[i]->GetXaxis()->CenterTitle();
		ToF_vs_Energy[i]->GetYaxis()->SetTitle("Energy (keV)");
		ToF_vs_Energy[i]->GetYaxis()->CenterTitle();
		//ToF_vs_Energy[i]->SetStats(0);

		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Tdiff_vs_refE";
		ToF_vs_refEnergy[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),500,-50,50,500,0,2000);
		ToF_vs_refEnergy[i]->GetXaxis()->SetTitle("TDiff (ns)");
		ToF_vs_refEnergy[i]->GetXaxis()->CenterTitle();
		ToF_vs_refEnergy[i]->GetYaxis()->SetTitle("Energy (keV)");
		ToF_vs_refEnergy[i]->GetYaxis()->CenterTitle();
		//ToF_vs_Energy[i]->SetStats(0);
		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_Tdiff_uncorrected_vs_refE";
		ToF_vs_refEnergy_uncorrected[i]= new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),500,0,100,500,0,2000);
		ToF_vs_refEnergy_uncorrected[i]->GetXaxis()->SetTitle("TDiff (ns)");
		ToF_vs_refEnergy_uncorrected[i]->GetXaxis()->CenterTitle();
		ToF_vs_refEnergy_uncorrected[i]->GetYaxis()->SetTitle("Energy (keV)");
		ToF_vs_refEnergy_uncorrected[i]->GetYaxis()->CenterTitle();

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_E_vs_PSD";
		Energy_vs_PSD[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),350,0,1400,200,-.1,.3);
		Energy_vs_PSD[i]->GetXaxis()->SetTitle("Energy (keV)");
		Energy_vs_PSD[i]->GetXaxis()->CenterTitle();
		Energy_vs_PSD[i]->GetYaxis()->SetTitle("PSD");
		Energy_vs_PSD[i]->GetYaxis()->CenterTitle();

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_neutron_KE_30kevee";
		neutron_KE_30[i]= new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),40,0,10000);
		neutron_KE_30[i]->GetXaxis()->SetTitle("Energy (keV)");
		neutron_KE_30[i]->GetXaxis()->CenterTitle();

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_neutron_KE_40kevee";
		neutron_KE_40[i]= new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),40,0,10000);
		neutron_KE_40[i]->GetXaxis()->SetTitle("Energy (keV)");
		neutron_KE_40[i]->GetXaxis()->CenterTitle();

		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_neutron_KE_50kevee";
		neutron_KE_50[i]= new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),40,0,10000);
		neutron_KE_50[i]->GetXaxis()->SetTitle("Energy (keV)");
		neutron_KE_50[i]->GetXaxis()->CenterTitle();
		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_neutron_KE_60kevee";
		neutron_KE_60[i]= new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),40,0,10000);
		neutron_KE_60[i]->GetXaxis()->SetTitle("Energy (keV)");
		neutron_KE_60[i]->GetXaxis()->CenterTitle();
		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_neutron_KE_Am_peakkevee";
		neutron_KE_Am_peak[i]= new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),40,0,10000);
		neutron_KE_Am_peak[i]->GetXaxis()->SetTitle("Energy (keV)");
		neutron_KE_Am_peak[i]->GetXaxis()->CenterTitle();
		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_neutron_KE_70kevee";
		neutron_KE_70[i]= new TH1F(nameStream.str().c_str(),nameStream.str().c_str(),40,0,10000);
		neutron_KE_70[i]->GetXaxis()->SetTitle("Energy (keV)");
		neutron_KE_70[i]->GetXaxis()->CenterTitle();
		
		
		nameStream.str("");
		nameStream<<TheSettings->GetBarName(i)<<"_KE_vs_E";
		KE_vs_E[i] = new TH2F(nameStream.str().c_str(),nameStream.str().c_str(),200,0,4000,400,0,1600); //////was using 500,0,10000,500,0,2000
		KE_vs_E[i]->GetXaxis()->SetTitle("Tn (keV)");
		KE_vs_E[i]->GetXaxis()->CenterTitle();
		KE_vs_E[i]->GetYaxis()->SetTitle("Light output (keVee)");
		KE_vs_E[i]->GetYaxis()->CenterTitle();
		
  }

  // ///////////////////////////////
  // // new more histograms here  //
  ///////////////////////////////
  


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
  int refLSC = 999;
  int RMD = 999;
  bool ref_single = false;
  bool ref_bar = false;
  int Ref_single_BarId = 999;
  double refE_single = 999;
  string LSCname = "LSC1";
  string RMDname = "RMD6";
  double softToF = 999;
  double cubicToF = 999;
  double ref_single_softtime = 999;
  double ref_single_cubictime = 999;
  double soft_tof = 999;
  double cubic_tof = 999;
  double rmd_softtime = 999;
  double rmd_cubictime = 999;
  double ref_single_softime_cutoff = 999;
  double ref_single_cubictime_cutoff = 999;
  double lsq_softtime = 999;
  double lsq_cubictime = 999;
  double slope = 1;
  double intercept = 0;
  string LENDAname = "NL01";
  double rmd4_scale_factor = 0.017905004;
  string RMD1name = "RMD1";
  string RMD2name = "RMD2";
  string RMD3name = "RMD3";
  string RMD4name = "RMD4";
  string RMD5name = "RMD5";
  string RMD6name = "RMD6";
  string LSC1name = "LSC1";
  string LSC2name = "LSC2";
  string refdet = "LSC2";
  double refPSD = -999;
  double PSD = -999;
  double reftime = -999.0000;
  double tof = -999;
  

	TString runnum = "5103";
  
  for (int i = 0; i<NumBarsInEvent; i++) ///This is the loop for neutron efficiency calculations. 
  {
	  
	string BarName = lendaevent->Bars[i].GetBarName(); //bar name
    int NumTopsInBar = lendaevent->Bars[i].NumTops; //Number of tops in ith bar
    int NumBottomsInBar = lendaevent->Bars[i].NumBottoms;//Number of bottoms in ith bar
    int BarId = lendaevent->Bars[i].BarId; //The ith Bars Unique Bar Id	  
	
	
	if (BarName == RMD1name)    {slope = .0230337; intercept = -71.69704;} //Energy calibrations for the bars/samples. Mostly doing this here to plot calibrated energies
	if (BarName == RMD4name)	{slope = .023552; intercept = -61.9719;}
	if (BarName == RMD5name) 	{slope = .0553628; intercept = -154.29709;}
	if (BarName == RMD6name) 	{slope = .0535289; intercept = -156.9073;}
	if (BarName == LSC1name)	{slope = .021405; intercept = -56.53591;}		//{slope = .02474467; intercept = -103.72709;}
	if (BarName == LSC2name)	{slope = .072599; intercept = -38.97435;}		//{slope = .088313; intercept = -99.8487;}
	 
	 
	 
	 
	double calibrated_E = lendaevent->Bars[i].Tops[0].GetPulseHeight() * slope + intercept; //This gives calibrated energy in keVee (light output calibrated using 
	CalibratedEnergy[BarId]->Fill(calibrated_E);											//241Am (high energy photopeak), 22Na (both compton edges), 137Cs (compton edge)
	
	
	
	
	if (NumBarsInEvent >=2 && BarName == refdet) //gating this an easy way, gating on ref detector only if there are 2 or more bars in event
		{
			
			refPSD = (1-lendaevent->Bars[i].Tops[0].GetShortGate()/lendaevent->Bars[i].Tops[0].GetLongGate()); //calculate ref PSD
			
			
			long double reftime = (long double)lendaevent->Bars[i].Tops[0].GetCorseTime() + (long double)lendaevent->Bars[i].Tops[0].GetCubicCFD(); //for some reason this needs to be a long double
																																					//Had significant rounding issues w/ normal double
			long double refE = (long double)lendaevent->Bars[i].Tops[0].GetPulseHeight()*.072599-38.97435;							  //calculating reference energy to use later
			long double cutoff_psd = refE * -.0000574879 + .072246;																	//This PSD cut is energy dependent due to poor PSD resolution in ref detector
			long double correction_5101 = (5.15539398e-13* refE*refE*refE*refE -5.92791292e-10 * refE*refE*refE + -8.94786772e-07 * refE*refE + 9.85034725e-04 * refE + 2.40147449e+01);//ToF vs E corrections for ref bar run 5101 (RMD4 + LSC1 +ref)
			long double correction_5103 = (-2.03696789e-12* refE*refE*refE*refE +6.08906433e-09 * refE*refE*refE + -7.19884255e-06* refE*refE + 3.53678719e-03 * refE + 2.36492113e+01);//ToF vs E corrections for ref bar run 5103 (RMD4 + LSC1 +ref)
			
			long double correction = (5.15539398e-13* refE*refE*refE*refE -5.92791292e-10 * refE*refE*refE + -8.94786772e-07 * refE*refE + 9.85034725e-04 * refE + 2.40147449e+01); //ToF correction for ref LSC detector
																																													//ToF had energy dependencies
			
			
			if (refPSD < cutoff_psd && refE > 300) //This is PSD/Energy cut. Gating on gamma peak in the ref detector
			{
				//AutoHisto("LSC2 E VS PSD Cut", refE, refPSD,350,0,1400,200,-.1,.3); 

				for (int j = 0; j<NumBarsInEvent; j++) //Now that I've made the cut on ref detector, need to get the events in other detectors
				  {
					 string TempBarName = lendaevent->Bars[j].GetBarName();
					 int TempBarId = lendaevent->Bars[j].BarId;
					 if (TempBarName == RMD1name)    {slope = .0230337; intercept = -71.69704;}
					 if (TempBarName == RMD4name)	{slope = .023552; intercept = -61.9719;}
					 if (TempBarName == RMD5name) 	{slope = .0553628; intercept = -154.29709;}
					 if (TempBarName == RMD6name) 	{slope = .0535289; intercept = -156.9073;}
					 if (TempBarName == LSC1name)	{slope = .021405; intercept = -56.53591;}			/////  {slope = .02474467; intercept = -103.72709;}
					 if (TempBarName == LSC2name)	{slope = .072599; intercept = -38.97435;}     			/////{slope = .088313; intercept = -99.8487;} 

					 long double E= (long double)lendaevent->Bars[j].Tops[0].GetPulseHeight()*(long double)slope + (long double)intercept;
					 long double tof = 4*((long double)lendaevent->Bars[j].Tops[0].GetCorseTime() + (long double)lendaevent->Bars[j].Tops[0].GetCubicCFD() - reftime);
					 long double tof_ref_adjusted = -999.;
					 long double tof_corrected = -999.;
					 long double c = .299702547;
					 long double L;
					 long double KE = -999;
					 long double RMD1_correction = (3.36876939e-12*E*E*E*E -8.93186479e-09*E*E*E + 7.84291107e-06*E*E -2.60376056e-03*E + 2.33896104e+01) - 1.165322; //ToF correction; detector response is energy dependent
					
					
					
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////Timing Corrections //////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
					
					if(TempBarName == LSC1name) //Set the tof correction and the distance L from the source to middle of detection medium
					{
						if(runnum == "5101")
						{
							tof_corrected = tof - (1.13238958e-12*E*E*E*E - 4.41258169e-09*E*E*E + 5.96533986e-06*E*E -  2.94342706e-03*E + 2.45124661e+01) + 1.1971036;
							tof_ref_adjusted = tof - correction_5101 ;
						}
						
						else if(runnum == "5103")
						{
							tof_corrected = tof - (1.41114867e-12*E*E*E*E -4.18811263e-09*E*E*E + 4.87179649e-06*E*E -2.29122834e-03*E + 2.44275961e+01) + 1.1971036;
							tof_ref_adjusted = tof - correction_5103 ;
						}
						
						else
						{
							tof_corrected = tof;
							tof_ref_adjusted = tof;
						}
						L = .358775+.0254;
					}//end LSC1 corrections
					
					
					if(TempBarName == RMD4name) //set tof for RMD4
					{
						tof_corrected = tof - 23.54799454 + 1.165322;//(2.74815094e-12* E*E*E*E - 8.14037687e-09* E*E*E + 7.65274543e-06 * E*E - 2.59176933e-03* E + 2.37493074e+01);
						L = .34925 + .0254/2;
					}
					
					
					
					if(TempBarName == RMD1name) //set tof for RMD1
					{
						tof_corrected = tof - RMD1_correction;
						L = .34925 + .0254/2;
					}
					
					
					
					if (tof_corrected > -100) //cut out some nonsensical events
					{
						ToF_vs_Energy[TempBarId]->Fill(tof_corrected,E);
						TopTOFs[TempBarId]->Fill(tof_corrected);
					}
					
					
					/*
					 int bin_begin = -1;
					 int bin_end = 5;
					 int bins = 400;
					 int ref_bin_begin = -4;
					 int ref_bin_end = 4;
					 			
					 
					 if (TempBarName == LSC1name && tof_corrected >bin_begin  && tof_corrected < bin_end) //plotting timing bins, this was used to determine ToF vs E corrections
					 {
						 
						if (E > 0.00 && E <100) AutoHisto("Bin_01", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 100  && E <200) AutoHisto("Bin_02", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 200. && E <300) AutoHisto("Bin_03", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 300  && E <400) AutoHisto("Bin_04", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 400  && E <500) AutoHisto("Bin_05", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 500  && E <600) AutoHisto("Bin_06", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 600. && E <700) AutoHisto("Bin_07", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 700  && E <800) AutoHisto("Bin_08", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 800. && E <900) AutoHisto("Bin_09", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 900  && E <1000) AutoHisto("Bin_10", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 1000  && E <1100) AutoHisto("Bin_11", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 1100  && E <1200) AutoHisto("Bin_12", tof_corrected, bins/2,bin_begin,bin_end);
						if (E > 1200  && E <1300) AutoHisto("Bin_13", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 1300  && E <1345) AutoHisto("Bin_14", tof_corrected, bins/2,bin_begin,bin_end);
					 }
					 
					 
					 if (TempBarName == LSC1name && tof_ref_adjusted >ref_bin_begin  && tof_ref_adjusted < ref_bin_end) //
					 {
						 
						//cout << "LSC1 time is: " << lendaevent->Bars[j].Tops[0].GetCubicTime() << endl;
			 			if (refE > 300  && refE <400) AutoHisto("ref_Bin_01", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 400  && refE <500) AutoHisto("ref_Bin_02", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 500  && refE <600) AutoHisto("ref_Bin_03", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 600. && refE <700) AutoHisto("ref_Bin_04", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 700  && refE <800) AutoHisto("ref_Bin_05", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 800. && refE <900) AutoHisto("ref_Bin_06", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			if (refE > 900  && refE <1000) AutoHisto("ref_Bin_07", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			if (refE > 1000  && refE <1100) AutoHisto("ref_Bin_08", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			if (refE > 1100  && refE <1147) AutoHisto("ref_Bin_09", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
						//if (refE > 1200  && refE <1300) AutoHisto("ref_Bin_10", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			//if (refE > 1300  && refE <1345) AutoHisto("ref_Bin_11", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
					 }
					 
					 if (TempBarName == RMD1name && tof_ref_adjusted >ref_bin_begin  && tof_ref_adjusted < ref_bin_end) //
					 {
						 
						//cout << "LSC1 time is: " << lendaevent->Bars[j].Tops[0].GetCubicTime() << endl;
			 			if (refE > 300  && refE <400) AutoHisto("RMD1_ref_Bin_01", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 400  && refE <500) AutoHisto("RMD1_ref_Bin_02", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 500  && refE <600) AutoHisto("RMD1_ref_Bin_03", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 600. && refE <700) AutoHisto("RMD1_ref_Bin_04", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 700  && refE <800) AutoHisto("RMD1_ref_Bin_05", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 800. && refE <900) AutoHisto("RMD1_ref_Bin_06", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			if (refE > 900  && refE <1000) AutoHisto("RMD1_ref_Bin_07", tof_ref_adjusted, bins/4,ref_bin_begin,ref_bin_end);
			 			if (refE > 1000  && refE <1100) AutoHisto("RMD1_ref_Bin_08", tof_ref_adjusted, bins/4,ref_bin_begin,ref_bin_end);
			 			if (refE > 1100  && refE <1147) AutoHisto("RMD1_ref_Bin_09", tof_ref_adjusted, bins/4,ref_bin_begin,ref_bin_end);
						//if (refE > 1200  && refE <1300) AutoHisto("RMD_ref_Bin_10", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			//if (refE > 1300  && refE <1345) AutoHisto("RMD_ref_Bin_11", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
					 }
					 
					 
					 
					 if (TempBarName == RMD1name && tof_corrected >bin_begin  && tof_corrected < bin_end) //
					 {
						 bins = 400;
						if (E > 0.00 && E <100) AutoHisto("RMD1_Bin_01", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 100  && E <200) AutoHisto("RMD1_Bin_02", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 200. && E <300) AutoHisto("RMD1_Bin_03", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 300  && E <400) AutoHisto("RMD1_Bin_04", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 400  && E <500) AutoHisto("RMD1_Bin_05", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 500  && E <600) AutoHisto("RMD1_Bin_06", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 600. && E <700) AutoHisto("RMD1_Bin_07", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 700  && E <800) AutoHisto("RMD1_Bin_08", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 800. && E <900) AutoHisto("RMD1_Bin_09", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 900  && E <1000) AutoHisto("RMD1_Bin_10", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 1000  && E <1100) AutoHisto("RMD1_Bin_11", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 1100  && E <1400) AutoHisto("RMD1_Bin_12", tof_corrected, bins/4,bin_begin,bin_end);
					 }
					
					 
					 if (TempBarName == RMD4name && tof_ref_adjusted >ref_bin_begin  && tof_ref_adjusted < ref_bin_end) //
					 {
						 
						//cout << "LSC1 time is: " << lendaevent->Bars[j].Tops[0].GetCubicTime() << endl;
			 			if (refE > 300  && refE <400) AutoHisto("RMD_ref_Bin_01", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 400  && refE <500) AutoHisto("RMD_ref_Bin_02", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 500  && refE <600) AutoHisto("RMD_ref_Bin_03", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 600. && refE <700) AutoHisto("RMD_ref_Bin_04", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 700  && refE <800) AutoHisto("RMD_ref_Bin_05", tof_ref_adjusted, bins,ref_bin_begin,ref_bin_end);
			 			if (refE > 800. && refE <900) AutoHisto("RMD_ref_Bin_06", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			if (refE > 900  && refE <1000) AutoHisto("RMD_ref_Bin_07", tof_ref_adjusted, bins/4,ref_bin_begin,ref_bin_end);
			 			if (refE > 1000  && refE <1100) AutoHisto("RMD_ref_Bin_08", tof_ref_adjusted, bins/4,ref_bin_begin,ref_bin_end);
			 			if (refE > 1100  && refE <1147) AutoHisto("RMD_ref_Bin_09", tof_ref_adjusted, bins/4,ref_bin_begin,ref_bin_end);
						//if (refE > 1200  && refE <1300) AutoHisto("RMD_ref_Bin_10", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
			 			//if (refE > 1300  && refE <1345) AutoHisto("RMD_ref_Bin_11", tof_ref_adjusted, bins/2,ref_bin_begin,ref_bin_end);
					 }
					 
					 
					 if (TempBarName == RMD4name && tof_corrected >bin_begin  && tof_corrected < bin_end) //
					 {
						 bins = 200;
						//cout << "LSC1 time is: " << lendaevent->Bars[j].Tops[0].GetCubicTime() << endl;
						if (E > 0.00 && E <100) AutoHisto("RMD_Bin_01", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 100  && E <200) AutoHisto("RMD_Bin_02", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 200. && E <300) AutoHisto("RMD_Bin_03", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 300  && E <400) AutoHisto("RMD_Bin_04", tof_corrected, bins,bin_begin,bin_end);
			 			if (E > 400  && E <500) AutoHisto("RMD_Bin_05", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 500  && E <600) AutoHisto("RMD_Bin_06", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 600. && E <700) AutoHisto("RMD_Bin_07", tof_corrected, bins/2,bin_begin,bin_end);
			 			if (E > 700  && E <800) AutoHisto("RMD_Bin_08", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 800. && E <900) AutoHisto("RMD_Bin_09", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 900  && E <1000) AutoHisto("RMD_Bin_10", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 1000  && E <1100) AutoHisto("RMD_Bin_11", tof_corrected, bins/4,bin_begin,bin_end);
			 			if (E > 1100  && E <1200) AutoHisto("RMD_Bin_12", tof_corrected, bins/4,bin_begin,bin_end);
						if (E > 1200  && E <1470) AutoHisto("RMD_Bin_13", tof_corrected, bins/4,bin_begin,bin_end);
					 }
					
					
					
					*/
					
					
					
					
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					/////////////////////////////////////Neutron Efficiency Calculations /////////////////////////////////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
					
					
					if(tof_corrected > 10.0 && TempBarName != LSC2name) //cutting off gammas to look at just the neutrons
					{
						long double gamma_sqrd = TMath::Sqrt(1/(1-(L/(c*tof_corrected))*(L/(c*tof_corrected))));
						KE = 939.565 * (gamma_sqrd-1)*1000;
						long double E_cut_lsc1 = 7.45309853e-05 * KE * KE + 4.47899886e-02 * KE + 4.45891414e+01;  //Background cuts; light output is non-linear for neutron energies <~4 MeV or so.
						long double E_cut_rmd4 = 9.01652625e-05 * KE * KE + 7.34638467e-02 * KE + 3.42166658e+01;  //Therefore, cutting events that lie above these cuts
						long double E_cut_rmd1 = 6.04309726e-05 * KE * KE + 1.57759432e-01 * KE + 1.33836688e+01;
						
						if (TempBarName == LSC1name && E < E_cut_lsc1) //Instituting different keVee threholds for efficiency comparisons
						{
							
							if(E > 30 && E < 1340)
							{
								neutron_KE_30[TempBarId]->Fill(KE);
							}
							
							if(E > 40 && E < 1340)
							{
								neutron_KE_40[TempBarId]->Fill(KE);
							}
							
							if(E > 50 && E < 1340)
							{
								neutron_KE_50[TempBarId]->Fill(KE);
							}
							
							if(E > 60.00 && E < 1340)   /////peak is at 59.95
							{
								neutron_KE_60[TempBarId]->Fill(KE);
							}
							
							if(E > 59.95 && E < 1340)   /////Threshold at the Am peak
							{
								neutron_KE_Am_peak[TempBarId]->Fill(KE);
								//KE_vs_E[TempBarId]->Fill(KE, E);
							}
							
							if(E > 70 && E < 1340)
							{
								neutron_KE_70[TempBarId]->Fill(KE);
							}
							KE_vs_E[TempBarId]->Fill(KE, E);
							
						}//end LSC1
						
						else if (TempBarName == RMD4name && E < E_cut_rmd4) //Same as LSC1 but with RMD4
						{
							if(E>30 && E < 1470) //Note that the high energy cut is dependent on the detector; it is where the pixie cards were oversaturated.
							{						//This was used as the lowest keVee cut/highest efficiency case
								neutron_KE_30[TempBarId]->Fill(KE);
							}
							
							if(E>40 && E < 1470)
							{
								neutron_KE_40[TempBarId]->Fill(KE);
							}
							
							if(E>50 && E < 1470)
							{
								neutron_KE_50[TempBarId]->Fill(KE);
							}
							
							if(E>60 && E <1470)   //////61.9 is where the Am peak really sits here 
							{
								neutron_KE_60[TempBarId]->Fill(KE);
							}
							
							if(E > 61.9 && E < 1470)   /////At the Am peak
							{
								neutron_KE_Am_peak[TempBarId]->Fill(KE);
								//KE_vs_E[TempBarId]->Fill(KE, E);
							}
							
							if(E > 70 && E < 1470)
							{
								neutron_KE_70[TempBarId]->Fill(KE);
							}
							KE_vs_E[TempBarId]->Fill(KE, E);
						}//end RMD4
						
						else if (TempBarName == RMD1name && E < E_cut_rmd1) //same cut but with RMD1
						{
							if(E>30 && E < 1430) //This was used as the lowest keVee cut/highest efficiency case
							{
								neutron_KE_30[TempBarId]->Fill(KE);
							}
							
							if(E>40 && E < 1430)
							{
								neutron_KE_40[TempBarId]->Fill(KE);
							}
							
							if(E>50 && E < 1430)
							{
								neutron_KE_50[TempBarId]->Fill(KE);
							}
							
							if(E>60.00 && E <1430) ////////Am peak sits at 52.7645
							{
								neutron_KE_60[TempBarId]->Fill(KE);
							}
							
							if(E > 52.7645 && E < 1470)   /////At the Am peak, 
							{
								neutron_KE_Am_peak[TempBarId]->Fill(KE);
								//KE_vs_E[TempBarId]->Fill(KE, E);
							}
							
							if(E > 70 && E < 1430)
							{
								neutron_KE_70[TempBarId]->Fill(KE);
							}
							KE_vs_E[TempBarId]->Fill(KE, E);
						}//end RMD1
						
						
						//KE_vs_E[TempBarId]->Fill(KE, E);
					}//end ToF cut if statement
				
					 
				  }	//end inner loop over bars in event
				  
			}//end PSD cut on ref detector
			
		}  //end if statement checking for ref detector
	  
	  
	  
  } //

  
  
  
  
  
  /*
  for (int i = 0; i<NumBarsInEvent;i++) //this is set up to do timing for
  {                                    //the LSD/RMD detectors. Will now play with energy cuts      

		  
		   
    int NumTopsInBar = lendaevent->Bars[i].NumTops; //Number of tops in ith bar
    int NumBottomsInBar = lendaevent->Bars[i].NumBottoms;//Number of bottoms in ith bar
    int BarId = lendaevent->Bars[i].BarId; //The ith Bars Unique Bar Id
    string temprefname = lendaevent->Bars[i].GetBarName();
	double E = lendaevent->Bars[i].Tops[0].GetEnergy();
	double E_rmd = -999;
	double E_lsc = -999;
	double lG = -999;
	double sG = -999;
	double lg_ref = lendaevent->Bars[i].Tops[0].GetLongGate();
	double sg_ref = lendaevent->Bars[i].Tops[0].GetShortGate();
	double refPSD = (lg_ref-sg_ref)/sg_ref;
	double PH = -999; 
	double PH_calibrated_RMD = -999;
	double PH_lsc = -999;
	double RMD_PSD = -999;
	double PH_calibrated_LSC = -999;
	double bin_begin = 5.5;
	double bin_end = 11.5;
	double lenda_tof = -999;
	string detector = "RMD4";
	bool run = true;
	bool calibrated = true;
	//cout << " " << endl;
	
	
	
	if (temprefname == RMD6name){
	for (int t=0;t<NumTopsInBar;t++)
	{
		double E=lendaevent->Bars[i].Tops[t].GetEnergy();
		double PH=lendaevent->Bars[i].Tops[t].GetPulseHeight();
		if (PH>0){
			AutoHisto("RMD6 Am241", PH, 700, 3000, 10000);
			}
	}
	}


	
	
	
	
	if (temprefname == RMD4name) 
	{

		//RMD_Calibrated_Energy[0]->Fill(rmd_slope * lendaevent->Bars[i].Tops[0].GetPulseHeight() + rmd_intercept);
		for (int t=0;t<NumTopsInBar;t++)
		{ // Loop over all the TOPS in the Bar
			double E=lendaevent->Bars[i].Tops[t].GetEnergy();
			double PH = lendaevent->Bars[i].Tops[0].GetPulseHeight();
			double E_calibrated = PH * rmd_slope + rmd_intercept;
			double PSD = (1-lendaevent->Bars[i].Tops[0].GetShortGate()/lendaevent->Bars[i].Tops[0].GetLongGate());
			AutoHisto("PSD All",  PSD, 400, -.2,.6);
			Energy_vs_PSD[BarId]->Fill(E_calibrated, PSD);
			
			
			if (E_calibrated > 200 && E_calibrated < 1450 )
			{
				AutoHisto("RMD4 Ratio Above 200 keV", PSD, 400, -.2,.6);
				AutoHisto("RMD4 energy vs PSD >200 keV", E_calibrated, PSD, 2000,0,4000,200,-.2,.6);
				RMDPSD[0]->Fill(PSD);
			}
			
			if (E_calibrated > 500 && E_calibrated < 1450 )
			{
				AutoHisto("RMD4 all above 500 keV", PSD, 400, -.2,.6);
				AutoHisto("RMD4 energy vs PSD >500 keV", E_calibrated, PSD, 2000,0,4000,200,-.2,.6);
			}
			
			Energy_vs_PSD[BarId]->Fill(E_calibrated, PSD);
			if (PSD > .131){
				if (E_calibrated > 0 && E_calibrated < 1450)
				{
					AutoHisto("RMD4_above_threshold_neutrons", PSD, 400, -.2,.6);
				}
				if(E_calibrated > 200 && E_calibrated < 1450)
				{
					AutoHisto("RMD4_above_200_neutrons", PSD, 400, -.2,.6);
				}
				if(E_calibrated > 500 && E_calibrated < 1450)
				{
					AutoHisto("RMD4_above_500_neutrons", PSD, 400, -.2,.6);
				}
				//AutoHisto("RMD5_PSD", PSD, 800, -.2,.6);
				//TopEnergies[BarId]->Fill(E_calibrated); //Fill Top Energy if >0
				//TopPulseHeight[BarId]->Fill(E_calibrated);
				
				
					}
					
			else
			{
				if (E_calibrated > 0 && E_calibrated < 1450)
				{
					AutoHisto("RMD4_above_threshold_gammas", PSD, 400, -.2,.6);
				}
				if(E_calibrated > 200 && E_calibrated < 1450)
				{
					AutoHisto("RMD4_above_200_gammas", PSD, 400, -.2,.6);
				}
				if(E_calibrated > 500 && E_calibrated < 1450)
				{
					AutoHisto("RMD4_above_500_gammas", PSD, 400, -.2,.6);
				}
				
			}
			
		}//End for over tops
		
		
		
	} //End of if statement
	
	
	
		
	if (lendaevent->NumBars == 2 && temprefname == LENDAname) //This is for the timing measurements
	{
		double topE = lendaevent->Bars[i].Tops[0].GetPulseHeight();
		double botE = lendaevent->Bars[i].Bottoms[0].GetPulseHeight();
		double E_calibrated = -999;
		
		
		if (topE > 3000 && botE > 4000) /////making an energy cut, don't have top and bottom calibrated but this should be fine using 3k and 4k respectively
		{
			
			for(int j = 0; j < NumBarsInEvent;j++)
			{
				double E = lendaevent->Bars[j].Tops[0].GetEnergy();
				string tempname2 = lendaevent->Bars[j].GetBarName();
				if (tempname2 == LENDAname)
				{
				   topE = lendaevent->Bars[j].Tops[0].GetPulseHeight();
				   botE = lendaevent->Bars[j].Bottoms[0].GetPulseHeight();
				   
				   lenda_tof = .5 * (lendaevent->Bars[j].Tops[0].GetCubicTime() +
				   lendaevent->Bars[j].Bottoms[0].GetCubicTime());
				
				//	sg_ref = lendaevent->Bars[j].Tops[0].GetShortGate();
				//	lg_ref = lendaevent->Bars[j].Tops[0].GetLongGate();
				//	refPSD = (lg_ref-sg_ref)/sg_ref;
				//	E_lsc = E;
					//PH_lsc = lendaevent->Bars[j].Tops[0].GetPulseHeight();
					//PH_calibrated_LSC = 0.040665 * PH_lsc - 134.759;
				//	lsq_softtime = lendaevent->Bars[j].Tops[0].GetSoftTime();
				//	lsq_cubictime = lendaevent->Bars[j].Tops[0].GetCubicTime();
				} //end inner-inner-inner if
				else
				{
					PH = lendaevent->Bars[j].Tops[0].GetPulseHeight();
					PH_calibrated_RMD = PH * rmd_slope + rmd_intercept;
					lG = lendaevent->Bars[j].Tops[0].GetLongGate();
					sG = lendaevent->Bars[j].Tops[0].GetShortGate();
					RMD_PSD = (lG-sG)/lG;
					E_rmd = E;
					rmd_softtime = lendaevent->Bars[j].Tops[0].GetSoftTime();
					rmd_cubictime = lendaevent->Bars[j].Tops[0].GetCubicTime();
				}//end else
				
				
				
			}//end inner for

			cubicToF = (rmd_cubictime - lenda_tof);
			
			
			if (PH_calibrated_RMD > 0.00 && PH_calibrated_RMD < 50.0) AutoHisto("Bin_01", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 50.0 && PH_calibrated_RMD < 75.0) AutoHisto("Bin_02", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 75.0 && PH_calibrated_RMD < 100.0) AutoHisto("Bin_03", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 100.0 && PH_calibrated_RMD < 125.0) AutoHisto("Bin_04", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 125.0 && PH_calibrated_RMD < 150.0) AutoHisto("Bin_05", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 150.0 && PH_calibrated_RMD < 175.0) AutoHisto("Bin_06", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 175.0 && PH_calibrated_RMD < 200.0) AutoHisto("Bin_07", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 200.0 && PH_calibrated_RMD < 225.0) AutoHisto("Bin_08", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 225.0 && PH_calibrated_RMD < 250.0) AutoHisto("Bin_09", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 250.0 && PH_calibrated_RMD < 275.0) AutoHisto("Bin_10", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 275.0 && PH_calibrated_RMD < 300.0) AutoHisto("Bin_11", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 300.0 && PH_calibrated_RMD < 325.0) AutoHisto("Bin_12", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 325.0 && PH_calibrated_RMD < 350.0) AutoHisto("Bin_13", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 350.0 && PH_calibrated_RMD < 375.0) AutoHisto("Bin_14", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 375.0 && PH_calibrated_RMD < 400.0) AutoHisto("Bin_15", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 400.0 && PH_calibrated_RMD < 425.0) AutoHisto("Bin_16", cubicToF, 250,bin_begin,bin_end);
			if (PH_calibrated_RMD > 425.0 && PH_calibrated_RMD < 450.0) AutoHisto("Bin_17", cubicToF, 250,bin_begin,bin_end);

			
			//cout << "RMD Energy is: " << E_rmd<< endl;
			//cout << "LSC Energy is: " << E_lsc<< endl;

			//AutoHisto("RMD Calibrated PH vs RMD_PSD", PH_calibrated_RMD, RMD_PSD, 500,0,1000, 500,-.1,.4);
			
			//AutoHisto("Soft ToF", rmd_softtime - lsq_softtime, 35, -.4,.40);
			//AutoHisto("Cubic ToF", rmd_cubictime - lsq_cubictime,50,-2,2);
			//AutoHisto("RMD4 ToF vs PSD", cubicToF, RMD_PSD, 250,-10,10, 125,-.1,.4);
			//AutoHisto("RMD3 ToF Vs Energy", cubicToF, PH_calibrated_RMD, 400,5,12,500,0,500);
			//AutoHisto("RMD5 ToF Vs Energy", cubicToF, PH_calibrated_RMD, 250,bin_begin,bin_end,500,0,500);
			ToF_vs_Energy[0]->Fill(cubicToF,PH_calibrated_RMD);
			
			AutoHisto("NL01T ToF vs Energy", cubicToF, topE, 400, 5, 12, 500,0,20000);
			AutoHisto("NL01B ToF vs Energy", cubicToF, botE, 400, 5, 12, 500,0,20000);
		}//end of inner if statement
			
		
	}//end of outer if statement
	
		 
  }//End of for loop
  */
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


