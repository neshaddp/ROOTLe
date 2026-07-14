// A script to plot the tracesm, fast filter, and CFD arrays for an event consisting of one single-channel detector
// (e.g. RMD scintillator) and two PMTs channels from one LENDA bar

void PlotAll_Timing(bool isbar, int icanvas, TCanvas *c, LendaChannel *theChannel,int fl,int fg,int w,int d){

  // This is where I create the trace graph
  int size = theChannel->GetTrace().size();
  Double_t *x = (Double_t*)malloc(size*sizeof(Double_t));
  for (int i=0; i < size; i++){
    x[i] = i;
  }

  Double_t *y = (Double_t*)malloc(size*sizeof(Double_t));
  for (int i=0; i < size; i++){
    y[i] = theChannel->GetTrace()[i];
  } 
  TGraph * theTraceGraph = new TGraph(size,x,y);  
  
  if (isbar) {
    theTraceGraph->GetHistogram()->GetXaxis()->SetRangeUser(20,100);
  } else {
    theTraceGraph->GetHistogram()->GetXaxis()->SetRangeUser(20,1500);
    //theTraceGraph->GetHistogram()->GetXaxis()->SetRangeUser(0,size);
  }

  free(y);
  free(x);


  // This is where I create the CFD graph
  LendaFilter aFilter;
  vector <UShort_t> trace=theChannel->GetTrace();
  vector <Double_t> CFD  =aFilter.GetNewFirmwareCFD(trace,fl,fg,d,w);//6,0,3,3);
  int cfdsize = CFD.size();

  
  Double_t *x = (Double_t*)malloc(cfdsize*sizeof(Double_t));
  Double_t *y = (Double_t*)malloc(cfdsize*sizeof(Double_t));
  for (int i=0; i < cfdsize; i++) {
    x[i] = i;
    y[i] = CFD[i];
  }

  TGraph * theCFDGraph = new TGraph(cfdsize,x,y);  

  if (isbar) {
    theCFDGraph->GetHistogram()->GetXaxis()->SetRangeUser(20,100);
  } else {
    //theCFDGraph->GetHistogram()->GetXaxis()->SetRangeUser(700,800);
    theCFDGraph->GetHistogram()->GetXaxis()->SetRangeUser(600,900);
  }


  free(x);
  free(y);





  // This is where I create the FF graph
  vector <Double_t> FF;
  aFilter.FastFilter(trace,FF,fl,fg);
  int ffsize = FF.size();

  cout << "Trace size: " << size << endl;
  cout << "CFD size: " << cfdsize << endl;
  cout << "FF size: " << cfdsize << endl;

  Double_t *x = (Double_t*)malloc(ffsize*sizeof(Double_t));
  Double_t *y = (Double_t*)malloc(ffsize*sizeof(Double_t));
  for (int i=0; i < ffsize; i++) {
    x[i] = i;
    y[i] = FF[i];
  }

  TGraph * theFFGraph = new TGraph(ffsize,x,y);  

  if (isbar) {
    theFFGraph->GetHistogram()->GetXaxis()->SetRangeUser(20,100);
  } else {
    theFFGraph->GetHistogram()->GetXaxis()->SetRangeUser(600,900);
  }




  // And I plot here the graphs created
  c->cd(icanvas + 1);
  theTraceGraph->Draw("AL*");

  c->cd(icanvas + 4);
  theCFDGraph->Draw("AL*");

  c->cd(icanvas + 7);
  theFFGraph->Draw("AL*");





  // I calculate now some quantitites that need to be checked out
  int MaxSpotInTrace_temp = -1;
  int pulseheight = aFilter.GetMaxPulseHeight(trace,MaxSpotInTrace_temp);
  double cubiccfd = -999;
  double gate = -999;
  int start = -1;
  

  if (isbar) {
    gate = aFilter.GetGate(trace,start,30);
    cubiccfd = aFilter.GetZeroCubic(CFD,MaxSpotInTrace_temp);
  } else {
    gate = aFilter.GetGateRMD(trace,start,200);
    cubiccfd = aFilter.GetZeroCubicRMD(CFD,start);
  }
  cout << "Pulse height: " << pulseheight << " Maximum: " << MaxSpotInTrace_temp << " CFD crossing: " << cubiccfd << " Start gate: " << start << endl;












}

//void PlotCFD_ForCal(TTree* theTree,long entry,int fl,int fg,int w,int d){
void PlotAll_Timing(long entry,int fl,int fg,int w,int d){

  TTree* theTree =(TTree*)gDirectory->Get("caltree");

  LendaEvent * e = new LendaEvent();
  
  theTree->SetBranchAddress("lendaevent",&e);
  theTree->GetEntry(entry);
  int Bars = e->NumBars;
  cout << "------" << endl;
  cout << "Numberof bars: " << Bars << endl;

  cout << "Names " << e->Bars[0].Name << " and " << e->Bars[1].Name << endl;

  if ( (Bars == 2) && ((e->Bars[0].Name == "SL01") || (e->Bars[1].Name == "SL01")) ) {

    TCanvas* canv;
    canv = new TCanvas("cCompare");
    canv->Divide(3,3);

    for (int i=0; i<Bars; i++) {
	
	bool barflag = false;

	string name = e->Bars[i].Name;
	cout << "Name: " << name << " index: "<< i << endl;

	if (name == "SL01") {
	  
	  barflag = true;
	  
	  cout << "Plotting SL01T" << endl;
	  LendaChannel ct = e->Bars[i].Tops[0];
	  PlotAll_Timing (barflag,i,canv, &ct,2,2,4,7);//  fl,fg,w,d);

	  cout << "Plotting SL01B" << endl;
	  LendaChannel cb = e->Bars[i].Bottoms[0];
	  PlotAll_Timing (barflag,i+1,canv, &cb,2,2,4,7);//fl,fg,w,d);

	} else {

	  barflag = false;

	  cout << "Plotting RMD" << endl;
	  LendaChannel c = e->Bars[i].Tops[0];
	  PlotAll_Timing (barflag,i,canv, &c,fl,fg,w,d);

	}
	  //PlotAll_Timing (&c,fl,fg,w,d);
     }
  }
    
  //for (int i=0; i<Bars; i++){
  //  cout << "Number of tops: " << e->Bars[i].NumTops<<endl;
  //  cout << "Name " << e->Bars[i].Name << endl;
  //  LendaChannel c = e->Bars[i].Tops[0];
  //  PlotCFD_ForCal (&c,fl,fg,w,d);
  //}


  return;

}
