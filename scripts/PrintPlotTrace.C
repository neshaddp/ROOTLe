// Script to print and plot a trace.
// It needs a raw root file generated from R00TLe
//It only needs an entry number. An offset can be added if desired

#include <iostream>
#include <fstream>
void PrintPlotTrace(long entry, int offset=0) {

  DDASEvent *e = new DDASEvent();


  rawtree->SetBranchAddress("ddasevent",&e);
  rawtree->GetEntry(entry);
  
  ddaschannel *theChannel = e->GetData()[0]; //Extract first channel of entry event 
  int slot = theChannel->GetSlotID();
  int chan = theChannel->GetChannelID();

  cout << "Slot: " << slot << " Channel: " << chan
       << " for first ddas chanel in entry = " << entry << endl;
  cout << " " << endl;
  
  int size = theChannel->GetTrace().size();
  
  int *x = (int*)malloc(size*sizeof(int));
  int *y = (int*)malloc(size*sizeof(int));

  for (int i=0;i<size;i++){
    x[i]=i;
    y[i]=theChannel->GetTrace()[i]-offset;
  }
    
  TGraph * theGraph = new TGraph(size,x,y);
  
  
  stringstream ss;
  ss<<"slot "<<slot<<" channel "<<chan;
  theGraph->SetTitle(ss.str().c_str());
  //theGraph->GetHistogram()->GetXaxis()->SetRangeUser(40,80);
  theGraph->Draw("AL*");
 
  ofstream myfile;
  myfile.open("trace_kolata.txt");
  
  for (int i=0; i<size; i++) {
    myfile << x[i] << "   " << y[i] <<  endl;
  }

  
}
  
