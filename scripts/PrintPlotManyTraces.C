// Script to print and plot a trace.
// It needs a raw root file generated from R00TLe
//It only needs an entry number. An offset can be added if desired
#include <vector>
#include <iostream>
#include <fstream>
void PrintPlotManyTraces(long entrymin, long entrymax) {

  int sl = 2; // Select DDAS slot number
  int ch = 1; // Select DDAS channel number
    
  
  DDASEvent *e = new DDASEvent();


  rawtree->SetBranchAddress("ddasevent",&e);


  ofstream myfile;
  myfile.open("kolata_manytraces.txt");

  int num = 0; 
  for (int entry=entrymin; entry<=entrymax; entry++) { // Loop over entries

    //cout << "Entry: " << entry << endl; 
    
    rawtree->GetEntry(entry);
    
    ddaschannel *theChannel = e->GetData()[0]; //Extract first channel of entry event 
    int slot = theChannel->GetSlotID();
    int chan = theChannel->GetChannelID();

    //cout << "Slot: " << slot << "  Channel: " << chan << endl;

    
    if (slot == slot && chan == ch) { 

      num++;
      //cout << "Printing trace number: " << num << endl;
      
      int size = theChannel->GetTrace().size();

      Double_t *x = (Double_t*)malloc(size*sizeof(Double_t));
      Double_t *y = (Double_t*)malloc(size*sizeof(Double_t));
   
      for (int i=0; i<size; i++) {
	x[i] = i;
	y[i] = theChannel->GetTrace()[i];
	
      }


      TGraph * theGraph = new TGraph(size,x,y);
  
  
      stringstream ss;
      ss<<"slot "<<slot<<" channel "<<chan;
      theGraph->SetTitle(ss.str().c_str());
      //theGraph->GetHistogram()->GetXaxis()->SetRangeUser(40,80);
      theGraph->Draw("AL*");

      for (int j=0; j<size; j++) {
	myfile << x[j] << "   " << y[j] << endl;
	//cout << "Trace number: " << num << " x[" << j << "] = " << x[j] << "  y[" << j << "] = " << y[j] << endl;
      }

      free(x);
      free(y);
      
    }

  } // End loop over entries


  cout << "Processed traces: " << num << endl;
  
  myfile.close();
  
}
  
