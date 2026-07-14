// Script to plot traces of a neutron and a gamma detected by one RMD detector
#include <iostream>
#include <fstream>
void PlotPrintRMDTrace_neutronANDgamma(long entryneutron, long entrygamma) {

  int num = 1;
  
  LendaEvent * e1 = new LendaEvent();
  LendaEvent * e2 = new LendaEvent();

  caltree->SetBranchAddress("lendaevent",&e1);
  caltree->GetEntry(entryneutron);

  caltree->SetBranchAddress("lendaevent",&e2);
  caltree->GetEntry(entrygamma);


  
  int Bars1 = e1->NumBars;
  int Bars2 = e2->NumBars;

  cout << "Plotting neutron trace for RMD bar out of " << Bars1 << " bars" << endl;
  cout << "Plotting gamma trace for RMD bar out of " << Bars2 << " bars" << endl;
    
  string name1 = e1->Bars[num].Name;
  string name2 = e2->Bars[num].Name;

  cout << "Neutron trace for bar name: " << name1 << endl;;
  cout << "Gamma trace for bar name: " << name2 << endl;;

  LendaChannel ch1 = e1->Bars[num].Tops[0];
  double energy1 = ch1->GetEnergy();


  LendaChannel ch2 = e2->Bars[num].Tops[0];
  double energy2 = ch2->GetEnergy();



  
  cout << " " << endl;
  cout << "Uncalibrated neutron RMD energy: " << energy1 << " taken from Pulse-height minus baseline" << endl;
  cout << "Uncalibrated gamma RMD energy: " << energy2 << " taken from Pulse-height minus baseline" << endl;

  vector <UShort_t> trace1 = ch1->GetTrace();
  int size1 = trace1.size();

  vector <UShort_t> trace2 = ch2->GetTrace();
  int size2 = trace2.size();


  LendaFilter aFilter1;
  int start;
  double shortGate = aFilter1.GetGateRMD(trace1,start,5);
  double longGate = aFilter1.GetGateRMD(trace1,start,100);
  double psd1 = (longGate-shortGate) / (longGate);
  cout << "Short gate neutron: " << shortGate << " started at: " << start << endl; 
  cout << "Long gate neutron: " << longGate << " started at: " << start << endl; 
  cout << " " << endl;

  LendaFilter aFilter2;
  int start;
  double shortGate = aFilter1.GetGateRMD(trace2,start,5);
  double longGate = aFilter1.GetGateRMD(trace2,start,100);
  double psd2 = (longGate-shortGate)  / (longGate);
  cout << "Short gate gamma: " << shortGate << " started at: " << start << endl; 
  cout << "Long gate gamma: " << longGate << " started at: " << start << endl; 
  cout << " " << endl;


  cout << "psd neutron: " << psd1 << " psd gamma: " << psd2 << " Ratio: " << 100*TMath::Abs(psd1-psd2)/(psd1+psd2) << "%" << endl;
  

  // Do all the plotting stuff
  TCanvas* canv;
  canv = new TCanvas("cTrace");

  
  // Plot trace neutron  
  Double_t *x1 = (Double_t*)malloc(size1*sizeof(Double_t));
  for (int i=0; i < size1; i++){
    x1[i] = i;
  }

  Double_t *y1 = (Double_t*)malloc(size1*sizeof(Double_t));
  for (int i=0; i < size1; i++){
    y1[i] = trace1[i];
  }
  
  TGraph * theTraceGraph1 = new TGraph(size1,x1,y1);  
  
  free(y1);
  free(x1);


  // Plot trace gamma 
  Double_t *x2 = (Double_t*)malloc(size2*sizeof(Double_t));
  for (int i=0; i < size2; i++){
    x2[i] = i;
  }

  Double_t *y2 = (Double_t*)malloc(size2*sizeof(Double_t));
  for (int i=0; i < size2; i++){
    y2[i] = trace2[i];
  }
  
  TGraph * theTraceGraph2 = new TGraph(size2,x2,y2);  
  
  free(y2);
  free(x2);







  //theTraceGraph1->GetHistogram()->GetXaxis()->SetRangeUser(0,750);
  theTraceGraph1->GetHistogram()->GetXaxis()->SetRangeUser(320,480);
  //theTraceGraph1->GetHistogram()->GetXaxis()->SetRangeUser(350,410);
  theTraceGraph1->SetMarkerColor(2);
  theTraceGraph1->SetLineColor(2);
  theTraceGraph1->Draw("AL*");

  theTraceGraph2->SetMarkerColor(4);
  theTraceGraph2->SetLineColor(4);
  theTraceGraph2->Draw("L*");

  //canv->SetLogy();

  //ofstream myfile;
  //myfile.open("trace_neutron5.txt");
  
  //cout << "Timestamp and trace" << endl;
  //for (int i=0; i<size; i++) {
  //   myfile << i << "   " << trace[i] <<  endl;
  //}

  //myfile.close();



  
}
  
