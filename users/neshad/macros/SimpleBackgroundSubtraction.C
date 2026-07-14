#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <iostream>
#include <iomanip>

void SimpleBackgroundSubtraction()
{
    bool wait = true;

    int XMax = 4000;
    int XMin = 0;
	
	int YMax = 30000;
    int YMin = -1000;

    TString detector = "NL01";

    TString Source_runnum = "5380";
    TString Background_runnum = "5381";

    TString hname = detector + "_BottomE";

    double runtime_src = 1;
    double runtime_bg = 1;

    TFile* f_src = TFile::Open("histograms/HistogramsFromRun" + Source_runnum + ".root");
    if (!(f_src->GetListOfKeys()->Contains(hname))) {
        std::cout << "ERROR: Can't find histogram " << hname << std::endl;
        exit(0);
    }

    TFile* f_bg = TFile::Open("histograms/HistogramsFromRun" + Background_runnum + ".root");
    if (!(f_bg->GetListOfKeys()->Contains(hname))) {
        std::cout << "ERROR: Can't find histogram " << hname << std::endl;
        exit(0);
    }

    // Retrieve histograms
    TH1D* h_src = (TH1D*)f_src->Get(hname);
    h_src->Sumw2();

    TH1D* h_bg = (TH1D*)f_bg->Get(hname);
    h_bg->Sumw2();

    // Scale histograms
    TH1D* h_bg_scale = (TH1D*)h_bg->Clone("h_bg_scale");
    h_bg_scale->Scale(1/runtime_bg);

    TH1D* h_src_scale = (TH1D*)h_src->Clone("h_src_scale");
    h_src_scale->Scale(1/runtime_src);

    // Subtract background
    TH1D* h_subtracted = (TH1D*)h_src_scale->Clone("h_subtracted");
    h_subtracted->Add(h_bg_scale, -1);

    // Calculate total counts for each histogram
    double total_counts_src = h_src_scale->Integral();
    double total_counts_bg = h_bg_scale->Integral();
    double total_counts_subtracted = h_subtracted->Integral();

    // Print total counts to console
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total counts in source histogram: " << total_counts_src << std::endl;
    std::cout << "Total counts in background histogram: " << total_counts_bg << std::endl;
    std::cout << "Total counts in subtracted histogram: " << total_counts_subtracted << std::endl;

    // Draw subtracted histogram
    TCanvas* c = new TCanvas("c", "Subtracted Histogram");
    h_subtracted->SetTitle(hname);
    h_subtracted->GetXaxis()->SetRangeUser(XMin, XMax);
	h_subtracted->GetYaxis()->SetRangeUser(YMin, YMax);
	h_src_scale->SetLineColor(kBlack);
	h_bg_scale->SetLineColor(kBlue);
	h_subtracted->SetLineColor(kRed);
    h_subtracted->Draw("HIST");
	h_src_scale->Draw("HIST SAME");
	h_bg_scale->Draw("HIST SAME");
	
    gStyle->SetOptStat(0);
    // Add total counts as text on the canvas
    TLatex latex;
    latex.SetNDC(); // Use normalized device coordinates
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.5, 0.5, Form("Total counts (Source): %.2f", total_counts_src));
    latex.DrawLatex(0.5, 0.45, Form("Total counts (Background): %.2f", total_counts_bg));
    latex.DrawLatex(0.5, 0.4, Form("Total counts (Subtracted): %.2f", total_counts_subtracted));
	latex.DrawLatex(0.65, 0.3, Form("Run time: %.2f Hr",1.00));
	
	TLegend *legend = new TLegend(0.6, 0.65, 0.88, 0.88);
	legend->SetBorderSize(0);
	legend->SetFillColor(0);
	legend->SetTextSize(0.035);
	legend->AddEntry(h_src_scale, "Source Data", "l");
    legend->AddEntry(h_bg_scale, "Background Data", "l");
    legend->AddEntry(h_subtracted, "Subtracted Data", "l");
    legend->Draw();

    if (wait) { gPad->WaitPrimitive(); } // Wait for user interaction

    // Clean up
    delete c;
    f_bg->Close();
    f_src->Close();
}
