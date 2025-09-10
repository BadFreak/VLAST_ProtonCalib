#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <vector>

void DrawOne() {
    bool isHeliumFile = false;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    std::string filename = "../ProtonCalib_ooc_ecal_filter.root";
    const std::string histPath = "ProtonMIPStat/hProtonSig_CellID_12";
    if (filename.find("helium") != std::string::npos) {
        isHeliumFile = true;
    }
    TText *text0 = new TText(0.64, 0.32, "VLAST-P ECAL");
    TText *text1 = nullptr;
    TCanvas *c1 = new TCanvas("c1", "Multiple Histograms", 800, 600);
    if (isHeliumFile) {
        text1 = new TText(0.6, 0.26, "Helium MIP response");
        std::cout << "Is Helium File : " << isHeliumFile << std::endl;

    } else {
        text1 = new TText(0.6, 0.26, "Proton MIP reponse");
        std::cout << "Not Helium File : " << isHeliumFile << std::endl;
    }
    text0->SetNDC();
    text0->SetTextFont(52);
    text1->SetNDC();
    text1->SetTextFont(52);

    // TLegend *legend = new TLegend(0.5, 0.7, 0.9, 0.9);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetGrid();

    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
    }
    TH1F *hist = (TH1F *)file->Get(histPath.c_str());
    if (!hist) {
        file->Close();
    }

    TF1 *fitfun = hist->GetFunction("fProtonSig_CellID_12");
    TText *text2 = new TText(0.5, 0.8, Form("MPV : %.1f", fitfun->GetParameter(1)));
    text2->SetNDC();

    hist->SetLineColor(kBlue);
    hist->SetLineWidth(0);

    hist->SetStats(0);
    hist->SetFillColor(kBlue);
    hist->SetFillStyle(3001);
    hist->GetXaxis()->SetTitle("E_{Cell} [MeV]");
    if (isHeliumFile) {
        hist->GetXaxis()->SetRangeUser(300, 700);
    } else {
        hist->GetXaxis()->SetRangeUser(60, 240);
    }
    hist->GetYaxis()->SetTitle("Entries");
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->Draw("hist");
    text0->Draw("same");
    text1->Draw("same");
    text2->Draw("same");
    fitfun->SetLineColor(kRed);
    fitfun->SetLineWidth(3);
    fitfun->Draw("same");

    hist->SetTitle("");

    // legend->Draw();
    c1->Update();
    if (isHeliumFile) {
        c1->Print("Cell_12_helium.png");
    } else {
        c1->Print("Cell_12_proton.png");
    }
}
