#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <vector>

void DrawCellTogether() {
    std::vector<std::string> filenames = {
        "ProtonCalib_ooc_ecal_filter_5hits_trackfindfit_lengthcor.root",
        "ProtonCalib_ooc_ecal_filter_5hits_trackfindfit_nolencor.root"};
    const std::string histPath = "ProtonMIPStat/hProtonSig_CellID_15";
    for_each(filenames.begin(), filenames.end(),
             [](std::string &filename) { filename = "../" + filename; });

    TCanvas *c1 = new TCanvas("c1", "Multiple Histograms", 800, 600);
    c1->SetGrid();
    TLegend *legend = new TLegend(0.48, 0.66, 0.88, 0.86);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(3002);
    TText *text = new TText(0.6, 0.3, "VLAST-P ECAL");
    text->SetNDC();
    text->SetTextFont(52);

    bool first = true;
    int color = 1;

    for (size_t i = 0; i < filenames.size(); ++i) {
        TFile *file = TFile::Open(filenames[i].c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Cannot open file: " << filenames[i] << std::endl;
            continue;
        }

        TH1F *hist = (TH1F *)file->Get(histPath.c_str());
        if (!hist) {
            std::cerr << "Histogram not found in file: " << filenames[i] << std::endl;
            file->Close();
            continue;
        }

        hist->SetLineWidth(1);
        std::string histName = "file" + std::to_string(i + 1);
        switch (i) {
        case 0:
            legend->AddEntry(hist, "w/ Length Correction", "f");
            hist->SetFillColorAlpha(kBlue, 0.1);
            hist->SetLineColor(kBlue);
            hist->SetFillStyle(3007);

            break;
        case 1:
            legend->AddEntry(hist, "w/o Length Correction", "f");
            hist->SetFillColorAlpha(kRed, 0.3);
            hist->SetLineColor(kRed);
            hist->SetFillStyle(3005);

            break;
        default:
            break;
        }

        hist->SetStats(0);

        if (first) {
            hist->GetXaxis()->SetTitle("E_{Cell} [MeV]");
            hist->GetXaxis()->SetRangeUser(20, 250);
            hist->GetYaxis()->SetTitle("Entries");
            hist->GetYaxis()->SetRangeUser(0, 600);
            hist->GetYaxis()->SetTitleOffset(1);
            hist->Draw("hist");
            first = false;
        } else {
            hist->Draw("hist same");
        }
        hist->SetTitle("");

        // 注意：不要关闭文件，否则 hist 会被删除
        // file->Close();
    }
    text->Draw();
    legend->Draw();
    c1->Update();
    c1->Print("Cell_Comparison_15.png");
}
