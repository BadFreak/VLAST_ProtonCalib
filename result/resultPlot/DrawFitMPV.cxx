#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <string>

int ExtractTF1Param(std::string filename, double arr[]) {
    TFile *file = TFile::Open(TString(filename));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return -1;
    }

    // 进入文件夹 ProtonMIPStat
    file->cd("ProtonMIPStat");

    // 遍历 0-24 的 histogram
    for (int i = 0; i < 25; i++) {
        std::string histName = "hProtonSig_CellID_" + std::to_string(i);
        TH1F *hist = (TH1F *)gDirectory->Get(histName.c_str());

        if (!hist) {
            std::cerr << "Histogram " << histName << " not found!" << std::endl;
            continue;
        }

        // 获取拟合函数（默认是 "fit" 或第一个 TF1）
        std::string fitName = "fProtonSig_CellID_" + std::to_string(i);
        TF1 *fit = hist->GetFunction(fitName.c_str());
        if (!fit) {
            // 如果没有 "fit" 名字，尝试取第一个函数
            if (hist->GetListOfFunctions()->GetSize() > 0)
                fit = (TF1 *)hist->GetListOfFunctions()->At(0);
        }

        if (fit) {
            arr[i] = fit->GetParameter(1);
        } else {
            std::cerr << "No TF1 found for " << histName << std::endl;
        }
    }

    file->Close();
    return 0;
}

void DrawFitMPV() {
    // 打开 ROOT 文件
    double fit1[25] = {0.};
    double fit2[25] = {0.};
    ExtractTF1Param("../ProtonCalib_ooc_ecal_filter_5hits_trackfindfit_lengthcor.root", fit1);
    ExtractTF1Param("../ProtonCalib_ooc_ecal_filter_reverse_spectrum.root", fit2);

    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    for (int i = 0; i < 25; i++) {
        g1->SetPoint(i, i, fit1[i]);
        g2->SetPoint(i, i, fit2[i]);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    c1->cd();
    leg->AddEntry(g1, "west", "lp");
    leg->AddEntry(g2, "east", "lp");
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);
    g1->GetXaxis()->SetTitle("CellID");
    g1->GetYaxis()->SetTitle("MPV");
    g1->GetXaxis()->SetLimits(0, 24);
    g1->GetYaxis()->SetRangeUser(105, 125);
    g1->Draw("AP");

    g2->SetMarkerStyle(21);
    g2->SetMarkerColor(kBlue);
    g2->Draw("Psame");
    leg->Draw();
    c1->Update();
}