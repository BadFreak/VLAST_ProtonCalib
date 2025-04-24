void DrawMPVTogether() {
    // std::vector<TFile *> files =
    // {TFile::Open("/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib/result/"
    //                                           "ProtonCalib_ooc_ecal_filter_restriction12.root"),
    //                               TFile::Open("/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib/result/"
    //                                           "ProtonCalib_ooc_ecal_filter_restriction123.root"),
    //                               TFile::Open("/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib/result/"
    //                                           "ProtonCalib_ooc_ecal_filter_restriction124.root")};
    std::vector<TFile *> files = {
        TFile::Open("/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib/result/"
                    "ProtonCalib_ooc_ecal_filter_restriction12.root"),
        TFile::Open("/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib/result/"
                    "ProtonCalib_ooc_ecal_filter_restriction123.root"),
        TFile::Open("/mnt2/USTC/jxwang/VLAST-P/VLAST_ProtonCalib/result/"
                    "ProtonCalib_ooc_ecal_filter_helium_restriction123.root")};

    double mu[3][25] = {0};
    for (size_t fidx = 0; fidx < files.size(); ++fidx) {
        std::cout << ">>> File " << fidx << std::endl;
        // files.at(fidx)->ls();
        TDirectory *dir = (TDirectory *)files[fidx]->Get("ProtonMIPStat");
        for (int i = 0; i < 25; ++i) {
            TString hname = TString::Format("hProtonSig_CellID_%d", i);
            TH1F *h = (TH1F *)dir->Get(hname);
            if (!h) continue;
            TString fname = TString::Format("fProtonSig_CellID_%d", i);
            TF1 *f = h->GetFunction(fname);
            if (!f) continue;
            mu[fidx][i] = f->GetParameter(1);
            // cout << mu[fidx][i] << endl;
        }
    }

    TGraph *gr12 = new TGraph();
    TGraph *gr123 = new TGraph();
    TGraph *gr124 = new TGraph();
    TGraph *gr123_124 = new TGraph();
    for (int i = 0; i < 25; ++i) {
        gr12->SetPoint(i, i, mu[0][i]);
        gr123->SetPoint(i, i, mu[1][i]);
        gr124->SetPoint(i, i, mu[2][i]);
        gr123_124->SetPoint(i, i, (mu[1][i] - mu[2][i]) / mu[2][i]);
    }
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    TLegend *lg = new TLegend(0.7, 0.7, 0.9, 0.9);
    gr12->GetXaxis()->SetTitle("CellID");
    gr12->GetYaxis()->SetTitle("MPV");
    gr12->GetYaxis()->SetRangeUser(108, 120);
    gr123->GetXaxis()->SetTitle("CellID");
    gr123->GetYaxis()->SetTitle("MPV");
    gr123->GetYaxis()->SetRangeUser(108, 120);

    gr12->SetMarkerColor(kRed);
    gr123->SetMarkerColor(kBlue);
    gr124->SetMarkerColor(kGreen);
    gr12->SetMarkerStyle(21);
    gr123->SetMarkerStyle(22);
    gr124->SetMarkerStyle(23);
    gr12->SetMarkerSize(1.5);
    gr123->SetMarkerSize(1.5);
    gr124->SetMarkerSize(1.5);
    gr12->SetLineStyle(2);
    gr123->SetLineStyle(2);
    gr124->SetLineStyle(2);
    // gr12->Draw("ALP");
    // gr123->Draw("LPSAME");
    gr123->Draw("ALP");
    gr124->Draw("LPSAME");
    // lg->AddEntry(gr12, "Restriction 1,2", "lp");
    // lg->AddEntry(gr123, "Restriction 1,2,3", "lp");
    // lg->AddEntry(gr124, "Restriction 1,2,4", "lp");
    lg->AddEntry(gr123, "east", "lp");
    lg->AddEntry(gr124, "west", "lp");
    lg->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    gr123_124->GetXaxis()->SetTitle("CellID");
    // gr123_124->GetYaxis()->SetTitle("(GenFit-Truth)/Truth");
    gr123_124->GetYaxis()->SetTitle("east-west/west");
    // gr123_124->GetYaxis()->SetRangeUser(, 0.02);
    gr123_124->GetYaxis()->SetTitleOffset(1.5);

    gr123_124->SetMarkerColor(kBlack);
    gr123_124->SetMarkerStyle(20);
    gr123_124->SetMarkerSize(1.5);
    gr123_124->SetLineStyle(2);
    gr123_124->Draw("ALP");
}
