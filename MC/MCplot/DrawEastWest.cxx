void DrawEastWest() {
    TFile *f = TFile::Open("../ooc_ecal_filter.root");
    TTree *vtree = (TTree*)f->Get("vtree");
	TH1F * hEast = new TH1F("East","", 100, 0, 100);
	TH1F * hWest = new TH1F("West","", 100, 0, 100);
    // 定义变量绑定 branch
    float Px, Py, Pz;
    vtree->SetBranchAddress("init_Px", &Px);
    vtree->SetBranchAddress("init_Py", &Py);
    vtree->SetBranchAddress("init_Pz", &Pz);

    Long64_t nentries = vtree->GetEntries();

    // 遍历前 10 个事件作为示例
    for (Long64_t i = 0; i < 100000; ++i) {
        vtree->GetEntry(i);
        float p = sqrt(Px*Px + Py*Py + Pz*Pz);
        // std::cout << "Event " << i << ": Px=" << Px << ", Py=" << Py << ", Pz=" << Pz << ", |p|=" << p << std::endl;
		float deg = acos(Pz/p) * 180.0 / M_PI;
		// std::cout << "deg=" << deg << std::endl;
		float degthr = 20;
		if (deg < degthr && deg > 0) {
		 	hEast->Fill(p/1000);
		}
		if (deg > (180 - degthr) && deg < 180) {
			hWest->Fill(p/1000);
		}
	}

	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	c1->SetLogy();
	c1->SetTitle("");
	c1->SetMargin(0.12,0.1,0.12,0.1);
	c1->SetGrid();
	TLegend *leg = new TLegend(0.55,0.7,0.88,0.88);
	leg->SetBorderSize(0);  // 无边框
	leg->SetFillStyle(0);   // 透明背景
	leg->SetTextSize(0.04);

	leg->AddEntry(hEast, "East (cone < 20^{#circ})", "f");
	leg->AddEntry(hWest, "West (cone < 20^{#circ})", "f");
	leg->Draw();
	TText * t = new TText(0.1,0.9,"");
	hEast->Scale(1/hEast->Integral());
	hWest->Scale(1/hWest->Integral());
	hEast->SetStats(0);
	hWest->SetStats(0);
	hEast->SetLineColor(kRed);
	hWest->SetLineColor(kBlue);
	hEast->SetFillColor(2);
	hWest->SetFillColor(4);
	hEast->SetFillStyle(3003);
	hWest->SetFillStyle(3004);
	
	
	hEast->GetXaxis()->SetTitleSize(0.045);
	hEast->GetYaxis()->SetTitleSize(0.045);
	hEast->GetXaxis()->SetLabelSize(0.04);
	hEast->GetYaxis()->SetLabelSize(0.04);
	hEast->GetXaxis()->SetTitleOffset(1);
	hEast->GetYaxis()->SetTitleOffset(1.1);
	hEast->GetXaxis()->SetLabelOffset(0.01);
	hEast->GetYaxis()->SetLabelOffset(0.01);
	hEast->GetXaxis()->SetTitle("Rigidity [GV]");
	hEast->GetYaxis()->SetTitle("Fraction");
	hEast->Draw("hist");
	hWest->Draw("histsame");
	leg->Draw();
	c1->Update();
    c1->SaveAs("EastWest.png");

    // f->Close();
}
