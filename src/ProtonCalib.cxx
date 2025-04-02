#include <iostream>
#include <string>
#include "ProtonCalib.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TString.h"
#include "TLeaf.h"
#include <algorithm>
#include "TMath.h"

void ProtonCalib::ReadTreeAndSetBranches(std::string filename) {
    // 打开 ROOT 文件
    inFile = TFile::Open(TString(filename), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        return;
    }

    // 获取 vtree
    inTree = (TTree*)inFile->Get("vtree");
    if (!inTree) {
        std::cerr << "Error: Cannot find TTree named 'vtree'" << std::endl;
        inFile->Close();
        return;
    }

    std::cout << "Reading branches from vtree: " << inTree->GetName() << std::endl;

    // 获取所有 Branch
    TObjArray *branches = inTree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch *branch = (TBranch*)branches->At(i);
        std::string branchName = branch->GetName();
        TLeaf *leaf = branch->GetLeaf(branchName.c_str());

        if (!leaf) {
            std::cerr << "Warning: Cannot determine type of branch " << branchName << std::endl;
            continue;
        }

        std::string typeName = leaf->GetTypeName();
        // std::cout << "Setting Branch: " << branchName << " (Type: " << typeName << ")" << std::endl;

        // 处理不同类型的 Branch
        if (typeName == "Float_t") {
            floatBranches[branchName] = 0.0f;
            inTree->SetBranchAddress(branchName.c_str(), &floatBranches[branchName]);
        } 
        else if (typeName == "Int_t") {
            intBranches[branchName] = 0;
            inTree->SetBranchAddress(branchName.c_str(), &intBranches[branchName]);
        } 
		else if (typeName == "vector<int>") {
		    vecIBranches[branchName] = nullptr;
            inTree->SetBranchAddress(branchName.c_str(), &vecIBranches[branchName]);
		}
        else if (typeName == "vector<float>") {
            vecFBranches[branchName] = nullptr;
            inTree->SetBranchAddress(branchName.c_str(), &vecFBranches[branchName]);
        } 
        else {
            std::cerr << "Warning: Unsupported branch type: " << typeName << std::endl;
        }
    }
	for (int i = 0; i < nCry; i++) {
		hProtonSig[i] = new TH1F(Form("hProtonSig_CellID_%d", i), Form("hProtonSig_CellID_%d", i), 150, 0, 300);
		fProtonSig[i] = new TF1(Form("fProtonSig_CellID_%d", i), LandauConvGausFunc, 80, 150, 4);
	}
	hECALE = new TH1F("hECALE", "hECALE", 400, 0, 800);
	hProtonHit = new TH1I("hProtonHit","hProtonHit", 25, 0, 25);

}

void ProtonCalib::FitTH1F(TH1F* &h, TF1* &f) {
	double landauMPV = h->GetBinCenter(h->GetMaximumBin()); 
	double landauWidth = h->GetRMS() / sqrt(2);
	double gausSigma = h->GetRMS() / sqrt(2);
	f->SetParameters(3, 120, 3000, 8);
	f->SetParLimits(0, 0.5, 5);
	f->SetParLimits(1, 110, 130);
	f->SetParLimits(3, 0, 15);
	h->Fit(f, "R", "", 80, 150);
	// landauMPV = f->GetParameter(1);
	// landauWidth = f->GetParameter(0);
	// gausSigma = f->GetParameter(3);
	// f->SetParameters(landauWidth, 120, 10000, gausSigma);
	// h->Fit(f, "R", "", 50, 300);
}

void ProtonCalib::ProtonMIPStat(std::string filename) {
	std::cout << "======== Starting : Proton MIP Stat" << std::endl;
    int totalEntries = inTree->GetEntries();
	std::cout << "Total entries: " << totalEntries << std::endl;

	int eventSelection1 = 0;
	int eventSelection2 = 0;
	int eventSelection3 = 0;
	int eventSelection4 = 0;
	for (int i = 0; i < totalEntries; i++) {
	// for (int i = 0; i < 1000000; i++) {
        if (i % 1000000 == 0) { std::cout << "Processing entry " << i << std::endl; }
		inTree->GetEntry(i);
		// Read ecal_cellid and ecal_celle
		float totalECALE = 0.;
		
		// jiaxuan selection 1: trigger logic
		float acd_e_0 = vecFBranches["acd_e"]->at(0);
		float acd_e_1 = vecFBranches["acd_e"]->at(1);
		bool is_acd_not = (acd_e_0 < 0.6) || (acd_e_1 < 0.6);
		float max_conv_e = !vecFBranches["conv_e"]->empty() ? *std::max_element(vecFBranches["conv_e"]->begin(), vecFBranches["conv_e"]->end()) : 0.0f;
		float max_ecal_celle = !vecFBranches["ecal_celle"]->empty() ? *std::max_element(vecFBranches["ecal_celle"]->begin(), vecFBranches["ecal_celle"]->end()) : 0.0f;

		if (is_acd_not || max_conv_e < 2.52 || max_ecal_celle < 36) continue;		
		eventSelection1++;
		
		// jiaxuan selection 2:
		int totalECALHit = vecIBranches["ecal_cellid"]->size();
		if (!IsLessHit(1, totalECALHit)) continue;
		eventSelection2++;
		
		// jiaxuan selection 3:
		float momentum = floatBranches["init_Px"] * floatBranches["init_Px"] + floatBranches["init_Py"] * floatBranches["init_Py"] + floatBranches["init_Pz"] * floatBranches["init_Pz"];
		float norm_z = floatBranches["init_Pz"] / sqrt(momentum);
		if (norm_z < 0.9206) continue; // 20/sqrt(2*6^2+20^2)
		eventSelection3++;

		// jiaxuan selection 4: incident particle energy restriction
		// float init_e = floatBranches["init_E"];
		// if (init_e > 45000 || init_e < 35000) continue;
		// eventSelection4++;

		for (int j = 0; j < totalECALHit; j++) {
			int cellid = vecIBranches["ecal_cellid"]->at(j);
			float celle = vecFBranches["ecal_celle"]->at(j);
			hProtonSig[cellid]->Fill(celle);
			totalECALE += celle;
		}
		hECALE->Fill(totalECALE);
		hProtonHit->Fill(vecIBranches["ecal_cellid"]->size());
	}
	std::cout << "=== Total Event		:	" << totalEntries << " ===" << std::endl;
	std::cout << "=== Event Selection 1	:	" << eventSelection1 << " ===" << std::endl;
	std::cout << "=== Event Selection 2	:	" << eventSelection2 << " ===" << std::endl;
	std::cout << "=== Event Selection 3	:	" << eventSelection3 << " ===" << std::endl;
	std::cout << "=== Event Selection 4	:	" << eventSelection4 << " ===" << std::endl;
	hECALE->GetXaxis()->SetTitle("E_{cell} [MeV]");
	hECALE->GetYaxis()->SetTitle("Entries");
	
	TCanvas * cProtonSig = new TCanvas("cProtonSig", "cProtonSig", 1000, 1000);
	cProtonSig->Divide(5, 5);
	for (int i = 0; i < nCry; i++) {
	    cProtonSig->cd(GetPadID(i));
		// gPad->SetLogy();
		FitTH1F(hProtonSig[i], fProtonSig[i]);
		hProtonSig[i]->GetListOfFunctions()->Add(fProtonSig[i]);
		hProtonSig[i]->GetXaxis()->SetTitle("E_{cell} [MeV]");
		hProtonSig[i]->GetYaxis()->SetTitle("Entries");
		hProtonSig[i]->GetXaxis()->SetRangeUser(30, 300);
		// hProtonSig[i]->GetYaxis()->SetRangeUser(10, 10000);
		hProtonSig[i]->GetXaxis()->SetLabelSize(0.045);
		hProtonSig[i]->GetYaxis()->SetLabelSize(0.045);
		hProtonSig[i]->GetXaxis()->SetTitleSize(0.045);
		hProtonSig[i]->GetYaxis()->SetTitleSize(0.045);
		hProtonSig[i]->Draw();
		fProtonSig[i]->Draw("same");
	}

	outFile = new TFile(TString(filename), "RECREATE");
	outFile->cd();
	cProtonSig->Write();
	outFile->mkdir("ProtonMIPStat");
	outFile->cd("ProtonMIPStat");
	for (int i = 0; i < nCry; i++) {
		hProtonSig[i]->Write();
	}
	outFile->cd();
	hECALE->Write();
	hProtonHit->Write();
	outFile->Close();
	inFile->Close();
}

int ProtonCalib::GetPadID(int cellid) {
    int column = cellid / 5;
	int row = 4 - cellid % 5;
	return (row * 5 + column + 1);
}

bool ProtonCalib::IsLessHit(int hitStd, int hit) {
    return (hit <= hitStd) && (hit > 0);
}
// float ProtonCalib::AngleCorrection(float dx, float dy, float dz, float &celle) {
// 	float ratio = ;
//     return celle;
// }

// void ProtonCalib::Draw() {}