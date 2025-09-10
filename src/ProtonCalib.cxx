#include "ProtonCalib.h"
#include "AbsTrackRep.h"
#include "ConstField.h"
#include "FieldManager.h"
#include "KalmanFitter.h"
#include "MaterialEffects.h"
#include "MeasurementFactory.h"
#include "MeasurementOnPlane.h"
#include "PlanarMeasurement.h"
#include "RKTrackRep.h"
#include "SpacepointMeasurement.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "Track.h"
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TObjArray.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
// #include <genfit/Exception.h>
#include <iostream>
#include <string>
#include <tuple>

bool ProtonCalib::IsHeFile(std::string filename) {
    bool isHeliumFile = false;
    if (filename.find("helium") != std::string::npos) {
        isHeliumFile = true;
    }
    return isHeliumFile;
}

void ProtonCalib::ReadTreeAndSetBranches(std::string filename) {
    // 打开 ROOT 文件
    if (IsHeFile(filename)) {
        std::cout << "is helium file :" << IsHeFile(filename) << "   ";
    }
    inFile = TFile::Open(TString(filename), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        return;
    }
    // 获取 vtree
    inTree = (TTree *)inFile->Get("vtree");
    if (!inTree) {
        std::cerr << "Error: Cannot find TTree named 'vtree'" << std::endl;
        inFile->Close();
        return;
    }
    std::cout << "Reading branches from vtree: " << inTree->GetName() << std::endl;

    // 获取所有 Branch
    TObjArray *branches = inTree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch *branch = (TBranch *)branches->At(i);
        std::string branchName = branch->GetName();
        TLeaf *leaf = branch->GetLeaf(branchName.c_str());

        if (!leaf) {
            std::cerr << "Warning: Cannot determine type of branch " << branchName << std::endl;
            continue;
        }

        std::string typeName = leaf->GetTypeName();
        // std::cout << "Setting Branch: " << branchName << " (Type: " << typeName << ")" <<
        // std::endl;

        // 处理不同类型的 Branch
        if (typeName == "Float_t") {
            floatBranches[branchName] = 0.0f;
            inTree->SetBranchAddress(branchName.c_str(), &floatBranches[branchName]);
        } else if (typeName == "Int_t") {
            intBranches[branchName] = 0;
            inTree->SetBranchAddress(branchName.c_str(), &intBranches[branchName]);
        } else if (typeName == "vector<int>") {
            vecIBranches[branchName] = nullptr;
            inTree->SetBranchAddress(branchName.c_str(), &vecIBranches[branchName]);
        } else if (typeName == "vector<float>") {
            vecFBranches[branchName] = nullptr;
            inTree->SetBranchAddress(branchName.c_str(), &vecFBranches[branchName]);
        } else {
            std::cerr << "Warning: Unsupported branch type: " << typeName << std::endl;
        }
    }
}

void ProtonCalib::InitHist(std::string filename) {
    for (int i = 0; i < nCry; i++) {
        if (IsHeFile(filename)) {
            hProtonSig[i] = new TH1F(Form("hProtonSig_CellID_%d", i),
                                     Form("hProtonSig_CellID_%d", i), 400, 0, 1000);
            fProtonSig[i] =
                new TF1(Form("fProtonSig_CellID_%d", i), LandauConvGausFunc, 380, 700, 4);
        } else {
            hProtonSig[i] = new TH1F(Form("hProtonSig_CellID_%d", i),
                                     Form("hProtonSig_CellID_%d", i), 150, 0, 300);
            fProtonSig[i] =
                new TF1(Form("fProtonSig_CellID_%d", i), LandauConvGausFunc, 80, 300, 4);
        }
    }
    hAngle = new TH1F("hAngle", "hAngle", 60, -1, 14);
    hECALE = new TH1F("hECALE", "hECALE", 400, 0, 800);
    hProtonHit = new TH1I("hProtonHit", "hProtonHit", 25, 0, 25);
    hTrackerHit = new TH1I("hTrackerHit", "hTrackerHit", 15, 0, 15);
    hFitMPV = new TH1F("hFitMPV", "hFitMPV", 500, 50, 550);
    hTest = new TH1F("hAngleDiff", "hAngleDiff", 100, -1, 1);

    grFrontSur = new TGraph();
    grBackSur = new TGraph();
}

double ProtonCalib::PointLineDistance(const TVector3 &point, const TVector3 &line_point,
                                      const TVector3 &line_dir) {
    // d = |(P - A) x dir| / |dir|
    TVector3 diff = point - line_point;
    TVector3 cross = diff.Cross(line_dir);
    return cross.Mag() / line_dir.Mag();
}

double ProtonCalib::CalculChi2(std::vector<TVector3> tempTrackPos) {
    TVector3 A = tempTrackPos[0];
    TVector3 B = tempTrackPos[1];
    TVector3 C = tempTrackPos[2];
    TVector3 center = (A + B + C) * (1.0 / 3.0);
    TVector3 dir = ((B - A) + (C - A)).Unit(); // 简单方向
    double chi2 = 0;
    for (auto &p : tempTrackPos) {
        double d = PointLineDistance(p, center, dir);
        chi2 += d * d; // 如果有误差，可以除以 sigma^2
    }
    return chi2;
}

bool ProtonCalib::TrackFind(std::vector<TVector3> &trackPos) {
    if (trackPos.size() > 5 || trackPos.size() < 3) return false;
    // 3 Tracker Layers
    std::vector<std::vector<TVector3>> layerTrackPos(3);
    double nowZpos = trackPos[0].Z();
    int nowLayer = 0;
    std::sort(trackPos.begin(), trackPos.end(),
              [](const TVector3 &a, const TVector3 &b) { return a.Z() < b.Z(); });

    for (int i = 0; i < trackPos.size(); i++) {
        if (trackPos[i].Z() > nowZpos) {
            nowLayer++;
            nowZpos = trackPos[i].Z();
        }
        layerTrackPos[nowLayer].push_back(trackPos[i]);
    }
    if (layerTrackPos[0].size() == 0 || layerTrackPos[1].size() == 0 ||
        layerTrackPos[2].size() == 0)
        return false;

    double bestChi2 = 1e9;
    std::vector<TVector3> bestTrack;
    for (int i = 0; i < layerTrackPos[0].size(); ++i) {
        for (int j = 0; j < layerTrackPos[1].size(); ++j) {
            for (int k = 0; k < layerTrackPos[2].size(); ++k) {
                std::vector<TVector3> tempTrack = {layerTrackPos[0][i], layerTrackPos[1][j],
                                                   layerTrackPos[2][k]};
                double chi2 = CalculChi2(tempTrack);
                if (chi2 < bestChi2) {
                    bestChi2 = chi2;
                    bestTrack = tempTrack;
                }
            }
        }
    }
    trackPos = bestTrack;
    return true;
}

TVector3 ProtonCalib::TrackGenFit(TVector3 pos, TVector3 mom, std::vector<TVector3> &trackPos) {
    // TVector3 hitDirectionUnit = (trackPos[2] - trackPos[0]).Unit();
    // TVector3 momUnit = mom.Unit();
    // hTest->Fill(hitDirectionUnit.Dot(momUnit));
    // if (hitDirectionUnit.Dot(momUnit) < 0.5) {
    //     return TVector3(-999, -999, -999);
    // }
    using namespace genfit;
    // std::cout << "Init momentum: " << mom.X() << " " << mom.Y() << " " << mom.Z() <<
    // std::endl;
    genfit::MaterialEffects::getInstance()->setNoEffects();
    genfit::AbsBField *field = new genfit::ConstField(0., 0., 0.); // 1 Tesla 沿 Z 轴
    genfit::FieldManager::getInstance()->init(field);              // GenFit 会接管 field 的内存管理
    // 初始状态：RKTrackRep + Track
    AbsTrackRep *rep = new RKTrackRep(2212); // 2212: proton
    Track *track = new Track(rep, pos, mom);
    const int detId(0); // detector ID
    int planeId(0);     // detector plane ID
    int hitId(0);       // hit ID

    // 遍历 trackPos 中的每个点
    for (size_t i = 0; i < trackPos.size(); ++i) {
        // std::cout << "TrackPos: " << trackPos[i].X() << " " << trackPos[i].Y() << " "
        // << trackPos[i].Z() << std::endl;
        TVector3 mpos3D = trackPos[i];
        TVector3 u(1, 0, 0); // x方向
        TVector3 v(0, 1, 0); // y方向
        TVectorD mpos2D(2);
        mpos2D[0] = mpos3D.Dot(u);
        mpos2D[1] = mpos3D.Dot(v);

        TMatrixDSym measCov(2);
        measCov.UnitMatrix();
        measCov *= 0.01; // 设定误差平方 (e.g. 1 mm^2)
        int detId = i;
        PlanarMeasurement *meas = new PlanarMeasurement(mpos2D, measCov, detId, ++hitId, nullptr);
        meas->setPlane(
            genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0, 0, mpos3D.Z()), u, v)),
            ++planeId); // 平面ID(自动递增));
        track->insertPoint(new genfit::TrackPoint(meas, track));
    }

    KalmanFitter *fitter = new KalmanFitter();
    try {
        fitter->processTrack(track);
    } catch (genfit::Exception &e) {
        std::cerr << "[WARNING] Caught unknown exception during fitting. Skipping this track."
                  << std::endl;
        return TVector3(-999, -999, -999);
    }
    // track->getFittedState().Print();
    TVector3 fitMom = track->getFittedState().getMom();
    // std::cout << "Fitted momentum: " << fitMom.X() << " " << fitMom.Y() << " " << fitMom.Z()
    // << std::endl;
    delete fitter; // 释放内存
    return fitMom;
}

TVector3 ProtonCalib::FindHitPoint(TVector3 pos, TVector3 mom, double z) {
    mom = mom.Unit();
    TVector3 hitPoint = pos + (z - pos.Z()) / mom.Z() * mom;
    return hitPoint;
}

bool ProtonCalib::IsInDetectorSur(TVector3 vec3) {
    if (fabs(vec3.X()) > 15 || fabs(vec3.Y()) > 15)
        return false;
    else
        return true;
}

bool ProtonCalib::intersectBox(const TVector3 &rayOrigin, const TVector3 &rayDirInv,
                               const TVector3 &boxMin, const TVector3 &boxMax, double &tmin,
                               double &tmax) {
    double tx1 = (boxMin.X() - rayOrigin.X()) * rayDirInv.X();
    double tx2 = (boxMax.X() - rayOrigin.X()) * rayDirInv.X();
    double ty1 = (boxMin.Y() - rayOrigin.Y()) * rayDirInv.Y();
    double ty2 = (boxMax.Y() - rayOrigin.Y()) * rayDirInv.Y();
    double tz1 = (boxMin.Z() - rayOrigin.Z()) * rayDirInv.Z();
    double tz2 = (boxMax.Z() - rayOrigin.Z()) * rayDirInv.Z();

    tmin = std::max({std::min(tx1, tx2), std::min(ty1, ty2), std::min(tz1, tz2)});
    tmax = std::min({std::max(tx1, tx2), std::max(ty1, ty2), std::max(tz1, tz2)});

    return tmax >= tmin && tmax > 0;
}

void ProtonCalib::RayTraceCrystals(TVector3 pos, TVector3 mom) {
    const int N = 5;
    const double sizeX = 60.0;
    const double sizeY = 60.0;
    const double sizeZ = 200.0;
    const double gap = 62. + 1.5; // CryGap = 1.5 * mm; Hole_XY = 62. * mm;
    const double gapX = gap;
    const double gapY = gap;
    if (intersectResults.size()) intersectResults.clear();

    TVector3 dirNorm = mom.Unit();
    if (dirNorm.Z() == 1) return;
    TVector3 invDir(1.0 / dirNorm.X(), 1.0 / dirNorm.Y(), 1.0 / dirNorm.Z());

    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            double cx = (ix - 2) * gapX;
            double cy = (iy - 2) * gapY;
            double cz = 0.0;

            TVector3 center(cx, cy, cz);
            TVector3 boxMin = center - TVector3(sizeX / 2, sizeY / 2, sizeZ / 2);
            TVector3 boxMax = center + TVector3(sizeX / 2, sizeY / 2, sizeZ / 2);

            double tmin, tmax;
            if (intersectBox(pos, invDir, boxMin, boxMax, tmin, tmax)) {
                double length = (tmax - tmin) * dirNorm.Mag();
                intersectResults.emplace_back(ix, iy, length);
            }
        }
    }

    // std::cout << "Intersected Crystals:\n";
    // for (auto &[ix, iy, len] : intersectResults) {
    //     std::cout << "Crystal (" << ix << ", " << iy << ") - Path Length: " << len << "
    //     mm\n";
    // }
}
bool ProtonCalib::IsLessHit(int hitStd, int hit) {
    return (hit <= hitStd) && (hit > 0);
}

void ProtonCalib::ProtonSelection(double fileweight) {
    std::cout << "======== Starting : Proton MIP Stat" << std::endl;
    int totalEntries = inTree->GetEntries();
    std::cout << "Total entries: " << totalEntries << std::endl;
    for (int i = 0; i < totalEntries; i++) {
        // for (int i = 0; i < 1000000; i++) {
        if (i % 1000000 == 0) {
            std::cout << "Processing entry " << i << std::endl;
        }
        inTree->GetEntry(i);
        // Read ecal_cellid and ecal_celle
        float totalECALE = 0.;

        // jiaxuan selection 1: trigger logic: 3 detectors threshold cut
        float acd_e_0 = vecFBranches["acd_e"]->at(0);
        float acd_e_1 = vecFBranches["acd_e"]->at(1);
        bool is_acd_not = (acd_e_0 < 0.6) || (acd_e_1 < 0.6);
        float max_conv_e =
            !vecFBranches["conv_e"]->empty()
                ? *std::max_element(vecFBranches["conv_e"]->begin(), vecFBranches["conv_e"]->end())
                : 0.0f;
        float max_ecal_celle = !vecFBranches["ecal_celle"]->empty()
                                   ? *std::max_element(vecFBranches["ecal_celle"]->begin(),
                                                       vecFBranches["ecal_celle"]->end())
                                   : 0.0f;

        if (is_acd_not || max_conv_e < 2.52 || max_ecal_celle < 36) continue;
        eventSelection1++;

        // jiaxuan selection 2: ECAL hitNo cut > 0 && <= xx
        int totalECALHit = vecIBranches["ecal_cellid"]->size();
        TVector3 position(floatBranches["init_x"], floatBranches["init_y"],
                          floatBranches["init_z"]);
        TVector3 momentum(floatBranches["init_Px"], floatBranches["init_Py"],
                          floatBranches["init_Pz"]);
        TVector3 frontPoint = FindHitPoint(position, momentum, -100);
        TVector3 backPoint = FindHitPoint(position, momentum, 100);
        grFrontSur->AddPoint(frontPoint.X(), frontPoint.Y());
        grBackSur->AddPoint(backPoint.X(), backPoint.Y());
        hProtonHit->Fill(totalECALHit, fileweight);
        // if (IsInDetectorSur(frontPoint) && IsInDetectorSur(backPoint)) {
        //     hProtonHit->Fill(totalECALHit);
        // } else {
        //     continue;
        // }
        if (!IsLessHit(5, totalECALHit)) continue;
        RayTraceCrystals(position, momentum);
        eventSelection2++;

        // jiaxuan selection 3a(truth level): incident particle direction
        float momentumValue = momentum.Mag();
        float norm_z = floatBranches["init_Pz"] / momentumValue;
        if (norm_z < 0.7624) continue; // 20/sqrt(2*6^2+20^2)
        eventSelection3++;

        // jiaxuan selection 3b: Track find and fit
        // int trackPoint = vecFBranches["tracker_hitx"]->size();
        // if (trackPoint != vecFBranches["tracker_hity"]->size() ||
        //     trackPoint != vecFBranches["tracker_hitz"]->size() || trackPoint == 0) {
        //     continue;
        // }
        // hTrackerHit->Fill(trackPoint,fileweight);
        // if (hitPos.size() != 0) hitPos.clear();
        // for (int i_track_hit = 0; i_track_hit < trackPoint; i_track_hit++) {
        //     float posx = vecFBranches["tracker_hitx"]->at(i_track_hit);
        //     float posy = vecFBranches["tracker_hity"]->at(i_track_hit);
        //     float posz = vecFBranches["tracker_hitz"]->at(i_track_hit);
        //     TVector3 vecPos(posx, posy, posz);
        //     hitPos.push_back(vecPos);
        // }
        // if (TrackFind(hitPos)) {
        //     // if (1) {
        //     TVector3 nowMom(TrackGenFit(position, momentum, hitPos));
        //     double angle = 180 / TMath::Pi() * nowMom.Angle(momentum);
        //     hAngle->Fill(angle,fileweight);
        //     if (nowMom.Z() / nowMom.Mag() < 0.7624) continue;
        //     eventSelection3++;
        // } else {
        //     continue;
        // }

        // jiaxuan selection 4 : ECAL energy restriction
        double PL[25] = {0.};
        double PLCor[25] = {0.};
        for (const auto &[ix, iy, length] : intersectResults) {
            int cellid = 5 * ix + iy;
            PL[cellid] = length;
            PLCor[cellid] = 200.0 / length;
        }
        for (int j = 0; j < totalECALHit; j++) {
            int cellid = vecIBranches["ecal_cellid"]->at(j);
            float celle = vecFBranches["ecal_celle"]->at(j);
            if (PLCor[cellid] < 3.) { // short distance abandoned
                hProtonSig[cellid]->Fill(celle * PLCor[cellid], fileweight);
            }
            totalECALE += celle;
        }
        hECALE->Fill(totalECALE, fileweight);
        eventSelection4++;
    }

    std::cout << "=== Total Event		:	" << totalEntries << " ===" << std::endl;
    std::cout << "=== Event Selection 1	:	" << eventSelection1 << " ===" << std::endl;
    std::cout << "=== Event Selection 2	:	" << eventSelection2 << " ===" << std::endl;
    std::cout << "=== Event Selection 3	:	" << eventSelection3 << " ===" << std::endl;
    std::cout << "=== Event Selection 4	:	" << eventSelection4 << " ===" << std::endl;
}

void ProtonCalib::FitTH1F(TH1F *&h, TF1 *&f, bool isHeliumFile) {
    double landauMPV = h->GetBinCenter(h->GetMaximumBin());
    double landauWidth = h->GetRMS() / sqrt(2);
    double gausSigma = h->GetRMS() / sqrt(2);
    if (isHeliumFile) {
        f->SetParameters(3, 450, 3000, 8);
        f->SetParLimits(0, 2, 30);
        f->SetParLimits(1, 420, 480);
        f->SetParLimits(3, 0, 40);
        h->Fit(f, "QR", "", 380, 700);
    } else {
        f->SetParameters(3, 120, 3000, 8);
        f->SetParLimits(0, 0.5, 5);
        f->SetParLimits(1, 110, 130);
        f->SetParLimits(3, 0, 15);
        h->Fit(f, "QR", "", 80, 300);
    } // landauMPV = f->GetParameter(1);
    // landauWidth = f->GetParameter(0);
    // gausSigma = f->GetParameter(3);
    // f->SetParameters(landauWidth, 120, 10000, gausSigma);
    // h->Fit(f, "R", "", 50, 300);
}

void ProtonCalib::SetHistStyle(TH1F *&h, TF1 *&f, bool isHeliumFile) {
    h->GetListOfFunctions()->Add(f);
    h->GetXaxis()->SetTitle("E_{cell} [MeV]");
    h->GetYaxis()->SetTitle("Entries");
    if (isHeliumFile) {
        h->GetXaxis()->SetRangeUser(100, 800);
    } else {
        h->GetXaxis()->SetRangeUser(30, 300);
    }
    // h->GetYaxis()->SetRangeUser(10, 10000);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->Draw();
}
int ProtonCalib::GetPadID(int cellid) {
    int column = cellid / 5;
    int row = 4 - cellid % 5;
    return (row * 5 + column + 1);
}

int ProtonCalib::CalibFileWrite(std::string filename) {
    bool isHeliumFile = IsHeFile(filename); // true: He file, false: proton file
    TCanvas *cProtonSig = new TCanvas("cProtonSig", "cProtonSig", 1000, 1000);
    cProtonSig->Divide(5, 5);
    for (int i = 0; i < nCry; i++) {
        cProtonSig->cd(GetPadID(i));
        // gPad->SetLogy();
        FitTH1F(hProtonSig[i], fProtonSig[i], isHeliumFile);
        hFitMPV->Fill(fProtonSig[i]->GetParameter(1));
        SetHistStyle(hProtonSig[i], fProtonSig[i], isHeliumFile);
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
    hAngle->GetXaxis()->SetTitle("Angle [deg]");
    hAngle->GetYaxis()->SetTitle("Count");
    hAngle->GetXaxis()->SetTitleSize(0.045);
    hAngle->GetYaxis()->SetTitleSize(0.045);
    hAngle->GetXaxis()->SetLabelSize(0.045);
    hAngle->GetYaxis()->SetLabelSize(0.045);
    hAngle->GetYaxis()->SetRangeUser(500, 1000000);
    hAngle->Write();
    hECALE->GetXaxis()->SetTitle("E_{cell} [MeV]");
    hECALE->GetYaxis()->SetTitle("Entries");
    hECALE->Write();
    hFitMPV->Write();
    hTrackerHit->GetXaxis()->SetTitle("Tracker Hit Number");
    hTrackerHit->GetYaxis()->SetTitle("Count");
    hTrackerHit->GetXaxis()->SetTitleSize(0.045);
    hTrackerHit->GetYaxis()->SetTitleSize(0.045);
    hTrackerHit->GetXaxis()->SetLabelSize(0.045);
    hTrackerHit->GetYaxis()->SetLabelSize(0.045);
    hTrackerHit->GetYaxis()->SetRangeUser(500, 1000000);
    hTrackerHit->Write();
    hProtonHit->GetXaxis()->SetTitle("ECAL Hit Number");
    hProtonHit->GetYaxis()->SetTitle("Count");
    hProtonHit->GetXaxis()->SetTitleSize(0.045);
    hProtonHit->GetYaxis()->SetTitleSize(0.045);
    hProtonHit->GetXaxis()->SetLabelSize(0.045);
    hProtonHit->GetYaxis()->SetLabelSize(0.045);
    hProtonHit->GetYaxis()->SetRangeUser(500, 500000);
    hProtonHit->Write();
    hTest->Write();
    grFrontSur->SetName("FrontSurHitProfile");
    grBackSur->SetName("BackSurHitProfile");
    grFrontSur->Write();
    grBackSur->Write();
    delete hAngle;
    delete hECALE;
    delete hProtonHit;
    outFile->Close();
    inFile->Close();
    return 0;
}
