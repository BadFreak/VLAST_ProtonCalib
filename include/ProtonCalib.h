#ifndef _PROTONCALIB_H_
#define _PROTONCALIB_H_

#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TMath.h>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

static double LandauConvGausFunc(double *x, double *par) {
    Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
    Double_t mpshift = -0.22278298;      // Landau最大值偏移
    Double_t np = 100.0;
    Double_t sc = 5.0;

    Double_t mpc = par[1] - mpshift * par[0]; // 修正 Landau MPV
    Double_t xlow = x[0] - sc * par[3];
    Double_t xupp = x[0] + sc * par[3];
    Double_t step = (xupp - xlow) / np;

    Double_t sum = 0.0;
    for (Double_t i = 1.0; i <= np / 2; i++) {
        Double_t xx1 = xlow + (i - 0.5) * step;
        Double_t xx2 = xupp - (i - 0.5) * step;
        Double_t fland1 = TMath::Landau(xx1, mpc, par[0]) / par[0];
        Double_t fland2 = TMath::Landau(xx2, mpc, par[0]) / par[0];

        sum += fland1 * TMath::Gaus(x[0], xx1, par[3]);
        sum += fland2 * TMath::Gaus(x[0], xx2, par[3]);
    }

    return par[2] * step * sum * invsq2pi / par[3];
}

class ProtonCalib {
  public:
    static const int nCry = 25;
    void ReadTreeAndSetBranches(std::string filename);
    void ProtonMIPStat(std::string filename);
    int GetPadID(int cellid);
    bool IsLessHit(int hitStd, int hit);
    TVector3 TrackGenFit(TVector3 pos, TVector3 mom, std::vector<TVector3> &trackPos);
    void FitTH1F(TH1F *&h, TF1 *&f);

  private:
    TFile *inFile;
    TFile *outFile;
    TTree *inTree;
    // 存储不同类型的变量
    std::map<std::string, float> floatBranches;               // 单值 float 类型分支
    std::map<std::string, int> intBranches;                   // 单值 int 类型分支
    std::map<std::string, std::vector<int> *> vecIBranches;   // vector<int> 类型分支
    std::map<std::string, std::vector<float> *> vecFBranches; // vector<float> 类型分支

    std::vector<TVector3> hitPos;
    std::vector<float> hitW;

    TH1F *hProtonSig[nCry];
    TH1F *hECALE;
    TH1I *hProtonHit;
    TF1 *fProtonSig[nCry];
};

#endif // _PROTONCALIB_H_