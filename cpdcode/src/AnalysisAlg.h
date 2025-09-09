#ifndef AnalysisAlg_H_
#define AnalysisAlg_H_

#include "SniperKernel/AlgBase.h"
#include "PodioDataSvc/DataHandle.hh"
#include<fstream>
#include "EventDataModel/MCParticleCollection.h"
#include "EventDataModel/TrackingSimHitCollection.h"
#include "EventDataModel/CaloSimCellCollection.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "json/json.h"

using namespace edm;

class AnalysisAlg : public AlgBase
{
   public :

	AnalysisAlg(const std::string& name);
	virtual ~AnalysisAlg();

	virtual bool initialize();
	virtual bool execute();
	virtual bool finalize();
    private :

	int   m_iEvt;
  string   seed_str;
  int      seed_int;
  const MCParticleCollection     *MCParticles;
  const CaloSimCellCollection    *CaloSimHits;
  const TrackingSimHitCollection *SCDSimHits;
  unsigned int scd_cellcode;
  int scd_trackID;
  Vector3f mcparts_momentum;
  Vector3f mcparts_vertex;
  double Ek;
  int Spectrum_N=74; 
  double *binReco;
  

  TRandom3 *rangen;
  double rate;
  double dt;
  double t;
  TFile *Orbit_file;
  TTree *Orbit;
  Double_t lat_geo, lon_geo, rad_geo, ra_scz, dec_scz, ra_scx, dec_scx, ra_scy, dec_scy, mjd,L;
  Double_t pos[3],vel[3];
  Int_t insaa;
  Double_t azimuth,zenith,is_allowed;
  double theta,phi;
  vector<TH1D>* E_spectrum_v;
  vector<TH1D>* E_spectrum_fortime_v;
  vector<TH1D>* E_spectrum_cut_v;
  vector<TH1D>* E_spectrum_fortime_cut_v;
  
  TH1D* E_spectrum;
  TH1D* E_spectrum0;
  TH1D* E_spectrum1;
  TH1D* E_spectrum2;
  TH1D* E_spectrum3;
  TH1D* E_spectrum4;
  TH1D* E_spectrum_fortime;
  TH1D* E_spectrum_fortime0;
  TH1D* E_spectrum_fortime1;
  TH1D* E_spectrum_fortime2;
  TH1D* E_spectrum_fortime3;
  TH1D* E_spectrum_fortime4;
  TH1D* E_spectrum_cut;
  TH1D* E_spectrum0_cut;
  TH1D* E_spectrum1_cut;
  TH1D* E_spectrum2_cut;
  TH1D* E_spectrum3_cut;
  TH1D* E_spectrum4_cut;
  TH1D* E_spectrum_fortime_cut;
  TH1D* E_spectrum_fortime0_cut;
  TH1D* E_spectrum_fortime1_cut;
  TH1D* E_spectrum_fortime2_cut;
  TH1D* E_spectrum_fortime3_cut;
  TH1D* E_spectrum_fortime4_cut;
  double flux_dampe(double x);
  double flux_ams(double x);

  //backtracing
  void vec_in_detector_to_json(double vx,double vy,double vz,double mjd,double *dec_pos,double *dec_vel,double lat,double lon,double &theta,double &arc);
double FT_Modulus(double arg1, double arg2);
double FT_GMST_rad(double timeUnix);
void	FT_GTOD2Equat(double &x, double &y, double &z, double time);
void	FT_Equat2GTOD(double &x, double &y, double &z, double time);
void	FT_GTOD2Equat(double &x, double &y, double &z, double &vx, double &vy, double &vz, double time);
void	FT_Equat2GTOD(double &x, double &y, double &z, double &vx, double &vy, double &vz, double time);
void xyz2polar(double *xyz,double *r_theta_phi);
void polar2xyz(double *xyz,double *r_theta_phi);
void rotatexyz(double *xyz,double theta,double phi,double *xyz_new);



};

#endif
