#include "AnalysisAlg.h"

#include "RootWriter/RootWriter.h"
#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperPtr.h"
#include "json/json.h"
// #include "yaml-cpp/yaml.h"
#include <iomanip>
#define pi 3.1415926535897932384626433832795028841971693993751
using namespace edm;

DECLARE_ALGORITHM(AnalysisAlg);

AnalysisAlg::AnalysisAlg(const std::string &name) : AlgBase(name) {
    m_iEvt = 0;
    declProp("Seed", seed_str);
}

AnalysisAlg::~AnalysisAlg() {}

bool AnalysisAlg::initialize() {
    LogInfo << " initialized successfully" << std::endl;

    Orbit_file = new TFile("/home/chenpd/lustre/backtracing/bin/forHERD/high_resolution/"
                           "HERDOrbitSimulation_20220405.root");
    Orbit = (TTree *)Orbit_file->Get("OrbitTree");
    Orbit->SetBranchAddress("lat_geo", &lat_geo); // latitude
    Orbit->SetBranchAddress("lon_geo", &lon_geo); // longitutde
    Orbit->SetBranchAddress("rad_geo", &rad_geo); // radius
    Orbit->SetBranchAddress("dec_scz", &dec_scz); // Declination of SC-Z
    Orbit->SetBranchAddress("ra_scx", &ra_scx);   // Right Ascension of SC-X
    Orbit->SetBranchAddress("dec_scx", &dec_scx); // Declination of SC-X
    Orbit->SetBranchAddress("ra_scy", &ra_scy);   // Right Ascension of SC-Y
    Orbit->SetBranchAddress("dec_scy", &dec_scy); // Declination of SC-Y
    Orbit->SetBranchAddress("mjd", &mjd);         // 修正儒略日（Modified Julian Date）
    Orbit->SetBranchAddress("L", &L);             // 磁壳参数（McIlwain L-parameter）
    Orbit->SetBranchAddress("insaa", &insaa);     // 南大西洋异常区标志
    Orbit->SetBranchAddress("pos[3]", &pos);      // 探测器在地心坐标系（ECEF）中的位置
    Orbit->SetBranchAddress("vel[3]", &vel);      // 探测器在地心坐标系（ECEF）中的速度

    seed_int = std::stoi(seed_str);
    rangen = new TRandom3(seed_int);
    rate = 10000.0 / 86400.0; // s^-1
    t = 59674;
    dt = 0;
    double EEL_reco = 0.5;
    double EEH_reco = 1000;
    double EEL_reco_Log = TMath::Log10(EEL_reco);
    double EEH_reco_Log = TMath::Log10(EEH_reco);
    double bstep_reco = (EEH_reco_Log - EEL_reco_Log) / Spectrum_N; // 74?
    binReco = new double[Spectrum_N + 1];
    for (int ib = 0; ib < (Spectrum_N + 1); ib++) {
        binReco[ib] = EEL_reco_Log + ib * bstep_reco;
        binReco[ib] = TMath::Power(10, binReco[ib]);
    }
    E_spectrum = new TH1D("E_spectrum", "E_spectrum", Spectrum_N, binReco);
    E_spectrum_fortime = new TH1D("E_spectrum_fortime", "E_spectrum_fortime", Spectrum_N, binReco);
    E_spectrum_cut = new TH1D("E_spectrum_cut", "E_spectrum_cut", Spectrum_N, binReco);
    E_spectrum_fortime_cut = new TH1D("E_spectrum_fortime_cut", "E_spectrum_fortime_cut", Spectrum_N, binReco);
    for (int isp = 0; isp < 5; ++isp) {
        E_spectrum_v[isp] = new TH1DF(Form("E_spectrum%d", isp), Form("E_spectrum%d", isp), Spectrum_N, binReco);
        E_spectrum_fortime_v[isp] =
            new TH1D(Form("E_spectrum_fortime%d", isp), Form("E_spectrum_fortime%d", isp), Spectrum_N, binReco);
        E_spectrum_cut_v[isp] =
            new TH1D(Form("E_spectrum_cut%d", isp), Form("E_spectrum_cut%d", isp), Spectrum_N, binReco);
        E_spectrum_fortime_cut_v[isp] =
            new TH1D(Form("E_spectrum_fortime_cut%d", isp), Form("E_spectrum_fortime_cut%d", isp), Spectrum_N, binReco);
    }
    return true;
}

bool AnalysisAlg::execute() {
    //    LogDebug << "Processing event " << m_iEvt << std::endl;
    ++m_iEvt;
    if (m_iEvt % 5000 == 0) cout << m_iEvt << endl;

    CaloSimHits = getROColl(CaloSimCellCollection, "calohits");
    MCParticles = getROColl(MCParticleCollection, "mcparts");
    SCDSimHits = getROColl(TrackingSimHitCollection, "scdhits");
    scd_trackID = -1;
    if (SCDSimHits) {
        for (size_t i = 0; i < SCDSimHits->size(); i++) {
            scd_cellcode = SCDSimHits->at(i).getCellCode();
            scd_trackID = SCDSimHits->at(i).getTrackID();
            if (scd_trackID == 1) break;
        }
    }
    if (MCParticles) {
        mcparts_momentum = MCParticles->at(0).getMomentum();
        mcparts_vertex = MCParticles->at(0).getVertex();
        Ek = sqrt(mcparts_momentum[0] * mcparts_momentum[0] + mcparts_momentum[1] * mcparts_momentum[1] +
                  mcparts_momentum[2] * mcparts_momentum[2] + 0.511e-3 * 0.511e-3) -
             0.511e-3;
        theta = acos(mcparts_momentum[2] /
                     sqrt(mcparts_momentum[0] * mcparts_momentum[0] + mcparts_momentum[1] * mcparts_momentum[1] +
                          mcparts_momentum[2] * mcparts_momentum[2]));
        phi = atan2(mcparts_momentum[1], mcparts_momentum[0]);
    }
    /*  if(Ek>55&&Ek<=100){
            E_spectrum->Fill(Ek,4*pow(pi,2)*pow(1.8,2)*86400*365*2*Ek*(TMath::Log(100)-TMath::Log(2))*flux_dampe(Ek)/200000000);
            }
            if(Ek>100&&Ek<=1000){
            E_spectrum->Fill(Ek,4*pow(pi,2)*pow(1.8,2)*86400*365*2*Ek*(TMath::Log(1000)-TMath::Log(100))*flux_dampe(Ek)/5000000);
            }
            if(Ek>1000&&Ek<=20000){
            E_spectrum->Fill(Ek,4*pow(pi,2)*pow(1.8,2)*86400*365*2*Ek*(TMath::Log(20000)-TMath::Log(1000))*flux_dampe(Ek)/3000000);
            }*/
    double Ek_weight = 4 * pow(pi, 2) * pow(1.8, 2) * 86400 * 365 * 2 * Ek * (TMath::Log(100) - TMath::Log(2)) *
                       flux_ams(Ek) / 3000000;
    if (Ek > 0.5 && Ek <= 100) {
        E_spectrum->Fill(Ek, Ek_weight);
        if (theta > asin(6471. / 6771.)) {
            E_spectrum_cut->Fill(Ek, Ek_weight);
        }

        if (scd_trackID == 1) {
            int cell_code = scd_cellcode / (int)1e6; // 计算 cell_code 只做一次
            bool valid_condition = false;

            switch (cell_code) {
            case 0: valid_condition = (mcparts_momentum[2] < 0); break;
            case 1: valid_condition = (mcparts_momentum[0] < 0); break;
            case 2: valid_condition = (mcparts_momentum[1] < 0); break;
            case 3: valid_condition = (mcparts_momentum[0] > 0); break;
            case 4: valid_condition = (mcparts_momentum[1] > 0); break;
            default: break;
            }

            if (valid_condition) {
                E_spectrum[cell_code]->Fill(Ek, Ek_weight);
                if (theta > asin(6471. / 6771.)) {
                    E_spectrum_cut[cell_code]->Fill(Ek, Ek_weight);
                }
            }
        }
    }

    double uu = rangen->Uniform();
    double theta_in_velocity;
    double arc_in_velocity;
    dt = -log(1 - uu) / rate;
    t = t + dt / 86400;
    if (t <= 59675) {
        if (m_iEvt % 5000 == 0) cout << setprecision(20) << "t=" << t << endl;
        int orbit_entries = (int)((t - 59674) * 8640000); // faster
        Orbit->GetEntry(orbit_entries);
        vec_in_detector_to_json(mcparts_momentum[0], mcparts_momentum[1], mcparts_momentum[2], mjd, pos, vel, lat_geo,
                                lon_geo, theta_in_velocity, arc_in_velocity);
        double lon;
        double lat;
        for (int i = -45; i <= 45; i = i + 5) { // lat
            if (abs(lat_geo - i) <= 2.5) {
                lat = i;
            }
        }
        if (lon_geo > 180) {
            lon = lon_geo - 360;
        } else {
            lon = lon_geo;
        }

        if (lon < -175 || lon > 175) {
            lon = 180;
        } else {
            lon = round(lon / 10.0) * 10; // 直接四舍五入到最近的10的倍数
        }

        double bin_edge[35] = {2,       2.16272, 2.42661, 2.7227,  3.05492, 3.42768, 3.84592, 4.31519, 4.84172,
                               5.4325,  6.09537, 6.83912, 7.67361, 8.60994, 9.66051, 10.8393, 12.1619, 13.6458,
                               15.3109, 17.1791, 19.2752, 21.6272, 24.2661, 27.227,  30.5492, 34.2768, 38.4592,
                               43.1519, 48.4172, 54.325,  60.9537, 68.3912, 76.7361, 86.0994, 100};
        ifstream theta_arc("/home/chenpd/lustre/herd_code/herdos_cutoff/src/theta_arc.txt", ios::in);
        double data[3];
        int theta_line_number;
        double dis = 999;
        for (int i = 0; i < 13963; i++) {
            for (int j = 0; j < 3; j++) {
                theta_arc >> data[j];
            }
            if ((pow((theta_in_velocity - data[1]), 2) + pow((arc_in_velocity - data[2]), 2)) < dis) {
                dis = pow((theta_in_velocity - data[1]), 2) + pow((arc_in_velocity - data[2]), 2);
                theta_line_number = i;
            }
        }
        theta_arc.close();
        Json::Reader reader;
        Json::Value value;
        string str_lon = to_string((int)lon);
        string str_lat = to_string((int)lat);
        int line_number = -1;
        for (int i = 0; i < 34; i++) {
            if (Ek < bin_edge[i + 1] && Ek >= bin_edge[i]) {
                line_number = i;
                break;
            }
        }
        string str_line_number = to_string(line_number);
        string str_json = "/home/chenpd/lustre/backtracing_+_jsonfile/small/result_lon_" + str_lon + "_lat_" + str_lat +
                          "_" + str_line_number + ".json";
        ifstream is(str_json, ios::in);
        if (!is.is_open()) {
            cout << "fail to open file" << endl;
            cout << str_json << endl;
            cout << "lat_geo=" << lat_geo << endl;
            cout << "lon_geo=" << lon_geo << endl;
        }
        if (!reader.parse(is, value)) cout << "error" << endl;
        is_allowed = value["allowed_data"][theta_line_number]["is_allowed"].asDouble();
        if (is_allowed == 1 && Ek > 0.5 && Ek <= 100) {
            double Ek_weight = 4 * pow(pi, 2) * pow(1.8, 2) * 86400 * 365 * 2 * Ek * (TMath::Log(100) - TMath::Log(2)) *
                               flux_ams(Ek) / 3000000;

            E_spectrum_fortime->Fill(Ek, Ek_weight);
            if (theta > asin(6471. / 6771.)) {
                E_spectrum_fortime_cut->Fill(Ek, Ek_weight);
            }

            if (scd_trackID == 1) {
                int cell_code = scd_cellcode / (int)1e6;
                bool valid_condition = false;

                switch (cell_code) {
                case 0: valid_condition = (mcparts_momentum[2] < 0); break;
                case 1: valid_condition = (mcparts_momentum[0] < 0); break;
                case 2: valid_condition = (mcparts_momentum[1] < 0); break;
                case 3: valid_condition = (mcparts_momentum[0] > 0); break;
                case 4: valid_condition = (mcparts_momentum[1] > 0); break;
                default: break;
                }

                if (valid_condition) {
                    E_spectrum_fortime[cell_code]->Fill(Ek, Ek_weight);
                    if (theta > asin(6471. / 6771.)) {
                        E_spectrum_fortime_cut[cell_code]->Fill(Ek, Ek_weight);
                    }
                }
            }
        }
        is.close();
    }
    return true;
}

bool AnalysisAlg::finalize() {
    SniperPtr<RootWriter> ws(getParent(), "wSvc");
    if (ws.valid()) {
        ws->attach("Fkey", E_spectrum);
        ws->attach("Fkey", E_spectrum0);
        ws->attach("Fkey", E_spectrum1);
        ws->attach("Fkey", E_spectrum2);
        ws->attach("Fkey", E_spectrum3);
        ws->attach("Fkey", E_spectrum4);
        ws->attach("Fkey", E_spectrum_fortime);
        ws->attach("Fkey", E_spectrum_fortime0);
        ws->attach("Fkey", E_spectrum_fortime1);
        ws->attach("Fkey", E_spectrum_fortime2);
        ws->attach("Fkey", E_spectrum_fortime3);
        ws->attach("Fkey", E_spectrum_fortime4);
        ws->attach("Fkey", E_spectrum_cut);
        ws->attach("Fkey", E_spectrum0_cut);
        ws->attach("Fkey", E_spectrum1_cut);
        ws->attach("Fkey", E_spectrum2_cut);
        ws->attach("Fkey", E_spectrum3_cut);
        ws->attach("Fkey", E_spectrum4_cut);
        ws->attach("Fkey", E_spectrum_fortime_cut);
        ws->attach("Fkey", E_spectrum_fortime0_cut);
        ws->attach("Fkey", E_spectrum_fortime1_cut);
        ws->attach("Fkey", E_spectrum_fortime2_cut);
        ws->attach("Fkey", E_spectrum_fortime3_cut);
        ws->attach("Fkey", E_spectrum_fortime4_cut);
        return true;
    } else {
        LogError << "Failed to attach histogram" << std::endl;
        return false;
    }
    LogInfo << " finalized successfully" << std::endl;
    return true;
}

double AnalysisAlg::flux_ams(double x) {
    double electron = pow(x, 2) / pow(x + 0.87, 2) * pow(1 + pow((x + 0.87) / 3.94, -2.14), -1) *
                      (1.13e-2 * pow((x + 0.87) / 20, -4.31) + 3.96e-6 * pow((x + 0.87) / 300, -3.14));
    double positron =
        pow(x, 2) / pow(x + 1.1, 2) *
        (6.51e-2 * pow((x + 1.1) / 7, -4.07) + 6.80e-5 * pow((x + 1.1) / 60, -2.58) * TMath::Exp(-(x + 1.1) / 810));
    return electron + positron;
}
double AnalysisAlg::flux_dampe(double x) {
    return 1.62e-4 * pow((x / 100), -3.09) * pow(1 + pow((x / 914), -(3.09 - 3.92) / 0.1), -0.1);
}

void AnalysisAlg::vec_in_detector_to_json(double vx, double vy, double vz, double mjd, double *dec_pos, double *dec_vel,
                                          double lat, double lon, double &theta, double &arc) {
    vx = -vx; // v in detector,reverse
    vy = -vy;
    vz = -vz;
    double x_pointing[3], y_pointing[3], z_pointing[3]; // detector xyz in eci
    for (int i = 0; i < 3; i++) {
        x_pointing[i] = dec_vel[i] / sqrt(dec_vel[0] * dec_vel[0] + dec_vel[1] * dec_vel[1] + dec_vel[2] * dec_vel[2]);
        z_pointing[i] = dec_pos[i] / sqrt(dec_pos[0] * dec_pos[0] + dec_pos[1] * dec_pos[1] + dec_pos[2] * dec_pos[2]);
    }
    y_pointing[0] = z_pointing[1] * x_pointing[2] - z_pointing[2] * x_pointing[1];
    y_pointing[1] = z_pointing[2] * x_pointing[0] - z_pointing[0] * x_pointing[2];
    y_pointing[2] = z_pointing[0] * x_pointing[1] - z_pointing[1] * x_pointing[0];
    double v_in_eci[3];
    v_in_eci[0] = x_pointing[0] * vx + y_pointing[0] * vy + z_pointing[0] * vz;
    v_in_eci[1] = x_pointing[1] * vx + y_pointing[1] * vy + z_pointing[1] * vz;
    v_in_eci[2] = x_pointing[2] * vx + y_pointing[2] * vy + z_pointing[2] * vz;
    double v_in_gtod[3] = {v_in_eci[0], v_in_eci[1], v_in_eci[2]};
    FT_Equat2GTOD(v_in_gtod[0], v_in_gtod[1], v_in_gtod[2], mjd);
    double v_final_xyz[3], v_final_polar[3];
    rotatexyz(v_in_gtod, -(90 - lat) * pi / 180, -lon * pi / 180, v_final_xyz);
    xyz2polar(v_final_xyz, v_final_polar);
    theta = v_final_polar[1];
    arc = v_final_polar[2];
    if (arc < 0) arc = arc + 2 * pi;
}
double AnalysisAlg::FT_Modulus(double arg1, double arg2) {
    /* Returns arg1 mod arg2 */
    int i;
    double ret_val;
    ret_val = arg1;
    i = (int)(ret_val / arg2);
    ret_val -= i * arg2;
    if (ret_val < 0.0) ret_val += arg2;
    return ret_val;
}
double AnalysisAlg::FT_GMST_rad(double timeUnix) {
    double time = timeUnix + 2400000.5;
    /* Greenwich mean sidereal time  in radians */
    double T = (time - 2451545.0) / 36525.0; /*Interval of time, measured in Julian centuries of
                                                36525 days of UT (mean solar day), elapsed since the
                                                epoch 2000 Jan 1d12hUT*/
    double d = (time - 2451545.0);
    double GMST_deg; /*greenwich mean sidereal time in degree (i.e. the
                        Greenwich hour angle of the mean equinox of date)*/
    /* GMST from Meeus, J., 2000. Astronomical Algorithms. Willman-Bell,
       Richmond,VA, 2nd ed.  p 84 (eq.11-4) adapdet from IDL procedure
       http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro */
    /*file:///C:/Users/19104/Downloads/astronomical-algorithms_compress.pdf*/
    GMST_deg = FT_Modulus((280.46061837 + 360.98564736629 * d + T * T * (0.000387933 - T / 38710000.0)), 360);
    // GMST_deg = (int)(280.46061837 + 360.98564736629 * d + T * T *
    // (0.000387933 -T /38710000.0)) % 360;  // AL - possible variant
    double GMSTrad = GMST_deg * pi / 180.; /* greenwich mean sidereal time  in radians*/
    return GMSTrad;
}
void AnalysisAlg::FT_GTOD2Equat(double &x, double &y, double &z, double time) {
    /* Conversion from GTOD to ECI coordinates */
    double oldX = x;
    double oldY = y;
    double oldZ = z;
    double alpha_g = 0;

    alpha_g = FT_GMST_rad(time);

    x = oldX * cos(alpha_g) - oldY * sin(alpha_g);
    y = oldX * sin(alpha_g) + oldY * cos(alpha_g);
    z = oldZ;
}

void AnalysisAlg::FT_Equat2GTOD(double &x, double &y, double &z, double time) {
    /* Conversion from ECI to GTOD coordinates */
    double oldX = x;
    double oldY = y;
    double oldZ = z;
    double alpha_g = 0;

    alpha_g = FT_GMST_rad(time);

    x = oldX * cos(alpha_g) + oldY * sin(alpha_g);
    y = -oldX * sin(alpha_g) + oldY * cos(alpha_g);
    z = oldZ;
}

void AnalysisAlg::FT_GTOD2Equat(double &x, double &y, double &z, double &vx, double &vy, double &vz, double time) {
    /* Conversion from ECI to GTOD coordinates */
    double oldX = x;
    double oldY = y;
    double oldZ = z;
    double alpha_g = 0;

    alpha_g = FT_GMST_rad(time);

    x = oldX * cos(alpha_g) - oldY * sin(alpha_g);
    y = oldX * sin(alpha_g) + oldY * cos(alpha_g);
    z = oldZ;

    double oldVX = vx;
    double oldVY = vy;
    double oldVZ = vz;
    double Omega = 2 * pi / 86400.; /* a complete round during a mean solar day*/

    vx = oldVX * cos(alpha_g) - oldVY * sin(alpha_g) + Omega * (-oldY * cos(alpha_g) - oldX * sin(alpha_g));
    vy = oldVX * sin(alpha_g) + oldVY * cos(alpha_g) + Omega * (oldX * cos(alpha_g) - oldY * sin(alpha_g));
    vz = oldVZ;
}

void AnalysisAlg::FT_Equat2GTOD(double &x, double &y, double &z, double &vx, double &vy, double &vz, double time) {
    /* Conversion from ECI to GTOD coordinates */
    double oldX = x;
    double oldY = y;
    double oldZ = z;
    double alpha_g = 0;

    alpha_g = FT_GMST_rad(time);

    x = oldX * cos(alpha_g) + oldY * sin(alpha_g);
    y = -oldX * sin(alpha_g) + oldY * cos(alpha_g);
    z = oldZ;

    /* to transform the velocity
     *
     *        |cos(alpha_g)   sin(alpha_g)   0  |    |vx|           | 0  1  0 |
     * |cos(alpha_g)   sin(alpha_g)   0  | |x| V_GTOD=|-sin(alpha_g)
     * cos(alpha_g) 0  |  * |vy| + Omega * |-1  0  0 |*|-sin(alpha_g)
     * cos(alpha_g)   0  |*|y| |    0              0          1  |    |vz| | 0
     * 0  0 | | 0              0          1  | |z|
     *
     *                                                              |y*cos(alpha_g)-x*sin(alpha_g)
     * |
     *                   ....                             + Omega *
     * |-x*cos(alpha_g)-y*sin(alpha_g)| |             0                |
     */
    double oldVX = vx;
    double oldVY = vy;
    double oldVZ = vz;
    double Omega = 2 * pi / 86400.; /* a complete round during a mean solar day*/

    vx = oldVX * cos(alpha_g) + oldVY * sin(alpha_g) + Omega * (oldY * cos(alpha_g) - oldX * sin(alpha_g));
    vy = -oldVX * sin(alpha_g) + oldVY * cos(alpha_g) + Omega * (-oldX * cos(alpha_g) - oldY * sin(alpha_g));
    vz = oldVZ;
}
void AnalysisAlg::xyz2polar(double *xyz, double *r_theta_phi) {
    r_theta_phi[0] = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
    r_theta_phi[1] = acos(xyz[2] / r_theta_phi[0]);
    r_theta_phi[2] = atan2(xyz[1], xyz[0]);
}
void AnalysisAlg::polar2xyz(double *xyz, double *r_theta_phi) {
    xyz[0] = r_theta_phi[0] * sin(r_theta_phi[1]) * cos(r_theta_phi[2]);
    xyz[1] = r_theta_phi[0] * sin(r_theta_phi[1]) * sin(r_theta_phi[2]);
    xyz[2] = r_theta_phi[0] * cos(r_theta_phi[1]);
}
void AnalysisAlg::rotatexyz(double *xyz, double theta, double phi, double *xyz_new) {
    // rotate along z axis for phi first, then rotate along y axis for theta
    // from gtod to velocity
    xyz_new[0] = xyz[0] * cos(theta) * cos(phi) - xyz[1] * sin(phi) * cos(theta) + xyz[2] * sin(theta);
    xyz_new[1] = xyz[0] * sin(phi) + xyz[1] * cos(phi);
    xyz_new[2] = -xyz[0] * sin(theta) * cos(phi) + xyz[1] * sin(theta) * sin(phi) + xyz[2] * cos(theta);
}
