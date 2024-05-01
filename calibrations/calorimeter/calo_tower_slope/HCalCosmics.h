// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HCALCOSMICS_H
#define HCALCOSMICS_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

// Forward declarations
class PHCompositeNode;
class TFile;
class TF1;
class TH1;
class TH2;

class HCalCosmics : public SubsysReco
{
 public:
  //! constructor
  HCalCosmics(const std::string &, const std::string &);
  //! destructor
  ~HCalCosmics() override = default;

  //! Processing
  int Init(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int process_towers(PHCompositeNode *);
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *) override;

  void set_tower_threshold(float fac) { tower_threshold = fac; }
  void set_vert_threshold(float fac) { vert_threshold = fac; }
  void set_veto_threshold(float fac) { veto_threshold = fac; }

  void Detector(const std::string &name) { detector = name; }
  void TowerPrefix(const std::string &name) { prefix = name; }

  static double gamma_function(double *x, double *par);
  void fitChannels(const std::string &infile, const std::string &outfile2);
  TF1 *fitHist(TH1 *);

 private:
  // HCal geometry
  static const int n_etabin{24};
  static const int n_phibin{64};

  TFile *outfile{nullptr};
  TH1 *h_channel_hist[n_etabin][n_phibin]{{nullptr}};
  TH2 *h_waveformchi2{nullptr};
  TH2 *h_time_energy{nullptr};
  TH1 *h_mip{nullptr};

  // Cut threshold
  int tower_threshold{500};
  int vert_threshold{500};
  int veto_threshold{350};
  int event{0};

  float m_peak[n_etabin][n_phibin]{};
  float m_chi2[n_etabin][n_phibin]{};

  //  bool debug {false};

  std::string prefix{"TOWERS_"};
  std::string detector{"HCALIN"};
  std::string outfilename;
};

#endif  // HCALCOSMICS_H
