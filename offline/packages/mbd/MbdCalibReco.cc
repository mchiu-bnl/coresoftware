#include "MbdCalibReco.h"
#include "MbdCalib.h"
#include "MbdDefs.h"
//#include "MbdRawContainer.h"
//#include "MbdRawHit.h"
#include "MbdPmtContainer.h"
#include "MbdPmtHit.h"
#include "MbdOut.h"

#include <ffamodules/CDBInterface.h>
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>
#include <ffarawobjects/Gl1Packet.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

MbdCalibReco::MbdCalibReco(const std::string &name)
  : SubsysReco(name)
{
}

int MbdCalibReco::Init(PHCompositeNode * /*topNode*/)
{
  _mbdcal = new MbdCalib();
  _mbdcal->Verbosity( Verbosity() );
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::InitRun(PHCompositeNode *topNode)
{
  _runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!_runheader)
  {
    std::cout << PHWHERE << " RunHeader node not found, will use run number 0" << std::endl;
  }

  _runnumber = _runheader ? _runheader->get_RunNumber() : 0;

  // Build run directory path and create it
  std::ostringstream oss;
  oss << _caldir << "/" << _runnumber;
  _rundir = oss.str();
  gSystem->Exec(("mkdir -p " + _rundir).c_str());

  if (!_cdbtag.empty())
  {
    // Download baseline calibrations from CDB
    recoConsts::instance()->set_StringFlag("CDB_GLOBALTAG", _cdbtag);
    CDBInterface* cdb = CDBInterface::instance();
    std::string url;

    url = cdb->getUrl("MBD_SAMPMAX");
    if (!url.empty()) { _mbdcal->Download_SampMax(url); }

    url = cdb->getUrl("MBD_PED");
    if (!url.empty()) { _mbdcal->Download_Ped(url); }

    url = cdb->getUrl("MBD_TIMECORR");
    if (!url.empty()) { _mbdcal->Download_TimeCorr(url); }

    url = cdb->getUrl("MBD_SLEWCORR");
    if (!url.empty()) { _mbdcal->Download_SlewCorr(url); }

    std::cout << Name() << ": loaded calibrations from CDB tag " << _cdbtag << std::endl;
  }
  else
  {
    // Load baseline calibrations from local files if they exist
    std::string calfile = _rundir + "/mbd_sampmax.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_SampMax(calfile);
    }
    calfile = _rundir + "/mbd_ped.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_Ped(calfile);
    }

    // Load slew correction for subpass >= 2
    if (_subpass >= 2)
    {
      calfile = _rundir + "/mbd_slewcorr.calib";
      if (gSystem->AccessPathName(calfile.c_str()) == 0)
      {
        _mbdcal->Download_SlewCorr(calfile);
        std::cout << Name() << ": loaded " << calfile << std::endl;
      }
      else
      {
        std::cout << Name() << ": WARNING: " << calfile << " not found" << std::endl;
      }
    }
  }

  // Load t0 offsets for subpass >= 1 (always from local files — outputs of previous subpass)
  if (_subpass >= 1)
  {
    std::string prevpass = "pass" + std::to_string(_subpass - 1) + "_";

    std::string calfile = _rundir + "/" + prevpass + "mbd_tq_t0.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_TQT0(calfile);
      std::cout << Name() << ": loaded " << calfile << std::endl;
    }
    else
    {
      std::cout << Name() << ": WARNING: " << calfile << " not found" << std::endl;
    }

    calfile = _rundir + "/" + prevpass + "mbd_tt_t0.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_TTT0(calfile);
      std::cout << Name() << ": loaded " << calfile << std::endl;
    }
    else
    {
      std::cout << Name() << ": WARNING: " << calfile << " not found" << std::endl;
    }
  }

  // Build bitmask of scaled triggers whose names begin with "MBD N&S"
  _mbias_trigger_mask = 0xfc00;

  InitHistos();

  // Open output ROOT file
  std::string outfname = _rundir + "/calmbdpass2." + std::to_string(_subpass);
  if (_subpass == 0)
  {
    outfname += "_time-" + std::to_string(_runnumber) + ".root";
  }
  else if (_subpass == 1 || _subpass == 2)
  {
    outfname += "_slew-" + std::to_string(_runnumber) + ".root";
  }
  else
  {
    outfname += "_q-" + std::to_string(_runnumber) + ".root";
  }
  _outfile = std::make_unique<TFile>(outfname.c_str(), "RECREATE");
  std::cout << Name() << ": output file " << outfname << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::getNodes(PHCompositeNode *topNode)
{
  _evtheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!_evtheader)
  {
    std::cout << PHWHERE << " EvtHeader not found, will use run number 0" << std::endl;
  }

  _gl1packet = findNode::getClass<Gl1Packet>(topNode,14001);
  if (!_gl1packet)
  {
    _gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " GL1Packet not found" << std::endl;
    }
  }

  /*
  _mbdraws = findNode::getClass<MbdRawContainer>(topNode, "MbdRawContainer");
  if (!_mbdraws)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdRawContainer not found" << std::endl;
    }
  }
  */

  _mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!_mbdpmts)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdPmtContainer not found" << std::endl;
    }
  }

  _mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  if (!_mbdout)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdOut not found" << std::endl;
    }
  }

  _mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
  if (!_mbdgeom)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdGeom not found" << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void MbdCalibReco::InitHistos()
{
  // Histograms must not be associated with the output TFile at creation
  // time (InitRun happens before _outfile is opened above, but we call
  // InitHistos before opening the file, so ROOT's current directory is
  // gROOT or whichever file is current from the framework).
  gROOT->cd();

  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    std::string sn = std::to_string(ipmt);

    h_tt[ipmt] = new TH1F(("h_tt" + sn).c_str(), ("tt" + sn).c_str(), 7000, -30., 30.);
    h_tt[ipmt]->SetXTitle("ns");

    h_tq[ipmt] = new TH1F(("h_tq" + sn).c_str(), ("tq" + sn).c_str(), 7000, -150., 31. * 17.7623);
    h_tq[ipmt]->SetXTitle("ns");

    h_qp[ipmt] = new TH1F(("h_q" + sn).c_str(), ("q" + sn).c_str(), 3000, -100., 14900.);
    h_qp[ipmt]->SetXTitle("ADC");

    if (_subpass >= 1)
    {
      h2_slew[ipmt] = new TH2F(("h2_slew" + sn).c_str(), ("slew curve, ch " + sn).c_str(), 4000, -0.5, 16000. - 0.5, 1100, -5., 6.);
      h2_slew[ipmt]->SetXTitle("ADC");
      h2_slew[ipmt]->SetYTitle("#Delta T (ns)");
    }
    else
    {
      h2_slew[ipmt] = nullptr;
    }
  }

  h2_tt = new TH2F("h2_tt", "ch vs tt", 900, -150., 150., MbdDefs::MBD_N_PMT, -0.5, MbdDefs::MBD_N_PMT - 0.5);
  h2_tt->SetXTitle("tt [ns]");
  h2_tt->SetYTitle("pmt ch");

  h2_tq = new TH2F("h2_tq", "ch vs tq", 900, -150., 150., MbdDefs::MBD_N_PMT, -0.5, MbdDefs::MBD_N_PMT - 0.5);
  h2_tq->SetXTitle("tq [ns]");
  h2_tq->SetYTitle("pmt ch");
}

int MbdCalibReco::process_event(PHCompositeNode *topNode)
{
  getNodes(topNode);

  // Require a scaled "MBD N&S" trigger
  if (_mbias_trigger_mask != 0)
  {
    uint64_t strig = _gl1packet->getScaledVector();  // scaled trigger only
    if ( (strig&_mbias_trigger_mask)==0 )
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Per-event arrays for corrected times
  /*
  std::array<float, MbdDefs::MBD_N_PMT> ttcorr{};
  std::array<float, MbdDefs::MBD_N_PMT> tqcorr{};
  std::array<float, MbdDefs::MBD_N_PMT> adc_arr{};
  ttcorr.fill(std::numeric_limits<float>::quiet_NaN());
  tqcorr.fill(std::numeric_limits<float>::quiet_NaN());
  adc_arr.fill(0);
  */

  std::array<Float_t, MbdDefs::MBD_N_ARMS> armtime{};
  armtime.fill(0);
  std::array<Float_t, MbdDefs::MBD_N_ARMS> nhit{};
  nhit.fill(0);

  Float_t zvtx = _mbdout->get_zvtx();
  // Vertex cut for subpass >= 1
  if ( _subpass >= 1 )
  {
    if (std::abs(zvtx) > 60.)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  for (int iarm=0; iarm<2; iarm++)
  {
    armtime[iarm] = _mbdout->get_time(iarm);
    nhit[iarm] = _mbdout->get_npmt(iarm);
  }

  for (int ipmt=0; ipmt < _mbdpmts->get_npmt(); ipmt++)
  {
    MbdPmtHit *pmt = _mbdpmts->get_pmt(ipmt);
    if ( !pmt )
    {
      continue;
    }

    Short_t pmtno = pmt->get_pmt();
    Float_t q = pmt->get_q();
    Float_t tt  = pmt->get_tt();
    Float_t tq  = pmt->get_tq();

    h_tt[pmtno]->Fill( tt );
    h2_tt->Fill( tt, pmtno );
    h_tq[pmtno]->Fill( tq );
    h2_tq->Fill( tq, pmtno );

    // Fill charge histogram for in-time hits
    if ( std::abs(tt)<26.0 && q > 0.)
    {
      h_qp[pmtno]->Fill( q );
    }

    int arm = _mbdgeom->get_arm( pmtno );

    // Fill slew histogram for subpass >= 1
    if (_subpass >= 1 && h2_slew[pmtno])
    {
      if (nhit[arm] >= 2. && q > 0.)
      {
        float dt = tt - armtime[arm];
        h2_slew[pmtno]->Fill(q, dt);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::End(PHCompositeNode * /*topNode*/)
{
  if (!_outfile)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Write histograms to output file
  _outfile->cd();
  h2_tt->Write();
  h2_tq->Write();
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    h_tt[ipmt]->Write();
    h_tq[ipmt]->Write();
    h_qp[ipmt]->Write();
    if (h2_slew[ipmt])
    {
      h2_slew[ipmt]->Write();
    }
  }

  // Always fit and write t0 (done at every subpass from the accumulated histograms)
  //FitAndWriteT0();

  /*
  if (_subpass == 1 || _subpass == 2)
  {
    FitAndWriteSlew();
  }
  */

  _outfile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::getRunType() const
{
  // Run number → collision system (mirrors get_runtype() from get_runstr.h)
  if (_runnumber <= 30000)
  {
    return 3;  // SIMAUAU200
  }
  if (_runnumber <= 53880)
  {
    return 1;  // PP200 (Run24)
  }
  if (_runnumber <= 54962)
  {
    return 0;  // AUAU200 (Run24)
  }
  if (_runnumber <= 78954)
  {
    return 0;  // AUAU200 (Run25)
  }
  if (_runnumber <= 81667)
  {
    return 1;  // PP200 (Run25)
  }
  if (_runnumber <= 82703)
  {
    return 2;  // OO200 (Run25)
  }
  return -1;
}

// ---------------------------------------------------------------------------
// FitAndWriteT0 — Gaussian fit to h_tt and h_tq, write *_t0.calib files
// ---------------------------------------------------------------------------
void MbdCalibReco::FitAndWriteT0()
{
  std::string passprefix = "pass" + std::to_string(_subpass) + "_";
  std::string tt_fname = _rundir + "/" + passprefix + "mbd_tt_t0.calib";
  std::string tq_fname = _rundir + "/" + passprefix + "mbd_tq_t0.calib";

  std::ofstream tt_file(tt_fname);
  std::ofstream tq_file(tq_fname);
  if (!tt_file.is_open() || !tq_file.is_open())
  {
    std::cout << Name() << "::FitAndWriteT0 ERROR cannot open calib files" << std::endl;
    return;
  }

  TF1 gaussian("mbdcal_gaus", "gaus", -25., 25.);
  gaussian.SetLineColor(2);

  double min_twindow = -25.;
  double max_twindow =  25.;

  // --- tt_t0 ---
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    if (ipmt == 0 || ipmt == 64)
    {
      h_tt[ipmt]->SetAxisRange(-25., 25.);
    }
    else
    {
      h_tt[ipmt]->SetAxisRange(min_twindow, max_twindow);
    }

    int    peakbin = h_tt[ipmt]->GetMaximumBin();
    double mean    = h_tt[ipmt]->GetBinCenter(peakbin);
    double peak    = h_tt[ipmt]->GetMaximum();

    gaussian.SetParameters(peak, mean, 5.);
    gaussian.SetRange(mean - 3., mean + 3.);
    h_tt[ipmt]->Fit(&gaussian, "RQ");

    mean              = gaussian.GetParameter(1);
    double meanerr    = gaussian.GetParError(1);
    double sigma      = gaussian.GetParameter(2);
    double sigmaerr   = gaussian.GetParError(2);

    if (ipmt == 0 || ipmt == 64)
    {
      min_twindow = mean - 3. * sigma;
      max_twindow = mean + 3. * sigma;
    }

    tt_file << ipmt << "\t" << mean << "\t" << meanerr << "\t"
            << sigma << "\t" << sigmaerr << "\n";

    // Normalise h2_tt row by fit peak amplitude
    double fitpeak = gaussian.GetParameter(0);
    if (fitpeak != 0.)
    {
      int nbinsx = h2_tt->GetNbinsX();
      for (int ibinx = 1; ibinx <= nbinsx; ibinx++)
      {
        float bc = h2_tt->GetBinContent(ibinx, ipmt + 1);
        h2_tt->SetBinContent(ibinx, ipmt + 1, bc / fitpeak);
      }
    }
  }
  tt_file.close();

  // Write canonical CDB ROOT file
  {
    MbdCalib tmpcal;
    tmpcal.Download_TTT0(tt_fname);
    std::string cdb_fname = tt_fname;
    cdb_fname.replace(cdb_fname.rfind(".calib"), 6, ".root");
    tmpcal.Write_CDB_TTT0(cdb_fname);
  }

  // --- tq_t0 ---
  min_twindow = -25.;
  max_twindow =  25.;

  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    if (ipmt == 0 || ipmt == 64)
    {
      h_tq[ipmt]->SetAxisRange(-25., 25.);
    }
    else
    {
      h_tq[ipmt]->SetAxisRange(min_twindow, max_twindow);
    }

    int    peakbin = h_tq[ipmt]->GetMaximumBin();
    double mean    = h_tq[ipmt]->GetBinCenter(peakbin);
    double peak    = h_tq[ipmt]->GetMaximum();

    gaussian.SetParameters(peak, mean, 5.);
    gaussian.SetRange(mean - 3., mean + 3.);
    h_tq[ipmt]->Fit(&gaussian, "RQ");

    mean              = gaussian.GetParameter(1);
    double meanerr    = gaussian.GetParError(1);
    double sigma      = gaussian.GetParameter(2);
    double sigmaerr   = gaussian.GetParError(2);

    if (ipmt == 0 || ipmt == 64)
    {
      min_twindow = mean - 3. * sigma;
      max_twindow = mean + 3. * sigma;
    }

    tq_file << ipmt << "\t" << mean << "\t" << meanerr << "\t"
            << sigma << "\t" << sigmaerr << "\n";

    // Normalise h2_tq row
    double fitpeak = gaussian.GetParameter(0);
    if (fitpeak != 0.)
    {
      int nbinsx = h2_tq->GetNbinsX();
      for (int ibinx = 1; ibinx <= nbinsx; ibinx++)
      {
        float bc = h2_tq->GetBinContent(ibinx, ipmt + 1);
        h2_tq->SetBinContent(ibinx, ipmt + 1, bc / fitpeak);
      }
    }
  }
  tq_file.close();

  {
    MbdCalib tmpcal;
    tmpcal.Download_TQT0(tq_fname);
    std::string cdb_fname = tq_fname;
    cdb_fname.replace(cdb_fname.rfind(".calib"), 6, ".root");
    tmpcal.Write_CDB_TQT0(cdb_fname);
  }

  std::cout << Name() << ": wrote " << tt_fname << " and " << tq_fname << std::endl;
}

// ---------------------------------------------------------------------------
// FindTH2Ridge — column-by-column Gaussian fits to find slew-correction ridge
// ---------------------------------------------------------------------------
void MbdCalibReco::FindTH2Ridge(const TH2 *h2, TGraphErrors *&gridge,
                                 TGraphErrors *&grms) const
{
  int nbinsx = h2->GetNbinsX();
  double min_yrange = h2->GetYaxis()->GetBinLowEdge(1);
  double max_yrange = h2->GetYaxis()->GetBinLowEdge(h2->GetNbinsY() + 1);

  gridge = new TGraphErrors();
  gridge->SetName("gridge");
  gridge->SetTitle("ridge");
  grms   = new TGraphErrors();
  grms->SetName("grms");
  grms->SetTitle("rms of ridge");

  TH1 *h_projx = h2->ProjectionX("_projx_tmp");
  TF1  gaussian("_slew_gaus", "gaus", min_yrange, max_yrange);
  gaussian.SetLineColor(4);

  TH1  *h_projy = nullptr;
  double adcmean = 0.;
  double adcnum  = 0.;

  for (int ibin = 1; ibin <= nbinsx; ibin++)
  {
    std::string projname = "_hproj_" + std::to_string(ibin);
    if (!h_projy)
    {
      h_projy  = h2->ProjectionY(projname.c_str(), ibin, ibin);
      adcmean  = h_projx->GetBinCenter(ibin);
      adcnum   = 1.;
    }
    else
    {
      TH1 *hadd = h2->ProjectionY(projname.c_str(), ibin, ibin);
      h_projy->Add(hadd);
      delete hadd;
      adcmean += h_projx->GetBinCenter(ibin);
      adcnum  += 1.;
    }

    if (h_projy->Integral() > 2000. || ibin == nbinsx)
    {
      adcmean /= adcnum;

      int    maxbin = h_projy->GetMaximumBin();
      double xmax_g = h_projy->GetBinCenter(maxbin);
      double ymax_g = h_projy->GetBinContent(maxbin);
      gaussian.SetParameter(0, ymax_g);
      gaussian.SetParameter(1, xmax_g);
      gaussian.SetRange(xmax_g - 0.6, xmax_g + 0.6);
      h_projy->Fit(&gaussian, "RWWQ");

      double mean    = gaussian.GetParameter(1);
      double meanerr = gaussian.GetParError(1);
      double rms     = gaussian.GetParameter(2);
      double rmserr  = gaussian.GetParError(2);

      if (meanerr < 1.0)
      {
        int n = gridge->GetN();
        gridge->SetPoint(n, adcmean, mean);
        gridge->SetPointError(n, 0., meanerr);
      }
      if (rmserr < 0.01)
      {
        int n = grms->GetN();
        grms->SetPoint(n, adcmean, rms);
        grms->SetPointError(n, 0., rmserr);
      }

      delete h_projy;
      h_projy  = nullptr;
      adcmean  = 0.;
      adcnum   = 0.;
    }
  }

  gridge->SetBit(TGraph::kIsSortedX);
  grms->SetBit(TGraph::kIsSortedX);
  delete h_projx;
}

// ---------------------------------------------------------------------------
// FitAndWriteSlew — build slew-correction LUT from h2_slew ridge
// ---------------------------------------------------------------------------
void MbdCalibReco::FitAndWriteSlew()
{
  const int NPOINTS = 16000;
  const int MINADC  = 0;
  const int MAXADC  = 15999;

  std::string scorr_fname = _rundir + "/mbd_slewcorr.calib";
  std::ofstream scorr_file(scorr_fname);
  if (!scorr_file.is_open())
  {
    std::cout << Name() << "::FitAndWriteSlew ERROR cannot open " << scorr_fname << std::endl;
    return;
  }

  std::string trms_fname = _rundir + "/mbd_timerms.calib";
  std::ofstream trms_file(trms_fname);

  // Arrays of slew/trms graph pointers, indexed by feech (only T-channels used)
  std::array<TGraphErrors *, MbdDefs::MBD_N_FEECH> g_slew{};
  std::array<TGraphErrors *, MbdDefs::MBD_N_FEECH> g_trms{};
  g_slew.fill(nullptr);
  g_trms.fill(nullptr);

  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    if (!h2_slew[ipmt])
    {
      continue;
    }

    int feech_t = (ipmt / 8) * 16 + ipmt % 8;

    TGraphErrors *gr = nullptr;
    TGraphErrors *grms_tmp = nullptr;
    FindTH2Ridge(h2_slew[ipmt], gr, grms_tmp);

    g_slew[feech_t] = gr;
    g_trms[feech_t] = grms_tmp;

    if (gr)
    {
      gr->SetName(("g_slew" + std::to_string(ipmt)).c_str());
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.25);
    }
    if (grms_tmp)
    {
      grms_tmp->SetName(("g_trms" + std::to_string(ipmt)).c_str());
      grms_tmp->SetMarkerStyle(20);
      grms_tmp->SetMarkerSize(0.25);
    }
  }

  // Write slew correction LUT (one T-channel feech at a time)
  for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
  {
    // Only T-channels: type = (feech/8) % 2 == 0
    if ((ifeech / 8) % 2 == 1)
    {
      continue;
    }

    if (!g_slew[ifeech])
    {
      continue;
    }

    scorr_file << ifeech << "\t" << NPOINTS << "\t" << MINADC << "\t" << MAXADC << "\n";
    int step = (MAXADC - MINADC) / (NPOINTS - 1);
    for (int iadc = MINADC; iadc <= MAXADC; iadc += step)
    {
      scorr_file << g_slew[ifeech]->Eval(iadc) << " ";
      if (iadc % 10 == 9)
      {
        scorr_file << "\n";
      }
    }
  }
  scorr_file.close();

  // Write time-RMS LUT
  if (trms_file.is_open())
  {
    for (int ifeech = 0; ifeech < MbdDefs::MBD_N_FEECH; ifeech++)
    {
      if ((ifeech / 8) % 2 == 1)
      {
        continue;
      }
      if (!g_trms[ifeech])
      {
        continue;
      }

      trms_file << ifeech << "\t" << NPOINTS << "\t" << MINADC << "\t" << MAXADC << "\n";
      int step = (MAXADC - MINADC) / (NPOINTS - 1);
      for (int iadc = MINADC; iadc <= MAXADC; iadc += step)
      {
        trms_file << g_trms[ifeech]->Eval(iadc) << " ";
        if (iadc % 10 == 9)
        {
          trms_file << "\n";
        }
      }
    }
    trms_file.close();
  }

  // Write graphs to ROOT file and create CDB ROOT file
  _outfile->cd();
  for (auto *g : g_slew)
  {
    if (g)
    {
      g->Write();
    }
  }
  for (auto *g : g_trms)
  {
    if (g)
    {
      g->Write();
    }
  }

  {
    MbdCalib tmpcal;
    tmpcal.Download_SlewCorr(scorr_fname);
    std::string cdb_fname = scorr_fname;
    cdb_fname.replace(cdb_fname.rfind(".calib"), 6, ".root");
    tmpcal.Write_CDB_SlewCorr(cdb_fname);
  }

  // Clean up
  for (auto *g : g_slew)
  {
    delete g;
  }
  for (auto *g : g_trms)
  {
    delete g;
  }

  std::cout << Name() << ": wrote " << scorr_fname << std::endl;
}

