// Discribtion:This code is used to add Fermimotion p_F to spectator neutrons

// include the header file here
#include "FermiMotion.h"

#include <phool/phool.h>

#include <gsl/gsl_rng.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>   // for GenParticle
#include <HepMC/GenVertex.h>     // for GenVertex, GenVertex::part...
#include <HepMC/HeavyIon.h>      // for HeavyIon
#include <HepMC/SimpleVector.h>  // for FourVector

#include <CLHEP/Vector/LorentzVector.h>

#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>

//____________________________________________________________________________..

// this method is use to find out the spectator neutron loss prob
// using the parameterization in the PHENIX Glauber
// Monte Carlo code written by Klaus Reygers to model
// the loss of forward
// neutrons into the ZDC due to larger fragments

namespace
{
  // Assume Au for now
  // make sure b is in fm
  double ploss(double b)
  {
    // para
    double p0 = 0.3305;
    double p1 = 0.0127;
    double p2 = 17.;
    double p3 = 2.;
    double ploss = p0 + (b * p1) + exp((b - p2) / p3);

    return ploss;
  }

  // this method is use to generate and random p_F
  // along a random direction and add it to the momentum
  // assume Au for now
  CLHEP::HepLorentzVector pwithpF(CLHEP::HepLorentzVector p, gsl_rng *RandomGenerator, int id, double pTspec, double bphi)
  {
    // id should be either 2112 or 2212
    if (!((id == 2112) || (id == 2212)))
    {
      std::cout << "wrong pid" << std::endl;
      return p;
    }
    // find pF max using Thomas-Fermi model, assume using Au.
    double pFmax = 0.28315;
    if (id == 2212)
    {
      pFmax = 0.23276;
    }
    // now generate the random p assuming probability is propotional to p^2dp
    // CLHEP::RandGeneral seems to be a better way to do it
    double pF = pFmax * pow(gsl_rng_uniform_pos(RandomGenerator), 1.0 / 3.0);
    double cotheta = (gsl_rng_uniform_pos(RandomGenerator) - 0.5) * 2;
    double phi = gsl_rng_uniform_pos(RandomGenerator) * 2 * M_PI;
    double pFx = pF * sqrt(1 - (cotheta * cotheta)) * cos(phi);
    double pFy = pF * sqrt(1 - (cotheta * cotheta)) * sin(phi);
    double pFz = pF * cotheta;
    double pSx = pTspec * cos(bphi);
    double pSy = pTspec * sin(bphi);

    if (p.pz() < 0)
    {
      pSx *= -1;
      pSy *= -1;
    }

    // now add the pF to p
    double px = p.px() + pFx + pSx;
    double py = p.py() + pFy + pSy;
    double pz = p.pz() + pFz;
    // calculate the total energy
    double const nrm = 0.938;
    double e = sqrt((px * px) + (py * py) + (pz * pz) + (nrm * nrm));

    CLHEP::HepLorentzVector pwithpF(px, py, pz, e);
    return pwithpF;
  }
}  // namespace
int FermiMotion(HepMC::GenEvent *event, gsl_rng *RandomGenerator, double pTspec)
{
  // find ploss
  // std::cout<<"getting b"<<std::endl;
  HepMC::HeavyIon *hi = event->heavy_ion();
  if (!hi)
  {
    std::cout << PHWHERE << ": Fermi Motion Afterburner needs the Heavy Ion Event Info, GenEvent::heavy_ion() returns NULL" << std::endl;
    exit(1);
  }
  double b = hi->impact_parameter();
  double bphi = hi->event_plane_angle();
  double pnl = ploss(b);
  // now loop over all particles and find spectator neutrons

  for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(), prev = event->particles_end(); p != event->particles_end(); prev = p, ++p)
  {
    // if not final state continue
    if ((*p)->status() != 1)
    {
      continue;
    }
    int id = (*p)->pdg_id();
    // if not neutron, skip
    if (!((id == 2112) || (id == 2212)))
    {
      continue;
    }

    // spectator neutron should have px==0&&py==0
    HepMC::GenParticle *n = (*p);
    double p_x = n->momentum().px();
    double p_y = n->momentum().py();
    if (!(p_x == 0 && p_y == 0))
    {
      continue;
    }

    if (id == 2112)
    {
      // std::cout<<"after: "<<n->barcode()<<std::endl;
      if (pnl > gsl_rng_uniform_pos(RandomGenerator))
      {
        // remove particle here
        delete ((*p)->production_vertex())->remove_particle(*p);
        // std::cout<<"removing: "<<n->barcode()<<std::endl;
        p = prev;
        continue;
      }
    }

    // add pF to the remaining

    CLHEP::HepLorentzVector p0(n->momentum().px(), n->momentum().py(), n->momentum().pz(), n->momentum().e());
    CLHEP::HepLorentzVector newp = pwithpF(p0, RandomGenerator, id, pTspec, bphi);
    (*p)->set_momentum(newp);
  }

  return 0;
}
