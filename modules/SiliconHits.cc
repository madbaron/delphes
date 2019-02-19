/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class SiliconHits
 *
 *  Finds hits in idealised tracking detector
 *
 *  \author K. Albertsson - CERN LTU/EISLAB
 *  \author F. Meloni - DESY
 *
 */

#include "modules/SiliconHits.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

using point_t = std::tuple<double, double, double>;
const double c_light = 2.99792458E8;

struct Helix {
  Helix(point_t centre_xyz, double R, double omega, double phi_0, double pz,
        double gammam) : centre_xyz(centre_xyz), R(R), omega(omega),
    phi_0(phi_0), pz(pz), gammam(gammam) {}
  Helix(const Candidate * candidate, double Bz) {
    TLorentzVector candidatePosition = candidate->InitialPosition;
    TLorentzVector candidateMomentum = candidate->Momentum;

    q = candidate->Charge;
    pt = candidateMomentum.Pt();
    pz = candidateMomentum.Pz();
    phi_0 = TMath::ATan2(candidateMomentum.Py(), candidateMomentum.Px() ); // [rad] in [-pi, pi]

    gammam = candidateMomentum.E() * 1.0E9 / (c_light * c_light);  // gammam in [eV/c^2]
    omega = candidate->Charge * Bz / (gammam);  // omega is here in [89875518/s]
    R = pt / (q * Bz) * 1.0E9 / c_light;      // in [m]

    // helix initial coordinates
    double x_c = candidatePosition.X() + R * TMath::Sin(phi_0) * 1.0E3;
    double y_c = candidatePosition.Y() - R * TMath::Cos(phi_0) * 1.0E3;
    double z_c = candidatePosition.Z();

    centre_xyz = {x_c, y_c, z_c};
    start_xyz = {candidatePosition.X(),
                 candidatePosition.Y(),
                 candidatePosition.Z()
                };
    start_rphiz = {TMath::Hypot(candidatePosition.X(), candidatePosition.Y()),
                   TMath::ATan2(candidatePosition.Y(), candidatePosition.X()),
                   candidatePosition.Z()
                  };
  }

  point_t operator()(double t) const {
    return xyz(t);
  }

  point_t xyz(double t) const {
    double x_t = std::get<0>(centre_xyz) + R * TMath::Sin(omega * t - phi_0) * 1.0E3; // in [mm]
    double y_t = std::get<1>(centre_xyz) + R * TMath::Cos(omega * t - phi_0) * 1.0E3; // in [mm]
    double z_t = std::get<2>(centre_xyz) + pz * 1.0E9 / c_light / gammam * t * 1.0E3; // in [mm]
    return {x_t, y_t, z_t};
  }

  double x(double t) const {
    return std::get<0>(centre_xyz) + R * TMath::Sin(omega * t - phi_0) * 1.0E3; // in [mm]
  }

  double y(double t) const {
    return std::get<1>(centre_xyz) + R * TMath::Cos(omega * t - phi_0) * 1.0E3; // in [mm]
  }

  double z(double t) const {
    return std::get<2>(centre_xyz) + pz * 1.0E9 / c_light / gammam * t * 1.0E3; // in [mm]
  }

  point_t rphiz(double t) const {
    return {r(t), phi(t), z(t)};
  }

  double r(double t) const {
    return TMath::Hypot(x(t), y(t));
  }

  double phi(double t) const {
    return TMath::ATan2(y(t), x(t));
  }

public:
  // Centre of helix
  point_t centre_xyz;
  point_t start_xyz;
  point_t start_rphiz;

  double R;
  double omega;
  double phi_0;
  double q;
  double pt;
  double pz;

  double gammam;
};

double bisect_barrel(const Helix & helix, double t_start, double t_end,
                     double barrel_r, double precision = 1e-6, size_t max_iter = 100) {

  if (t_start > t_end) {
    std::swap(t_start, t_end);
  }

  auto is_before = [](double r, double barrel_r) {return r < barrel_r;};
  auto is_after = [](double r, double barrel_r) {return r >= barrel_r;};
  auto is_close = [precision](double r, double barrel_r) {
    return fabs(r - barrel_r) < precision;
  };

  for (size_t iIter = 0; iIter < max_iter; ++iIter) {
    double t = t_start + (t_end - t_start) / 2.0;
    double r = helix.r(t);

    // if (fDebug > 1) {
    //   std::cout << "iter: " << iIter
    //             << " t: " << t
    //             << " r_start:" << helix.r(t_start)
    //             << " r_end:" << helix.r(t_end)
    //             << " r: " << r
    //             << " barrel_r: " << barrel_r
    //             << " diff: " << fabs(r - barrel_r)
    //             << std::endl;
    // }

    if (is_close(r, barrel_r)) {return t;}
    else if (is_before(r, barrel_r)) {t_start = t;}
    else if (is_after(r, barrel_r)) {t_end = t;}
  }

  return t_start + (t_end - t_start) / 2.0;
}

double bisect_disk(const Helix & helix, double t_start, double t_end,
                   double disk_z, double precision = 1e-6, size_t max_iter = 100) {
  if (t_start > t_end) {
    std::swap(t_start, t_end);
  }

  auto is_before = [](double z, double disk_z) {return z < disk_z;};
  auto is_after = [](double z, double disk_z) {return z >= disk_z;};
  auto is_close = [precision](double z, double disk_z) {
    return fabs(z - disk_z) < precision;
  };

  for (size_t iIter = 0; iIter < max_iter; ++iIter) {
    double t = t_start + (t_end - t_start) / 2.0;
    double z = helix.z(t);

    if (is_close(z, disk_z)) {return t;}
    else if (is_before(z, disk_z)) {t_start = t;}
    else if (is_after(z, disk_z)) {t_end = t;}
  }

  return t_start + (t_end - t_start) / 2.0;
}


//------------------------------------------------------------------------------

SiliconHits::SiliconHits() : fDebug(0), fItInputArray(0)
{
  std::vector<double> pix_barrel_r {33.4, 50.5, 88.5, 122.5};
  std::vector<double> sct_barrel_r {299., 371., 443., 514.};

  std::vector<double> pix_barrel_z {330.15, 401., 401., 401.};
  std::vector<double> sct_barrel_z {746.00, 746., 746., 746.};

  std::vector<double> pix_disk_rmin = {89., 89., 89.,
                                       89., 89., 89.
                                      };
  std::vector<double> pix_disk_rmax = {150., 150., 150.,
                                       150., 150., 150.
                                      };
  std::vector<double> pix_disk_z = { -650., -580., -495.,
                                     +495., +580., +650.
                                   };

  std::vector<double> sct_disk_rmin = {338., 270., 270., 270., 270., 270., 338., 408., 439.,
                                       338., 270., 270., 270., 270., 270., 338., 408., 439.
                                      };
  std::vector<double> sct_disk_rmax = {560., 560., 560., 560., 560., 560., 560., 560., 560.,
                                       560., 560., 560., 560., 560., 560., 560., 560., 560.
                                      };
  std::vector<double> sct_disk_z = { -847.5, -934., -1084., -1262., -1377., -1747., -2072., -2500., -2713.,
                                     +847.5, +934., +1084., +1262., +1377., +1747., +2072., +2500., +2713.
                                   };

  // Add pixel detector planes
  barrel_r.insert(barrel_r.end(), pix_barrel_r.begin(), pix_barrel_r.end());
  barrel_z.insert(barrel_z.end(), pix_barrel_z.begin(), pix_barrel_z.end());
  disk_rmin.insert(disk_rmin.end(), pix_disk_rmin.begin(), pix_disk_rmin.end());
  disk_rmax.insert(disk_rmax.end(), pix_disk_rmax.begin(), pix_disk_rmax.end());
  disk_z.insert(disk_z.end(), pix_disk_z.begin(), pix_disk_z.end());

  // Add sct detector planes
  barrel_r.insert(barrel_r.end(), sct_barrel_r.begin(), sct_barrel_r.end());
  barrel_z.insert(barrel_z.end(), sct_barrel_z.begin(), sct_barrel_z.end());
  disk_rmin.insert(disk_rmin.end(), sct_disk_rmin.begin(), sct_disk_rmin.end());
  disk_rmax.insert(disk_rmax.end(), sct_disk_rmax.begin(), sct_disk_rmax.end());
  disk_z.insert(disk_z.end(), sct_disk_z.begin(), sct_disk_z.end());
}

//------------------------------------------------------------------------------

SiliconHits::~SiliconHits()
{

}

//------------------------------------------------------------------------------

void SiliconHits::Init()
{
  // read parameters

  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "TrackMergerAll/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "hits"));
}

//------------------------------------------------------------------------------

void SiliconHits::Finish()
{
  if (fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------
// Simulates a candidate track through the geometry and adds any intersections
// to the list of output hits.

void SiliconHits::Track(const Candidate * candidate, size_t partIdx)
{
  DelphesFactory * factory = GetFactory();

  // initial transverse momentum p_{T0}: Part->pt
  // initial transverse momentum direction phi_0 = -atan(p_X0/p_Y0)
  // relativistic gamma: gamma = E/mc^2; gammam = gamma * m
  // gyration frequency omega = q/(gamma m) Bz
  // helix radius r = p_{T0} / (omega gamma m)

  const double Bz = 2.0; // in Tesla

  const Helix helix {candidate, Bz};

  double cumul_distance = 0;
  const double stop_criterion = 1E4;  // in [mm]
  const double hit_precision = 1E0;  // in [mm]
  double step_size = hit_precision * 1.0E-3 / c_light; // in [s]

  if (fDebug > 0) {
    point_t p = helix.start_xyz;
    std::cout << std::endl;
    std::cout << "Origin at " << std::get<0>(p) << " " << std::get<1>(p) << " " << std::get<2>(p) << std::endl;
    std::cout << "pt " << helix.pt << " pz " << helix.pz << std::endl;
    std::cout << "R " << helix.R << std::endl;
  }

  for (size_t iStep = 0; true; ++iStep) {

    double t = iStep * step_size;
    double prev_t = (iStep == 0) ? (0) : (iStep - 1) * step_size;
    cumul_distance += step_size * c_light;

    // compute position
    double r_t = helix.r(t);
    double z_t = helix.z(t);
    double prev_r_t = helix.r(prev_t);
    double prev_z_t = helix.z(prev_t);

    // Interactions with barrel
    for (UInt_t nlayer = 0; nlayer < barrel_r.size(); nlayer++) {

      if (fabs(z_t) > barrel_z[nlayer]) continue;

      double barrel_rval = barrel_r[nlayer];

      if (fDebug > 2) {
        std::cout << std::endl;
        std::cout << "Barrel check: " << barrel_rval
                  << " r : " << r_t
                  << " pr: " << prev_r_t
                  << std::endl;
      }

      bool prev_was_before = (prev_r_t < barrel_rval);
      bool curr_is_after = (r_t >= barrel_rval);
      if (prev_was_before and curr_is_after) {

        double t_exact = bisect_barrel(helix, prev_t, t, barrel_rval);
        // TODO: Verify that "exact" hit is valid i.e.
        //       z_t close to barrel_z

        Candidate * new_hit = factory->NewCandidate();
        new_hit->Xd = helix.r(t_exact);
        new_hit->Yd = helix.z(t_exact);
        new_hit->Zd = helix.phi(t_exact);
        new_hit->partIdx = partIdx;
        new_hit->vxTruth = candidate->vxTruth;
        fOutputArray->Add(new_hit);

        cumul_distance = 0;
      }
    }

    // Interactions with end cap discs
    for (UInt_t ndisk = 0; ndisk < disk_z.size(); ndisk++) {

      if (r_t < disk_rmin[ndisk] or r_t > disk_rmax[ndisk]) continue;

      double disk_zval = disk_z[ndisk];

      bool is_same_sign = (signbit(prev_z_t) == signbit(disk_zval) and
                           signbit(z_t) == signbit(prev_z_t));
      bool prev_was_before = (fabs(prev_z_t) < fabs(disk_zval));
      bool curr_is_after = (fabs(z_t) >= fabs(disk_zval));
      if (prev_was_before and curr_is_after and is_same_sign) {

        double t_exact = bisect_disk(helix, prev_t, t, disk_zval);
        // TODO: Verify that "exact" hit is valid i.e.
        //       r_t close to disk_r

        Candidate * new_hit = factory->NewCandidate();
        new_hit->Xd = helix.r(t_exact);
        new_hit->Yd = helix.z(t_exact);
        new_hit->Zd = helix.phi(t_exact);
        new_hit->partIdx = partIdx;
        new_hit->vxTruth = candidate->vxTruth;
        fOutputArray->Add(new_hit);

        cumul_distance = 0;
      }
    }

    if (fabs(z_t) > 3000. or fabs(r_t) > 1000.0) {
      if (fDebug > 0) {
        std::cout << "Particle leaves tracking volume @ r=" << r_t << " z=" << z_t << std::endl;
        std::cout << "Tracking for particle " << partIdx << " done at step " << iStep << std::endl;
      }
      break;
    }

    if (cumul_distance > stop_criterion) {
      if (fDebug > 0) {
        std::cout << "Max track length reached at r=" << r_t << " z=" << z_t << std::endl;
        std::cout << "Tracking for particle " << partIdx << " done at step " << iStep << std::endl;
      }
      break;
    }
  }
}

//------------------------------------------------------------------------------

void SiliconHits::Process()
{
  // loop over all input candidates
  fItInputArray->Reset();

  auto start = std::chrono::high_resolution_clock::now();

  size_t num_particles = (size_t)fInputArray->GetEntriesFast();
  for (size_t partIdx = 0; partIdx < num_particles; ++partIdx)
  {
    Track(static_cast<Candidate *>(fInputArray->operator[](partIdx)), partIdx);
    if (fDebug > 1) {
      std::cout << "Number of hits recorded so far: " << fOutputArray->GetEntriesFast() << std::endl;
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto secs = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

  std::cout << "Simulated " << num_particles << " particles in " << secs.count() * 1000. << " milliseconds (" << 1.*num_particles / secs.count() << " particles/s)" << std::endl;
}

//------------------------------------------------------------------------------
