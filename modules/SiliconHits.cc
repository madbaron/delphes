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
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
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
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

SiliconHits::SiliconHits() : fItInputArray(0)
{

  barrel_r = {33.4, 50.5, 88.5, 122.5, 299., 371., 443., 514.};
  barrel_z = {330.15, 401., 401., 401., 746., 746., 746., 746.}; //this is half of the layer lenght

  /*
  //Pixel Only
  barrel_r = {33.4, 50.5, 88.5, 122.5};
  barrel_z = {330.15, 400., 400., 400.}; //this is half of the layer lenght
  disk_rmin = 89.;
  disk_rmax = 150.;
  disk_z = {-650., -580., -495., 495., 580., 650.}; //this is the distance from 0.
  */

  disk_rmin = {89., 89., 89., 338., 270., 270., 270., 270., 270., 338., 408., 439., 89., 89., 89., 338., 270., 270., 270., 270., 270., 338., 408., 439.};
  disk_rmax = {150., 150., 150., 560., 560., 560., 560., 560., 560., 560., 560., 560., 150., 150., 150., 560., 560., 560., 560., 560., 560., 560., 560., 560.};
//  disk_z = { -2713., -2500., -2072., -1747., -1377., -1262., -1084., -934., -847.5, -650., -580., -495., 495., 580., 650., 847.5, 934., 1084., 1262., 1377., 1747., 2072., 2500., 2713.};
  disk_z = { -495., -580., -650., -847.5, -934., -1084., -1262., -1377., -1747., -2072., -2500., -2713., 495., 580., 650., 847.5, 934., 1084., 1262., 1377., 1747., 2072., 2500., 2713.};

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

void SiliconHits::Process()
{
  Candidate *candidate, *new_hit;
  TLorentzVector candidatePosition, candidateMomentum;
  DelphesFactory *factory;

  // loop over all input candidates
  fItInputArray->Reset();
  Int_t my_partIdx = 0;

  factory = GetFactory();

  while ((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->InitialPosition;
    candidateMomentum = candidate->Momentum;

    // initial transverse momentum p_{T0}: Part->pt
    // initial transverse momentum direction phi_0 = -atan(p_X0/p_Y0)
    // relativistic gamma: gamma = E/mc^2; gammam = gamma * m
    // gyration frequency omega = q/(gamma m) fBz
    // helix radius r = p_{T0} / (omega gamma m)

    const double c_light = 2.99792458E8;
    const double fBz = 2.0; // in Tesla

    double pt = candidateMomentum.Pt();
    double pz = candidateMomentum.Pz();
    double phi_0 = TMath::ATan2(candidateMomentum.Py(), candidateMomentum.Px() ); // [rad] in [-pi, pi]
    double q = candidate->Charge;

    double gammam = candidateMomentum.E() * 1.0E9 / (c_light * c_light);  // gammam in [eV/c^2]
    double omega = candidate->Charge * fBz / (gammam);  // omega is here in [89875518/s]
    double r = pt / (q * fBz) * 1.0E9 / c_light;      // in [m]

    // helix initial coordinates
    double x_c = candidatePosition.X() * 1.0E-3 + r * TMath::Sin(phi_0);
    double y_c = candidatePosition.Y() * 1.0E-3 - r * TMath::Cos(phi_0);
    double z_c = candidatePosition.Z() * 1.0E-3;

    double cumul_distance = 0;
    double t = 0;
    const double step = 1E-10;
    const double stop_criterion = 3E4;
    const double hit_precision = 1E-1;

    std::cout << "Origin at " << x_c << " " << y_c << " " << z_c << std::endl;
    std::cout << "pt " << pt << " pz " << pz << std::endl;
    std::cout << "R " << r << std::endl;

    while (1) {

      // compute position in terms of x(t), y(t), z(t)
      double x_t = x_c + r * TMath::Sin(omega * t - phi_0) * 1.0E-3;
      double y_t = y_c + r * TMath::Cos(omega * t - phi_0) * 1.0E-3;
      double z_t = z_c + pz * 1.0E9 / c_light / gammam * t * 1.0E-3;
      t += step;
      cumul_distance += step*c_light;

      double pos_r = TMath::Hypot(x_t, y_t);

      for (UInt_t nlayer = 0; nlayer < barrel_r.size(); nlayer++) {

        if(fabs(z_t) > barrel_z[nlayer]) continue;

        if(fabs(pos_r-barrel_r[nlayer]) < hit_precision){

          double hit_r = pos_r;
          double hit_phi = TMath::ATan2(y_t, x_t);
          double hit_z = z_t;

          new_hit = factory->NewCandidate();
          new_hit->Xd = hit_r;
          new_hit->Yd = hit_z;
          new_hit->Zd = hit_phi;
          new_hit->partIdx = my_partIdx;
          new_hit->vxTruth = candidate->vxTruth;
          fOutputArray->Add(new_hit);

          cumul_distance = 0;
        }
      }

      for (UInt_t ndisk = 0; ndisk < disk_z.size(); ndisk++) {

        if(pos_r < disk_rmin[ndisk] || pos_r > disk_rmax[ndisk]) continue;

        if(fabs(z_t-disk_z[ndisk]) < hit_precision){

          double hit_r = pos_r;
          double hit_phi = TMath::ATan2(y_t, x_t);
          double hit_z = z_t;

          new_hit = factory->NewCandidate();
          new_hit->Xd = hit_r;
          new_hit->Yd = hit_z;
          new_hit->Zd = hit_phi;
          new_hit->partIdx = my_partIdx;
          new_hit->vxTruth = candidate->vxTruth;
          fOutputArray->Add(new_hit);

          cumul_distance = 0;
        }
      }

      if(cumul_distance > stop_criterion){
        std::cout << "stopping at " << pos_r << " " << z_t << std::endl;
        break;
      }
    }

    my_partIdx++;
  }
}

//------------------------------------------------------------------------------
