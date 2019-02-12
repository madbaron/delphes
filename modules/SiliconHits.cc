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

    double track_R = 10000. * candidateMomentum.Pt() / (0.3 * 20.); // (R in mm, B in kGauss and p in GeV/c
    double track_theta = 2.*TMath::ATan(TMath::Exp(-1.*candidateMomentum.Eta()));
    double r_prime = sqrt(candidatePosition.X() * 1.0E-3 * candidatePosition.X() * 1.0E-3 + candidatePosition.Y() * 1.0E-3 * candidatePosition.Y() * 1.0E-3);

    /*
        std::cout << "***********************" << std::endl;
        std::cout << "     NEW candidate     " << std::endl;
        std::cout << "track_R " << track_R << std::endl;
        std::cout << "track_theta " << track_theta << std::endl;
        std::cout << "r_prime " << r_prime << std::endl;
    */
    bool spiraling = false;

    for (UInt_t nlayer = 0; nlayer < barrel_r.size(); nlayer++) {

      double hit_r = barrel_r[nlayer];
      double hit_z = (hit_r - r_prime) / TMath::Tan(track_theta) + candidatePosition.Z() * 1.0E-3;

      if (2 * track_R > barrel_r[nlayer]) {

        if (fabs(hit_z) > barrel_z[nlayer]) continue;

        double hit_phi = (hit_z - candidatePosition.Z() * 1.0E-3) / (2. * track_R * TMath::Tan(track_theta)) + candidateMomentum.Phi();
        if (candidate->Charge < 0) {
          hit_phi = -(hit_z - candidatePosition.Z() * 1.0E-3) / (2. * track_R * TMath::Tan(track_theta)) + candidateMomentum.Phi();
        }
        hit_phi = fmod(hit_phi + TMath::Pi(), 2.*TMath::Pi());
        if (hit_phi < 0) {
          hit_phi += 2.*TMath::Pi();
        }
        hit_phi -= TMath::Pi();

        /*
          double hit_phi = TMath::ASin(hit_r / (2 * track_R) ) - candidateMomentum.Phi();
          if (candidate->Charge < 0) {
            hit_phi = -TMath::ASin(hit_r / (2 * track_R) ) - candidateMomentum.Phi();
          }
          hit_phi = fmod(hit_phi + TMath::Pi(), 2.*TMath::Pi());
          if (hit_phi < 0) {
            hit_phi += 2.*TMath::Pi();
          }
          hit_phi -= TMath::Pi();
        */
        new_hit = factory->NewCandidate();
        new_hit->Xd = hit_r;
        new_hit->Yd = hit_z;
        new_hit->Zd = hit_phi;
        new_hit->partIdx = my_partIdx;
        new_hit->vxTruth = candidate->vxTruth;

        fOutputArray->Add(new_hit);

      }
      /*
      else {
        if (!spiraling) {

          spiraling = true;
          UInt_t layer_index = nlayer-1;
          int change = -1;

          if(layer_index<0) continue;

          while (fabs(hit_z) < barrel_z[layer_index]) {

              hit_z += (hit_r - r_prime) / TMath::Tan(track_theta);
              hit_r = barrel_r[layer_index];

              double hit_phi = (hit_z - candidatePosition.Z() * 1.0E-3) / (2. * track_R * TMath::Tan(track_theta)) + candidateMomentum.Phi();
              if (candidate->Charge < 0) {
                hit_phi = -(hit_z - candidatePosition.Z() * 1.0E-3) / (2. * track_R * TMath::Tan(track_theta)) + candidateMomentum.Phi();
              }
              hit_phi = fmod(hit_phi + TMath::Pi(), 2.*TMath::Pi());
              if (hit_phi < 0) {
                hit_phi += 2.*TMath::Pi();
              }
              hit_phi -= TMath::Pi();

              new_hit = factory->NewCandidate();
              new_hit->Xd = hit_r;
              new_hit->Yd = hit_z;
              new_hit->Zd = hit_phi;
              new_hit->partIdx = my_partIdx;
              new_hit->vxTruth = candidate->vxTruth;

              fOutputArray->Add(new_hit);

              layer_index += change;
              if(layer_index==0) change = 1;
              if(layer_index==nlayer-1) change = -1;

              if(layer_index<0) break;

            }
          }
        }*/
      }

      for (UInt_t ndisk = 0; ndisk < disk_z.size(); ndisk++) {

        double hit_z = disk_z[ndisk];
        if ( hit_z < candidatePosition.Z() * 1.0E-3 && ((track_theta < TMath::Pi() / 2.) && (track_theta > -TMath::Pi() / 2.)) ) continue;
        if ( hit_z > candidatePosition.Z() * 1.0E-3 && ((track_theta > TMath::Pi() / 2.) || (track_theta < -TMath::Pi() / 2.)) ) continue;

        double z_lenght = hit_z - candidatePosition.Z() * 1.0E-3;
        double hit_r = z_lenght * TMath::Tan(track_theta) + r_prime;
        if (fabs(hit_r) < disk_rmin[ndisk] || fabs(hit_r) > disk_rmax[ndisk]) continue;
        double hit_phi = (hit_z - candidatePosition.Z() * 1.0E-3) / (2. * track_R * TMath::Tan(track_theta)) + candidateMomentum.Phi();
        if (candidate->Charge < 0) {
          hit_phi = -(hit_z - candidatePosition.Z() * 1.0E-3) / (2. * track_R * TMath::Tan(track_theta)) + candidateMomentum.Phi();
        }
        hit_phi = fmod(hit_phi + TMath::Pi(), 2.*TMath::Pi());
        if (hit_phi < 0) {
          hit_phi += 2.*TMath::Pi();
        }
        hit_phi -= TMath::Pi();

        /*
              double hit_phi = TMath::ASin(hit_r / 2.*track_R) - candidateMomentum.Phi();

              if (candidate->Charge < 0) {
                hit_phi = -TMath::ASin(hit_r / 2.*track_R) - candidateMomentum.Phi();
              }
              hit_phi = fmod(hit_phi + TMath::Pi(), 2.*TMath::Pi());
              if (hit_phi < 0) {
                hit_phi += 2.*TMath::Pi();
              }
              hit_phi -= TMath::Pi();
        */

        new_hit = factory->NewCandidate();
        new_hit->Xd = fabs(hit_r);
        new_hit->Yd = hit_z;
        new_hit->Zd = hit_phi;
        new_hit->partIdx = my_partIdx;
        new_hit->vxTruth = candidate->vxTruth;

        fOutputArray->Add(new_hit);
      }

      my_partIdx++;
    }
  }

//------------------------------------------------------------------------------
