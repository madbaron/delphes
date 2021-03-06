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


/** \class TrackEfficiency
 *
 *  Selects candidates from the InputArray according to the Efficiency formula.
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "modules/TrackEfficiency.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesLongFormula.h"

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

TrackEfficiency::TrackEfficiency() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesLongFormula;
}

//------------------------------------------------------------------------------

TrackEfficiency::~TrackEfficiency()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void TrackEfficiency::Init()
{
  // read TrackEfficiency formula

  fFormula->Compile(GetString("EfficiencyFormula", "1.0"));

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void TrackEfficiency::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TrackEfficiency::Process()
{ 
  Candidate *candidate;
  Double_t pt, eta, phi, e, d0, dz, ctgTheta;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    d0 = candidate->D0;
    dz = candidate->DZ;
    ctgTheta = candidate->CtgTheta;

    // apply an efficency formula
    if(gRandom->Uniform() > fFormula->Eval(pt, eta, phi, e, d0, dz, ctgTheta)) continue;
    
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
