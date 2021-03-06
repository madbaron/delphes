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


/** \class DisplacedVertexDumper
 *
 *  Saves displaced vertices from long lived particle decays
 *
 *  \author K. Albertsson - CERN LTU/EISLAB
 *  \author F. Meloni - DESY
 *
 */

#include "modules/DisplacedVertexDumper.h"

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

DisplacedVertexDumper::DisplacedVertexDumper() :
  fItInputArray(0)
{

}

//------------------------------------------------------------------------------

DisplacedVertexDumper::~DisplacedVertexDumper()
{

}

//------------------------------------------------------------------------------

void DisplacedVertexDumper::Init()
{
  // read parameters

  // import input array(s)
  fInputArray = ImportArray(GetString("InputArray", "TrackMergerAll/tracks"));
  fItInputArray = fInputArray->MakeIterator();
  fInputParticles = ImportArray("Delphes/allParticles");

  // create output array(s)
  fOutputArray = ExportArray(GetString("OutputArray", "displacedvertices"));
  fItOutputArray = fOutputArray->MakeIterator();

}

//------------------------------------------------------------------------------

void DisplacedVertexDumper::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fItOutputArray) delete fItOutputArray;
}

//------------------------------------------------------------------------------

void DisplacedVertexDumper::Process()
{
  Candidate *candidate;
  Candidate *mother;
  Candidate *candidate_out;
  TLorentzVector candidatePosition;
  TLorentzVector candidate_outPosition;
  Int_t motherPID = 999999999;
  Int_t lastmotherPID;
  Int_t firstmotherPID = 0;

  // loop over all input candidates
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {

    candidatePosition = candidate->InitialPosition;

    //If the mother exists, get its PID
    UInt_t mpos = candidate->M1;
    if(mpos>=0 && mpos < (UInt_t)fInputParticles->GetSize()){
      mother = static_cast<Candidate*>(fInputParticles->At(mpos));
      lastmotherPID = mother->PID;
    }
    else{
      lastmotherPID = -1;
    }

    //Now navigate back the chain to the find the first outgoing SUSY particle
    while(mpos >= 0 && mpos < (UInt_t)fInputParticles->GetSize())
    {
      mother = static_cast<Candidate*>(fInputParticles->At(mpos));
      mpos = mother->M1;
      motherPID = mother->PID;
      if((fabs(motherPID)>1000000 && fabs(motherPID)<1000040) || (fabs(motherPID)>2000000 && fabs(motherPID)<2000016)){        
        if (mpos >= 0 && mpos < (UInt_t)fInputParticles->GetSize()){
          mother = static_cast<Candidate*>(fInputParticles->At(mpos));
          mpos = mother->M1;
          firstmotherPID = mother->PID;
          while (mpos >= 0 && mpos < (UInt_t)fInputParticles->GetSize() && firstmotherPID==motherPID){
            mother = static_cast<Candidate*>(fInputParticles->At(mpos));
            mpos = mother->M1;
            firstmotherPID = mother->PID;
          }
        }
        break;
      }
    }

    // Check radial displacement
    double r = TMath::Hypot(candidatePosition.X(), candidatePosition.Y()) * 1.0E3; // in [mm]

    // Save SUSY particles and everything above 4 mm (TODO: make this configurable)
    if( r > 4. || ((fabs(motherPID)>1000000 && fabs(motherPID)<1000040) || (fabs(motherPID)>2000000 && fabs(motherPID)<2000016)) ){

      double min_distance = 99999999999.;
      //Loop over output to avoid duplicates
      fItOutputArray->Reset();
      while((candidate_out = static_cast<Candidate *>(fItOutputArray->Next()))){
            
        candidate_outPosition = candidate_out->InitialPosition;

        double distance = sqrt(
          pow((candidatePosition.X()-candidate_outPosition.X()),2.) + 
          pow((candidatePosition.Y()-candidate_outPosition.Y()),2.) + 
          pow((candidatePosition.Z()-candidate_outPosition.Z()),2.)            
        ) * 1.0E3; // in [mm]

        if ( distance<min_distance ) min_distance = distance;
  
      }

      // Use 1 mm to discriminate between vertices
      if( min_distance>1. ){
        // Add to output array
        candidate = static_cast<Candidate *>(candidate->Clone());
        candidate->InitialPosition = candidatePosition;
        if( ((fabs(motherPID)>1000000 && fabs(motherPID)<1000040) || (fabs(motherPID)>2000000 && fabs(motherPID)<2000016)) ){
          candidate->PID = motherPID;
          candidate->PID_mom = firstmotherPID;
        }
        else{
          candidate->PID = lastmotherPID;
          candidate->PID_mom = 0;
        } 
        fOutputArray->Add(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------
