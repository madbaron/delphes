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

#ifndef DisplacedVertexDumper_h
#define DisplacedVertexDumper_h

/** \class DisplacedVertexDumper
 *
 *  Saves displaced vertices from long lived particle decays
 *
 *  \author K. Albertsson - CERN LTU/EISLAB
 *  \author F. Meloni - DESY
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;

class DisplacedVertexDumper: public DelphesModule
{
public:

  DisplacedVertexDumper();
  ~DisplacedVertexDumper();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!
  TIterator *fItOutputArray; //!
  
  const TObjArray *fInputArray; //!
  const TObjArray *fInputParticles; //!

  TObjArray *fOutputArray; //!

  ClassDef(DisplacedVertexDumper, 1)
};

#endif
