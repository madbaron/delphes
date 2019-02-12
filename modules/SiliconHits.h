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

#ifndef SiliconHits_h
#define SiliconHits_h

/** \class SiliconHits
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author F. Meloni, DESY
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class SiliconHits: public DelphesModule
{
public:

  SiliconHits();
  ~SiliconHits();

  void Init();
  void Process();
  void Finish();

private:

  std::vector<float> barrel_r;
  std::vector<float> barrel_z;

  std::vector<float> disk_rmin;
  std::vector<float> disk_rmax;
  std::vector<float> disk_z;


  TIterator *fItInputArray; //!
  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(SiliconHits, 1)
};

#endif
