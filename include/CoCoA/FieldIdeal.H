#ifndef CoCoA_FieldIdeal_H
#define CoCoA_FieldIdeal_H

//   Copyright (c)  2005-2007,2009  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include <vector>
// using std::vector;


namespace CoCoA
{

  class RingElem; // forward declaration -- defined in RingElem.H
  class ideal;    // forward declaration -- defined in ideal.H
  class ring;     // forward declaration -- defined in ring.H

  /*-----------------------------------------------------------------*/
  /**
     \brief FieldIdealImpl -- the simple ideals of a field.

     \include FieldIdeal.txt
  */
  /*-----------------------------------------------------------------*/
  ideal NewFieldIdeal(const ring& k, const std::vector<RingElem>& gens);


} // end of namespace CoCoA

#endif
