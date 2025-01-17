#ifndef CoCoA_F5_H
#define CoCoA_F5_H

//   Copyright (c) 2007 Alberto Arri

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

#include "CoCoA/ring.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/SparsePolyRing.H"

namespace CoCoA
{

  namespace F5ns{
    
  struct F5opt_t{
    bool incremental, skip_rows, prev_red, use_NR, GBLT2SYZLT, checkGB;
    F5opt_t(){
      incremental = true; //default = true
      skip_rows = true;   //default = true
      prev_red = true;    //default = true
      use_NR = false;     //makes sense only if skip_rows = true
      GBLT2SYZLT = true;  //default = true
      checkGB = false;    //default = false
    }
  };
    
  } // namespace F5ns

  std::vector<RingElem> F5_mat(const std::vector<RingElem>& GensI, 
                               const F5ns::F5opt_t& F5opt = F5ns::F5opt_t());
  std::vector<RingElem> F5_poly(const std::vector<RingElem>& GensI);  

  inline std::vector<RingElem> F5_mat(const ideal& I)  { return F5_mat(gens(I)); }
  inline std::vector<RingElem> F5_poly(const ideal& I)  { return F5_poly(gens(I)); }

} // end of namespace CoCoA

#endif
