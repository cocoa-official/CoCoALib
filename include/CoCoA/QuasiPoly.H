#ifndef CoCoA_QuasiPoly_H
#define CoCoA_QuasiPoly_H

//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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


#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;

namespace CoCoA
{

  class RingElem;   // fwd decl -- defined in ring.H
  class BigInt;     // fwd decl -- defined in BigInt.H
  class MachineInt; // fwd decl -- defined in MachineInt.H

  class QuasiPoly
  {
  public:
    explicit QuasiPoly(const std::vector<RingElem>& v);
    RingElem operator()(const MachineInt& n) const;
    RingElem operator()(const BigInt& N) const;
  private: // data members
    std::vector<RingElem> myConstituents;
    friend const std::vector<RingElem>& constituents(const QuasiPoly& p);
  };

  const std::vector<RingElem>& constituents(const QuasiPoly& p);
  std::ostream& operator<<(std::ostream& out, const QuasiPoly& p);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/QuasiPoly.H,v 1.6 2022/02/18 14:11:47 abbott Exp $
// $Log: QuasiPoly.H,v $
// Revision 1.6  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/07 14:21:51  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2015/04/16 16:20:14  abbott
// Summary: Removed useless include directive
// Author: JAA
//
// Revision 1.3  2014/07/14 14:19:51  abbott
// Summary: Removed a rogue "e" at start of the file
// Author: JAA
//
// Revision 1.2  2014/07/14 13:16:25  abbott
// Summary: Minor cleaning
// Author: JAA
//
// Revision 1.1  2014/07/14 10:04:29  abbott
// Summary: Added new impl of quasipolynomials
// Author: JAA
//

#endif
