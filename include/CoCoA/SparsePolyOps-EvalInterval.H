#ifndef CoCoA_SparsePolyOps_EvalInterval_H
#define CoCoA_SparsePolyOps_EvalInterval_H

//   Copyright (c)  2018  John Abbott,  Anna M. Bigatti

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
//using std::vector;

namespace CoCoA
{

  class BigRat;  // fwd decl -- defined in BigRat.H
  class BigRatInterval;  // fwd decl -- defined in BigRatInterval.H

  BigRatInterval eval(const std::vector<BigRat>& f, const BigRatInterval& I);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-EvalInterval.H,v 1.2 2022/02/18 14:11:49 abbott Exp $
// $Log: SparsePolyOps-EvalInterval.H,v $
// Revision 1.2  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2018/04/20 13:12:15  abbott
// Summary: Added new fn eval(poly, interval) -- still only a prototype
//
// Revision 1.6  2017/07/21 13:18:56  abbott
// Summary: Now copyright "template" is Abbott+Bigatti
//
// Revision 1.5  2010/03/22 11:49:01  abbott
// Updated copyright year
//
// Revision 1.4  2009/06/04 14:30:08  abbott
// Changed copyright year.
//
// Revision 1.3  2008/12/16 21:05:54  abbott
// Updated licensing notice from GPL2 to GPL3+ -- evidently I forgot to change these
// files when I updated the others.
//
// Revision 1.2  2008/03/12 14:41:28  abbott
// Updated copyright year
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2007/03/03 15:24:39  cocoa
// Changed 2006 to 2007.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//

#endif
