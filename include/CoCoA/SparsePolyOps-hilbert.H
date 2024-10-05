#ifndef CoCoA_hilbert_H
#define CoCoA_hilbert_H

//   Copyright (c)  2006  Anna Bigatti

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

#include "CoCoA/factorization.H"
#include "CoCoA/ring.H"
#include "CoCoA/BigInt.H"

namespace CoCoA
{

  //  class ConstRefRingElem; // forward declaration -- defined in ring.H
  //  class RingElem;     // forward declaration -- defined in ring.H
  class QuotientRing; // forward declaration -- defined in SparsePolyRing.H
  class ideal;        // forward declaration -- defined in ideal.H

  class HPSeries
  {
  public:
    HPSeries(ConstRefRingElem num, const factorization<RingElem>& den);
    HPSeries(const std::vector<BigInt>& DenseRepr, const std::vector<long>& DenExponents, long shift = 0);
//  HPSeries(const vector<BigInt>& DenseRepr, long DenPower);

  private:    // Data members
    RingElem myNum;
    factorization<RingElem> myDenFactors;

    friend const RingElem& num(const HPSeries& S);
    friend const factorization<RingElem>& DenFactors(const HPSeries& S);
    friend std::ostream& operator<<(std::ostream& out, const HPSeries& S);
  };


  RingElem HilbertNumQuot_C(const ideal& I);
  RingElem HilbertNumQuot(const ideal& I);
  //  RingElem MGHilbertNumQuot(const SparsePolyRing& HSRing, const ideal& I);
  RingElem MGHilbertNumQuot(const ideal& I);
  void EndPoincare_C(); // just for mempool complaints on global unipoly array

  HPSeries HilbertSeries(const ideal& I);
  HPSeries HilbertSeries(const QuotientRing& PModI);
  HPSeries HilbertSeriesQuot(const ideal& I);
  //  long dim(const HPSeries& HPS); // TODO
  //  long multiplicity(const HPSeries& HPS); // TODO
  RingElem HilbertPoly(const QuotientRing& PModI);
  RingElem HilbertPoly(const HPSeries& HPS);
  std::vector<BigInt> HVector(const QuotientRing& PModI);
  std::vector<BigInt> HVector(const HPSeries& HPS);
  long DimQuot(const ideal& I);
  BigInt MultiplicityQuot(const ideal& I);

  //----------------------------------------------------------------------
  // Inline functions
  inline const RingElem& num(const HPSeries& S)
  { return S.myNum; }

  inline const factorization<RingElem>& DenFactors(const HPSeries& S)
  { return S.myDenFactors; }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-hilbert.H,v 1.2 2022/02/18 14:11:49 abbott Exp $
// $Log: SparsePolyOps-hilbert.H,v $
// Revision 1.2  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2018/04/06 15:01:45  bigatti
// -- renamed TmpHilbert.H
//
// Revision 1.16  2016/10/27 13:55:13  bigatti
// -- added MultiplicityQuot, DimQuot
//
// Revision 1.15  2016/02/16 16:26:55  bigatti
// -- just some reminder for dim and multiplicity
//
// Revision 1.14  2015/09/03 10:17:11  bigatti
// -- changes by Christof Soeger (in Aarhus)
//
// Revision 1.13  2014/10/24 12:47:22  bigatti
// -- added HilbertPoly, HVector (to be refined)
//
// Revision 1.12  2014/01/30 10:01:28  bigatti
// -- removed DOS newlines
//
// Revision 1.11  2014/01/29 17:28:38  bigatti
// -- added HPSeries ctor for DenseRepr (by Christof Soeger)
//
// Revision 1.10  2013/07/30 15:02:37  bigatti
// -- added HPSeries class
// -- added HilbertSeriesQuot
//
// Revision 1.9  2013/07/23 14:26:06  bigatti
// -- fwd decl of SparsePolyRing instead of DenseUPolyRing (it compiled!?!)
//
// Revision 1.8  2013/06/18 12:27:37  abbott
// Renamed HibertSeriesPolyRing to RingQQt.
//
// Revision 1.7  2013/06/17 08:54:39  abbott
// Added HilbertSeriesPolyRing (untested).
//
// Revision 1.6  2011/04/26 10:13:04  bigatti
// -- added MGHilbertNumQuot
//
// Revision 1.5  2011/04/08 14:04:46  bigatti
// -- renamed HilbertNumeratorMod into HilbertNumQuot
//
// Revision 1.4  2010/10/29 11:59:12  bigatti
// -- added function for cleaning memory after using C global
//
// Revision 1.3  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/10 14:36:34  bigatti
// new: TmpHilbertDir/TmpPoincareCPP.[CH]: poincare code using
// 	C++ univariate polynomials (DenseUPolyRing)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.1  2007/01/17 17:38:11  bigatti
// -- moved all cocoa-4 code for hilbert into src/TmpHilbertDir
//
// Revision 1.3  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.2  2006/11/24 17:43:37  cocoa
// -- reorganized includes of header files
//
// Revision 1.1  2006/10/09 16:48:58  cocoa
// -- first import
//

#endif
