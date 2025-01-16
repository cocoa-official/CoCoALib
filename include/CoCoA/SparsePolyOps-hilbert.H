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

#endif
