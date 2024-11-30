#ifndef CoCoA_normaliz_H
#define CoCoA_normaliz_H

//   Copyright (c)  2010-2015 John Abbott, Anna M. Bigatti, Christof Soeger
//   Authors:  2010-2015 John Abbott, Anna M. Bigatti, Christof Soeger

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

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_NORMALIZ

//#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/cone_property.h"
#include "libnormaliz/input_type.h"

#include "CoCoA/SmartPtrIRC.H"
#include "CoCoA/SparsePolyOps-hilbert.H"

#include <map>
#include <vector>
#include <iosfwd>

namespace CoCoA
{
  class BigInt;          // forward declaration -- defined in BigInt.H 
  class BigRat;          // forward declaration -- defined in BigRat.H 
  class ConstMatrixView; // forward declaration -- defined in matrix.H
  class PPVector;        // forward declaration -- defined in TmpPPVector.H 
  class PPMonoid;        // forward declaration -- defined in PPMonoid.H 
  class symbol;          // forward declaration -- defined in symbol.H 
  class QuasiPoly;       // forward declaration -- defined in QuasiPoly.H 

namespace Normaliz
{
  class ConeImpl; // fwd decl for SmartPtrIRC
  class cone
  {
  public:

    cone(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m);
    cone(libnormaliz::Type::InputType InputType, const std::vector<std::vector<BigInt> >& v);
    cone(libnormaliz::Type::InputType InputType, ConstMatrixView M);

    explicit cone(const ConeImpl* ptr); //: mySmartPtr(ptr) {}
    // Default copy ctor and dtor works fine, but we implement them to ensure they are not made when including this .H file which will lead to errors because the declaration of ConeImpl is unknown
    ~cone();
    cone(const cone&);
    // Assignment disabled because SmartPtrIRC has no assignment.

    void myComputation(const libnormaliz::ConeProperties& CPs) const;
    void myComputation(libnormaliz::ConeProperty::Enum CP) const;
    void myComputation(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const;
    void myComputation() const;   //default: compute everything possible 
    bool isComputed(libnormaliz::ConeProperty::Enum CP) const;
    bool isComputed(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const;
    const ConeImpl* operator->() const; // { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    //    const ConeImpl* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions
  private: // data members
    SmartPtrIRC<const ConeImpl> mySmartPtr;
  };

  std::ostream& operator<< (std::ostream& out, const cone& C);

  void SetVerbosityLevel(long v); 
  void SetVerbosityLevel(cone& C, long v);
  long VerbosityLevel();
  long VerbosityLevel(cone& c);

  std::vector<std::vector<BigInt> > HilbertBasis(const cone& c);
  std::vector<std::vector<BigInt> > ModuleGenerators(const cone& c);
  std::vector<std::vector<BigInt> > TriangulationGenerators(const cone& c);
  std::vector<std::vector<BigInt> > ExtremeRays(const cone& c);
  std::vector<std::vector<BigInt> > VerticesOfPolyhedron(const cone& c);
  std::vector<std::vector<BigInt> > Deg1Elements(const cone& c);
  std::vector<std::vector<BigInt> > GeneratorsOfToricRing(const cone& c);
  std::vector<std::vector<BigInt> > SupportHyperplanes(const cone& c);
  std::vector<std::vector<BigInt> > Equations(const cone& c);
  std::vector<std::vector<BigInt> > Congruences(const cone& c);
  std::vector<std::vector<BigInt> > ExcludedFaces(const cone& c);
  //std::vector<BigInt> HVector(const cone& c);
  HPSeries HilbertSeries(const cone& c);
  RingElem HilbertPoly(const cone& c);
  QuasiPoly HilbertQuasiPoly(const cone& c);
  BigRat multiplicity(const cone& c);
  std::vector<BigRat> grading(const cone& c);
  bool IsPointed(const cone& c);
  bool IsInhomogeneous(const cone& c);
  bool IsIntegrallyClosed(const cone& c);
  bool IsDeg1HilbertBasis(const cone& c);
  long EmbeddingDim(const cone& c);
  long rank(const cone& c);
  long RecessionRank(const cone& c);
  long AffineDim(const cone& c);
  long ModuleRank(const cone& c);
  std::vector<BigInt> dehomogenization(const cone& c);

  // applications to monomials
  PPVector NormalToricRing(const PPVector& ppv);
  PPVector IntClosureToricRing(const PPVector& ppv);
  PPVector IntClosureMonIdeal(const PPVector& ppv);
  PPVector IntClosureMonIdeal(const PPVector& ppv, long var_pos);
  PPVector EhrhartRing(const PPVector& ppv, long var_pos);
  // further applications
  PPVector TorusInvariants(const std::vector< std::vector<BigInt> >& T, const PPMonoid& ppm);
  PPVector FiniteDiagInvariants(const std::vector< std::vector<BigInt> >& Cong, const PPMonoid& ppm);
  PPVector DiagInvariants(const std::vector< std::vector<BigInt> >& T, const std::vector< std::vector<BigInt> >& Cong, const PPMonoid& ppm);
  PPVector IntersectionValRings (const std::vector< std::vector<BigInt> >& V, const PPMonoid& ppm);

  // conversion functions used in the CoCoA5 Normaliz functions
  std::vector<std::vector<BigInt> > MatrixToVecVecBigInt(ConstMatrixView M);
  libnormaliz::InputType ToInputType(const std::string& TypeString);
  PPVector MonomialsToPPV(const std::vector<RingElem>& v);

 
} // end of namespace Normaliz
} // end of namespace CoCoA

#endif // CoCoA_WITH_NORMALIZ
#endif
