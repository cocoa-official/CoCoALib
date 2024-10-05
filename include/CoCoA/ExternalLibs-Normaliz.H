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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ExternalLibs-Normaliz.H,v 1.31 2022/02/25 10:36:37 abbott Exp $
// $Log: ExternalLibs-Normaliz.H,v $
// Revision 1.31  2022/02/25 10:36:37  abbott
// Summary: Removed unnecessary includes from header filed
//
// Revision 1.30  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.29  2020/10/02 18:41:46  abbott
// Summary: Renamed Generators to TriangulationGenerators; cleaned include directive
//
// Revision 1.28  2018/04/10 15:46:08  bigatti
// -- fixed includes
//
// Revision 1.27  2018/04/06 15:05:24  bigatti
// -- fixed include
//
// Revision 1.26  2017/04/28 18:20:58  bigatti
// ++ verbosity for Normaliz: now as for cocoa, NmzSetVerbosityLevel
//
// Revision 1.25  2017/03/13 12:17:08  abbott
// Summary: Include PREPROCESSOR_DEFNS.H (in header files which need the info)
//
// Revision 1.24  2015/09/03 10:17:11  bigatti
// -- changes by Christof Soeger (in Aarhus)
//
// Revision 1.23  2014/07/14 10:02:56  abbott
// Summary: Christof has added some fns which use quasi polys
// Author: JAA
//
// Revision 1.22  2014/07/11 10:12:25  abbott
// Summary: Christof added printing for cones
// Author: JAA
//
// Revision 1.21  2014/07/07 11:44:02  abbott
// Summary: Checking in Chritof's additions
// Author: JAA
//
// Revision 1.20  2014/05/12 14:33:51  bigatti
// -- updated from Christof Soeger
//
// Revision 1.19  2014/05/09 14:56:21  bigatti
// -- new fn by Christof Soeger
//
// Revision 1.18  2014/01/30 09:56:30  bigatti
// -- removed DOS newlines
//
// Revision 1.17  2014/01/29 17:38:00  bigatti
// -- added HilbertSeries (by Christof Soeger)
// -- added multiplicity (by Christof Soeger)
//
// Revision 1.16  2012/10/08 13:52:02  bigatti
// -- more cleaning and updates by Christof Soeger
//
// Revision 1.15  2012/10/05 10:17:04  bigatti
// by Christof Soeger
// * Made the NewCone functions to constructors of cone.
// * Introduced some abbreviation in the method names:
//  IntegralClosure -> IntClosure
//  MonomialIdeal -> MonIdeal
//  Normaliz -> Nmz  (the prefix for CoCoA5 functions)
// * New function IntClosureMonIdeal
//
// Revision 1.14  2012/08/03 16:31:22  bigatti
// -- changed: procedural --> functional (by C.Soeger)
//
// Revision 1.13  2012/07/25 12:47:39  bigatti
// -- merged with C.Soeger changes for NormalizComputation (init from map, etc)
//
// Revision 1.12  2012/07/19 17:12:05  abbott
// Added NewCone -- unified pseudo-ctor so user does not have to choose between long and BigInt.
//
// Revision 1.11  2012/06/19 14:42:30  bigatti
// -- unix-ified
//
// Revision 1.10  2012/06/19 14:41:35  bigatti
// -- Ht1 --> Deg1 (by C.Soeger)
//
// Revision 1.9  2011/11/07 11:09:32  bigatti
// -- new ctors taking ConstRefRingElem
// -- added PPVector NormalToricRing(const PPVector& ppv)
// -- removed void HilbertBasis(std::vector<std::vector<BigInt> >& v, const cone& c)
//
// Revision 1.8  2011/10/04 13:03:02  bigatti
// -- new logo for gui
//
// Revision 1.7  2011/09/30 12:55:44  bigatti
// -- introduced namespace "Normaliz" and removed Normaliz from function names
// -- input of Normaliz functions in CoCoA-5 is now a matrix instead of
//    vector<vector<BigInt>>
//
// Revision 1.6  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.5  2011/07/20 15:33:05  bigatti
// -- changed NORMALIZ_DIR to have no "libnormaliz/"
//
// Revision 1.4  2011/07/20 13:49:37  bigatti
// -- added "Normaliz" postfix to Normaliz function calls
//
// Revision 1.3  2011/07/20 12:45:12  bigatti
// -- new normaliz interface (not yet public)
//
// Revision 1.2  2011/02/17 16:50:04  bigatti
// -- getting ready for new official veson on Normaliz: added HVector, removed Triangulation
//
// Revision 1.1  2010/11/05 14:27:18  bigatti
// -- was TmpNormaliz**
//
// Revision 1.5  2010/10/12 11:22:54  bigatti
// -- TmpNormaliz.H simplified:
//    now cone is a smart pointer to ConeBase
//    and the concrete classes are entirely in the .C file
// -- added Ht1Elements, SupportHyperplanes, Triangulation
// -- added some text in ex-Normaliz1 and 2
//
// Revision 1.4  2010/10/08 16:34:08  bigatti
// -- activated ConeMPZ
//
// Revision 1.3  2010/10/08 10:39:32  bigatti
// -- extended interface for normaliz
//

#endif // CoCoA_WITH_NORMALIZ
#endif
