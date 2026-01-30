//   Copyright (c)  2005,2009-2018  John Abbott and Anna M. Bigatti

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


// Source code for functions and member function on monomial ideals

#include "CoCoA/SparsePolyOps-ideal-monomial.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/PPWithMask.H"  // for monomial ideals
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/RingDistrMPolyClean.H" // for NewPolyRing_DMP
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingTwinFloat.H" // for IsRingTwinFloat
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/TmpGOperations.H"  // for myGcd
#include "CoCoA/TmpPPVector.H"  // for interreduce(PPs) etc
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"     // for myGcd
#include "CoCoA/symbol.H"
#include "CoCoA/verbose.H"

#include <algorithm>
using std::swap;
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
#include <list>
// using std::list;
#include <string>
// using std::string;
//#include <vector>
using std::vector;

namespace CoCoA
{

  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasis_MonId() const
  {
    CoCoA_ASSERT(IhaveMonomialGens());
    if (IhaveGBasis()) return myGBasisValue;
    CoCoA_ASSERT(IsField(CoeffRing(myRing())));
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    // convert input into PPVector, operate, convert back
    const SparsePolyRing P = myRing();
    PPVector g(PPM(P), NewDivMaskEvenPowers());
    convert(g, myGens());
    interreduce(g);
    convert(myGBasisValue, P, g);
    IhaveGBasisValue = true;
    // myMinGensValue = myGBasisValue; // not necessary
    return myGBasisValue;
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myMinGens_MonId() const
  { return myGBasis_MonId(); }   // interreduced


  void SparsePolyRingBase::IdealImpl::myTestIsRadical_MonId() const
  {
    VerboseLog VERBOSE("myTestIsRadical_MonId");
    VERBOSE(1000) << " starting " << std::endl;
    if (!IhaveMonomialGens()) CoCoA_THROW_ERROR1(ERR::ReqMonomialGens);
    const std::vector<RingElem>& GB = myGBasis(NoCpuTimeLimit());
    const SparsePolyRing P = myRing();
    // *** assumes GB minimal/interreduced ***
    for (const RingElem& f: GB)
      if (!IsSqFree(LPP(f))) { myAssignRadicalFlag(false); return; }
    myAssignRadicalFlag(true);
  }


  ideal SparsePolyRingBase::IdealImpl::myRadical_MonId() const
  {
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(IsField(CoeffRing(myRing())));
    VerboseLog VERBOSE("myRadical_MonId");
    VERBOSE(1000) << " starting" << std::endl;
    std::vector<RingElem> RadGens;
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector RadGensPPV(PPM(P), NewDivMaskEvenPowers());
      for (const RingElem& f: myGens())  PushBack(RadGensPPV, radical(LPP(f)));
      interreduce(RadGensPPV);
      convert(RadGens, P, RadGensPPV);
    }
    // assign into new ideal RadI
    ideal RadI(myRing(), RadGens); // assignment, copy
    std::swap(ourGetPtr(RadI)->myGBasisValue, RadGens); // assignment, no copy
    ourGetPtr(RadI)->IhaveGBasisValue = true;
    ourGetPtr(RadI)->IhaveMonomialGens3Value = true3;
    ourGetPtr(RadI)->IhaveSqFreeMonomialGens3Value = true3;
    ourGetPtr(RadI)->myAssignRadicalFlag(true);
    return RadI;
  }  


  void SparsePolyRingBase::IdealImpl::myIntersect_MonId(const ideal& J)
  {
    CoCoA_ASSERT(IsField(CoeffRing(myRing())));
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(AreGensMonomial(J));
    VerboseLog VERBOSE("myIntersect_MonId");
    VERBOSE(1000) << " starting" << std::endl;
//    CoCoA_ASSERT(IhaveMonomialGens());
//    CoCoA_ASSERT(AreGensMonomial(J));
    if (IamZero()) return;
    std::vector<RingElem> res;
    if (!gens(J).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g1(PPM(P), NewDivMaskEvenPowers(), myGens());
      PPVector g2(PPM(P), DMR(g1), gens(J));
      PPVector g(PPM(P), DMR(g1));
      g.myLcms(g1, g2);
      interreduce(g);
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    if ((!IsTrue3(IhaveSqFreeMonomialGens3Value))
        || (!IsTrue3(ourGetPtr(J)->IhaveSqFreeMonomialGens3Value)))
      IhaveSqFreeMonomialGens3Value = uncertain3; // can do better than this..
    std::swap(myGensValue, res); // assignment
    myGBasisValue = myGensValue;
    IhaveGBasisValue = true;
    // myMinGensValue = myGBasisValue; // not necessary
  }


  void SparsePolyRingBase::IdealImpl::myColon_MonId(const ideal& J)
  {
    VerboseLog VERBOSE("myColon_MonId");
    VERBOSE(1000) << " starting" << std::endl;
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(AreGensMonomial(J));
    if (IamZero()) return;
    std::vector<RingElem> res;
    if (!gens(J).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g(PPM(P), NewDivMaskEvenPowers());
      PPVector g1(PPM(g), DMR(g), myGens());
      PPVector g2(PPM(g), DMR(g), gens(J));
      PPVector tmp(PPM(g), DMR(g));
      const long len1 = len(g1);
      const long len2 = len(g2);
      for (long i=0; i<len2; ++i)
      {
        tmp.myClear();
        for (long j=0; j<len1; ++j)
          PushBack(tmp, colon(PP(g1[j]),PP(g2[i])));
        interreduce(tmp);
        if (i==0) swap(g,tmp);
        else lcms(g, g, tmp);
        interreduce(g);
      }
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    std::swap(myGensValue, res); // assignment
    if (!IsTrue3(IhaveSqFreeMonomialGens3Value))
      IhaveSqFreeMonomialGens3Value = uncertain3; // can do better than this..
    myGBasisValue = myGensValue;
    IhaveGBasisValue = true;
    // myMinGensValue = myGBasisValue; // not necessary
  }


  void SparsePolyRingBase::IdealImpl::myMul_MonId(const ideal& J)
  {
    CoCoA_ASSERT(IsField(CoeffRing(RingOf(J)))); // ASSUMES CoeffRing is FIELD
    VerboseLog VERBOSE("myMul_MonId");
    VERBOSE(1000) << " starting" << std::endl;
    CoCoA_ASSERT(IhaveMonomialGens());
    CoCoA_ASSERT(AreGensMonomial(J));
    if (IamZero()) return;
    std::vector<RingElem> res;
    if (!gens(J).empty())
    {
      // convert input into PPVector, operate, convert back
      const SparsePolyRing P = myRing();
      PPVector g(PPM(P), NewDivMaskEvenPowers());
      PPVector g1(PPM(g), DMR(g), myGens());
      PPVector g2(PPM(g), DMR(g), gens(J));
      for (long i=len(g2)-1; i>=0; --i)
        for (long j=len(g1)-1; j>=0; --j)
          PushBack(g, PP(g1[j])*PP(g2[i]));
      // need an iterator on PPVector for doing this:
//       for (const auto& f1: g1)
//         for (const auto& f2: g2)
//           PushBack(g, PP(f1)*PP(f2));

      interreduce(g);
      convert(res, P, g);
    }
    // assign into this ideal
    // IhaveMonomialGens3Value = true; // stays constant
    IhaveSqFreeMonomialGens3Value = uncertain3; // can do better than this..
    std::swap(myGensValue, res); // assignment
    myGBasisValue = myGensValue;
    IhaveGBasisValue = true;
    // myMinGensValue = myGBasisValue; // not necessary
  }

  
  void SparsePolyRingBase::IdealImpl::myAssignElim_MonId(const ideal& I, const std::vector<RingElem>& ElimIndets)
  {
    VerboseLog VERBOSE("myAssignElim_MonId");
    VERBOSE(100) << "-- called --" << std::endl;
    CoCoA_ASSERT_ALWAYS(myRing()==RingOf(I));
    CoCoA_ASSERT_ALWAYS(!(IsZero(I)));
    CoCoA_ASSERT(AreGensMonomial(I));
    std::vector<RingElem> ElimGens;
    {
      // convert input into PPVectors, operate, convert back
      const SparsePolyRing P = RingOf(I);
      PPVector g(PPM(P), NewDivMaskEvenPowers());
      PPVector g1(PPM(g), DMR(g), gens(I));
      PPVector g2(PPM(g), DMR(g), ElimIndets);
      const long len1 = len(g1);
      for (long i=0; i<len1; ++i)
        if (!IsDivisible(g1[i], g2)) PushBack(g, g1[i]);
      interreduce(g);
      convert(ElimGens, P, g); // might throw
    }
    // now assign into this ideal
    myReset();
    std::swap(myGensValue, ElimGens); // assignment
    myGBasisValue = myGensValue;
    IhaveGBasisValue = true;
    IhaveMonomialGens3Value = true;
    //    if (IsTrue3(AreGensSqFreeMonomial3(I)))
    //      IhaveSqFreeMonomialGens3Value = true3; // and radical
  }


  ideal IndetsIdeal(const PolyRing& P, ConstRefPPMonoidElem pp)
  {
    vector<RingElem> g;
    for (long i=0 ; i < NumIndets(owner(pp)) ; ++i )
      if ( exponent(pp, i) != 0 )
      {
        if ( exponent(pp, i) != 1 )
          CoCoA_THROW_ERROR2(ERR::BadArg, "power-product must be square free");
        g.push_back(indet(P, i));
      }
    return ideal(P, g);
    // IhaveMonomialGens3Value = true;
    // IhaveSqFreeMonomialGens = true;
  }


  ideal AlexanderDual(const ideal& I)
  {
    VerboseLog VERBOSE("AlexanderDual");
    VERBOSE(1000) << " starting" << std::endl;
    const SparsePolyRing P = RingOf(I);
    if (!AreGensMonomial(I)) CoCoA_THROW_ERROR1(ERR::ReqMonomialGens);
    if (!AreGensSqFreeMonomial(I)) CoCoA_THROW_ERROR2(ERR::NYI, "non-square-free ideal");
    DivMaskRule DMR = NewDivMaskEvenPowers();
    PPVector g(PPM(P), DMR, gens(I));
    PPVector AD(PPM(P), DMR);
    AD.myAlexanderDual(g);
    vector<RingElem> res;
    convert(res, P, AD);
    return ideal(P, res);
  }


  vector<ideal> SparsePolyRingBase::IdealImpl::myPrimaryDecomposition_MonId() const
  {
    VerboseLog VERBOSE("myPrimaryDecomposition_MonId");
    VERBOSE(1000) << " starting" << std::endl;
    if (!IhaveMonomialGens())  CoCoA_THROW_ERROR1(ERR::ReqMonomialGens);
    const SparsePolyRing P = myRing();
    // assumes GB interreduced
    if (!IhaveSqFreeMonomialGens())  CoCoA_THROW_ERROR2(ERR::NYI, "non-square-free ideal");
    const ideal I(const_cast<IdealImpl*>(this));
    const ideal AD = AlexanderDual(I);
    vector<ideal> PD;
    for (const RingElem& f: gens(AD))  // by value
      PD.push_back(IndetsIdeal(P, LPP(f)));
    return PD;
  }



} // end of namespace CoCoA
