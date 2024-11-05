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
#include <functional>
using std::not1;    // for IsMonomial
using std::ptr_fun; // for IsMonomial
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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-ideal-monomial.C,v 1.15 2024/08/02 14:28:16 bigatti Exp $
// $Log: SparsePolyOps-ideal-monomial.C,v $
// Revision 1.15  2024/08/02 14:28:16  bigatti
// Summary: error message: minor improvement in string
//
// Revision 1.14  2024/07/23 16:40:56  bigatti
// Summary: fixed typo
//
// Revision 1.13  2024/07/23 16:38:30  bigatti
// Summary: error codes: ERR::ReqMonomialGens, and using new CoCoA_THROW_ERROR1/2
//
// Revision 1.12  2024/07/23 13:23:08  bigatti
// Summary: added myMinGens_MonId
//
// Revision 1.11  2024/03/15 15:44:50  bigatti
// Summary: changes myElim_MonId into myAssignElim_MonId (redmine #1824)
//
// Revision 1.10  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2021/07/19 11:15:02  abbott
// Summary: Layout
//
// Revision 1.8  2021/04/28 12:44:03  abbott
// Summary: Added assertions that coeffs are in a field (redmine 1381)
//
// Revision 1.7  2021/04/28 09:03:05  abbott
// Summary: Added assertion in myMul_MonId that coeff ring is a field. (redmine 1381)
//
// Revision 1.6  2021/03/04 17:00:49  abbott
// Summary: Minor change to for loops
//
// Revision 1.5  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2020/02/12 09:01:47  bigatti
// -- changed myTestIsMaximal etc to return void (and consequences)
//
// Revision 1.3  2019/10/21 16:44:52  bigatti
// -- using auto declaration for iterators
// -- using range based for loops
// -- using new ctor with vector<RingElem>
//
// Revision 1.2  2019/10/18 15:07:22  bigatti
// -- using auto decl instead of iterator (in for loops)
//
// Revision 1.1  2019/10/15 12:57:55  bigatti
// -- renamed files for ideals
//
// Revision 1.8  2019/10/03 13:33:56  bigatti
// -- implemented radical for monomial ideals (and used where useful)
//
// Revision 1.7  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.6  2018/05/18 16:38:52  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.5  2018/05/18 12:23:50  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.4  2018/04/20 16:19:36  bigatti
// -- fixed myTestIsRadical_MonId
//
// Revision 1.3  2018/04/10 14:20:10  bigatti
// -- started work on primary decomposition
//
// Revision 1.2  2018/04/09 16:26:26  bigatti
// -- minor cleaning
//
// Revision 1.1  2018/04/06 15:22:44  bigatti
// -- renaming TmpMonomialIdeal.C
//
// Revision 1.20  2018/03/29 09:36:39  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.19  2018/03/20 11:44:23  bigatti
// -- changed: ***MonId --> ***_MonId
//
// Revision 1.18  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.17  2016/11/07 13:54:21  bigatti
// ++ added AreGensSqFreeMonomial
// ++ replaced all SquareFree into SqFree
//
// Revision 1.16  2016/11/07 12:22:56  bigatti
// -- Changed myGBasisIsValid into IhaveGBasis
//
// Revision 1.15  2016/10/27 14:24:53  bigatti
// -- added myRadicalMonId
//
// Revision 1.14  2016/05/23 12:48:15  bigatti
// -- aesthetics
//
// Revision 1.13  2014/07/30 14:11:41  abbott
// Summary: Changed name AmbientRing --> RingOf; myAmbientRing --> myRing
// Author: JAA
//
// Revision 1.12  2014/07/14 15:08:59  abbott
// Summary: Removed include of tmp.H (no longer needed?)
// Author: JAA
//
// Revision 1.11  2014/07/08 15:23:53  abbott
// Summary: Updated comment
// Author: JAA
//
// Revision 1.10  2014/07/07 13:14:25  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.9  2014/03/27 14:58:05  bigatti
// -- just comments to remind that GBasis is minimal generators
//
// Revision 1.8  2012/05/30 13:44:45  bigatti
// -- renamed IhaveMonomialGensB3Value --> IhaveMonomialGens3Value
//
// Revision 1.7  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.6  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.5  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.4  2011/11/07 11:06:14  bigatti
// -- myGBasisMonId: using now convert(g, myGens()); instead of extended code
//
// Revision 1.3  2011/07/27 15:50:56  bigatti
// -- improved use of len in for loops
//
// Revision 1.2  2011/07/05 15:02:17  bigatti
// -- added AlexanderDual
// -- added ad-hoc functions for colon, elim on monomial ideals
//
// Revision 1.1  2011/06/27 13:26:51  bigatti
// -- first import (some functions were in SparsePolyRing)
//
