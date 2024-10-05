//   Copyright (c)  2006,2021  Anna M. Bigatti

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

#include "CoCoA/ReductionCog.H"

#include "CoCoA/PPMonoid.H"
#include "CoCoA/geobucket.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/assert.H"
#include "CoCoA/ring.H"

#include <cstddef>
//using std::size_t;
#include <iostream>
//using std::ostream;
#include <memory>
//using std::auto_ptr;
//#include <vector>
//using std::vector;


namespace CoCoA
{

  std::ostream& operator<<(std::ostream& out, const ReductionCog& F)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return F->myOutput(out);
  }


  namespace RedCog
  {

    class PolyFieldImpl: public ReductionCogBase
    {
    public:
      PolyFieldImpl(SparsePolyRing P);
      ~PolyFieldImpl()  {};
      void myAssignReset(RingElem& f) override;
      void myAssignReset(RingElem& f, long fLen) override;
      void myRelease(RingElem& f) override;
      ConstRefPPMonoidElem myActiveLPP() const override;
      void myMoveToNextLM() override;
      bool IamActiveZero() const override;
      void myReduce(ConstRefRingElem reducer, long RedLen=0) override;
      std::ostream& myOutput(std::ostream& out) const override;

    private:
      RingElem myActiveSummandsValue;
      RingElem myIgnoredPPsValue;
      RingElem myTmpLM;
    };


    class PolyGCDImpl: public ReductionCogBase
    {
    public:
      PolyGCDImpl(SparsePolyRing P);
      ~PolyGCDImpl()  {};
      void myAssignReset(RingElem& f) override;
      void myAssignReset(RingElem& f, long fLen) override;
      void myRelease(RingElem& f) override;
      ConstRefPPMonoidElem myActiveLPP() const override;
      void myMoveToNextLM() override;
      bool IamActiveZero() const override;
      void myReduce(ConstRefRingElem reducer, long RedLen=0) override;
      std::ostream& myOutput(std::ostream& out) const override;

    private:
      RingElem myActiveSummandsValue;
      RingElem myIgnoredPPsValue;
      RingElem myIgnoredPPsScaleValue;
      RingElem myTmpLM;
      RingElem myTmpScaleValue;
      long myReductionCount;
    };


    class GeobucketFieldImpl: public ReductionCogBase
    {
    public:
      GeobucketFieldImpl(SparsePolyRing P);
      ~GeobucketFieldImpl()  {};
      void myAssignReset(RingElem& f) override;
      void myAssignReset(RingElem& f, long fLen) override;
      void myRelease(RingElem& f) override;
      ConstRefPPMonoidElem myActiveLPP() const override;
      void myMoveToNextLM() override;
      bool IamActiveZero() const override;
      void myReduce(ConstRefRingElem reducer, long RedLen=0) override;
      std::ostream& myOutput(std::ostream& out) const override;

    private:
      geobucket myActiveSummandsValue;
      RingElem myIgnoredPPsValue;
      RingElem myTmpLM;
    };


    class GeobucketGCDImpl: public ReductionCogBase
    {
    public:
      GeobucketGCDImpl(SparsePolyRing P);
      ~GeobucketGCDImpl()  {};
      void myAssignReset(RingElem& f) override;
      void myAssignReset(RingElem& f, long fLen) override;
      void myRelease(RingElem& f) override;
      ConstRefPPMonoidElem myActiveLPP() const override;
      void myMoveToNextLM() override;
      bool IamActiveZero() const override;
      void myReduce(ConstRefRingElem reducer, long RedLen=0) override;
      std::ostream& myOutput(std::ostream& out) const override;

    private:
      geobucket myActiveSummandsValue;
      RingElem myIgnoredPPsValue;
      RingElem myIgnoredPPsScaleValue;
      RingElem myTmpLM;
      RingElem myTmpScaleValue;
      long myReductionCount;
    };

  }  // namespace RedCog




//--  PolyFieldImpl  ----------------------------------

//   RedCog::PolyFieldImpl::PolyFieldImpl(ConstRefRingElem f):
//     myActiveSummandsValue(f), myIgnoredPPsValue(owner(f))
//   {}


  RedCog::PolyFieldImpl::PolyFieldImpl(SparsePolyRing P):
    ReductionCogBase(),
    myActiveSummandsValue(P), myIgnoredPPsValue(P), myTmpLM(P)
  {}


  void RedCog::PolyFieldImpl::myAssignReset(RingElem& f)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    CoCoA_ASSERT( P->myIsValid(raw(f)) );
    P->myAssignZero(raw(myActiveSummandsValue));
    P->myAssignZero(raw(myIgnoredPPsValue));
    swap(f, myActiveSummandsValue);
  }


  void RedCog::PolyFieldImpl::myAssignReset(RingElem& f, long /*fLen*/)
  {
    myAssignReset(f);
  }


  void RedCog::PolyFieldImpl::myRelease(RingElem& f)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    P->myAddClear(raw(myIgnoredPPsValue), raw(myActiveSummandsValue));
    P->myAssignZero(raw(f));
    //    if ( !IsZero(myIgnoredPPsValue) )
    //      P->myDivByCoeff(raw(myIgnoredPPsValue), raw(LC(myIgnoredPPsValue)));
    swap(f, myIgnoredPPsValue);
  }


  ConstRefPPMonoidElem RedCog::PolyFieldImpl::myActiveLPP() const
  { return LPP(myActiveSummandsValue); }


  void RedCog::PolyFieldImpl::myMoveToNextLM()
  {
    SparsePolyRing(owner(myTmpLM))->myMoveLMToBack(raw(myIgnoredPPsValue), raw(myActiveSummandsValue));
  }


  bool RedCog::PolyFieldImpl::IamActiveZero() const
  {  return IsZero(myActiveSummandsValue); }


  void RedCog::PolyFieldImpl::myReduce(ConstRefRingElem g, long /*gLen*/)
  {
    CoCoA_ASSERT( !IamActiveZero() );
    const SparsePolyRing P = owner(g);
    CoCoA_ASSERT( IsField(CoeffRing(P)) );
    CoCoA_ASSERT( P == owner(myActiveSummandsValue) );
    P->myReductionStep(raw(myActiveSummandsValue), raw(g));
  }

  std::ostream& RedCog::PolyFieldImpl::myOutput(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "(" << myIgnoredPPsValue
               << ") + (" << myActiveSummandsValue << ")";
  }


//--  PolyGCDImpl  ----------------------------------

//   RedCog::PolyGCDImpl::PolyGCDImpl(ConstRefRingElem f):
//     myActiveSummandsValue(f), myIgnoredPPsValue(owner(f))
//   {}


  RedCog::PolyGCDImpl::PolyGCDImpl(SparsePolyRing P):
    ReductionCogBase(),
    myActiveSummandsValue(P),
    myIgnoredPPsValue(P),
    myIgnoredPPsScaleValue(one(CoeffRing(P))),
    myTmpLM(P),
    myTmpScaleValue(one(CoeffRing(P)))
  { myReductionCount = 0; }


  void RedCog::PolyGCDImpl::myAssignReset(RingElem& f)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    CoCoA_ASSERT( P->myIsValid(raw(f)) );
    P->myAssignZero(raw(myActiveSummandsValue));
    P->myAssignZero(raw(myIgnoredPPsValue));
    myIgnoredPPsScaleValue = 1;
    //    myTmpScaleValue = 1; // does not matter
    myReductionCount = 0;
    swap(f, myActiveSummandsValue);
  }


  void RedCog::PolyGCDImpl::myAssignReset(RingElem& f, long /*fLen*/)
  {
    myAssignReset(f);
  }


  void RedCog::PolyGCDImpl::myRelease(RingElem& f)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    P->myMulByCoeff(raw(myIgnoredPPsValue), raw(myIgnoredPPsScaleValue));
    //myIgnoredPPsScaleValue = 1;
    P->myAddClear(raw(myIgnoredPPsValue), raw(myActiveSummandsValue));
    P->myAssignZero(raw(f));
    if ( !IsZero(myIgnoredPPsValue) )
      P->myRemoveBigContent(raw(myIgnoredPPsValue));
    swap(f, myIgnoredPPsValue);
    myReductionCount = 0;
  }


  ConstRefPPMonoidElem RedCog::PolyGCDImpl::myActiveLPP() const
  { return LPP(myActiveSummandsValue); }


  void RedCog::PolyGCDImpl::myMoveToNextLM()
  {
    const SparsePolyRing P = owner(myTmpLM);

    P->myMoveLMToFront(raw(myTmpLM), raw(myActiveSummandsValue));
    P->myMulByCoeff(raw(myIgnoredPPsValue), raw(myIgnoredPPsScaleValue));
    myIgnoredPPsScaleValue = 1;
    P->myAppendClear(raw(myIgnoredPPsValue), raw(myTmpLM)); // myTmpLM is 0
  }


  bool RedCog::PolyGCDImpl::IamActiveZero() const
  {  return IsZero(myActiveSummandsValue); }


  void RedCog::PolyGCDImpl::myReduce(ConstRefRingElem g, long /*gLen*/)
  {
    CoCoA_ASSERT( !IamActiveZero() );
    const SparsePolyRing P = owner(g);
    CoCoA_ASSERT( IsTrueGCDDomain(CoeffRing(P)) );
    CoCoA_ASSERT( P == owner(myActiveSummandsValue) );

    ++myReductionCount;
    P->myReductionStepGCD(raw(myActiveSummandsValue), raw(g), myTmpScaleValue);
    if ( IamActiveZero() ) return;
    if ( !IsZero(myIgnoredPPsValue) )
    {
      if ( !IsOne(myTmpScaleValue) )
        myIgnoredPPsScaleValue *= myTmpScaleValue;
    }
    else
      if ( myReductionCount==50 )
      {
        P->myRemoveBigContent(raw(myActiveSummandsValue));
        myReductionCount = 0;
      }
  }


  std::ostream& RedCog::PolyGCDImpl::myOutput(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "(" << myIgnoredPPsScaleValue
               << ")*(" << myIgnoredPPsValue  << ")"
               << ") + (" << myActiveSummandsValue << ")";
  }


//--  GeobucketFieldImpl  ----------------------------------

//   RedCog::GeobucketFieldImpl::GeobucketFieldImpl(ConstRefRingElem f):
//     myActiveSummandsValue(f), myIgnoredPPsValue(owner(f))
//   {}


  RedCog::GeobucketFieldImpl::GeobucketFieldImpl(SparsePolyRing P):
    ReductionCogBase(),
    myActiveSummandsValue(P), myIgnoredPPsValue(P), myTmpLM(P)
  {}


  void RedCog::GeobucketFieldImpl::myAssignReset(RingElem& f)
  {
    myAssignReset(f, NumTerms(f));
  }


  void RedCog::GeobucketFieldImpl::myAssignReset(RingElem& f, long fLen)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    CoCoA_ASSERT( P->myIsValid(raw(f)) );
    //    P->myAssignZero(raw(myActiveSummandsValue));  // to be added to geobucket!!!!
    P->myAssignZero(raw(myIgnoredPPsValue));
    myActiveSummandsValue.myAddClear(f, fLen);
  }


  void RedCog::GeobucketFieldImpl::myRelease(RingElem& f)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    AddClear(myIgnoredPPsValue, myActiveSummandsValue);
    P->myAssignZero(raw(f));
    //    if ( !IsZero(myIgnoredPPsValue) )
    //      P->myDivByCoeff(raw(myIgnoredPPsValue), raw(LC(myIgnoredPPsValue)));
    swap(f, myIgnoredPPsValue);
  }


  ConstRefPPMonoidElem RedCog::GeobucketFieldImpl::myActiveLPP() const
  { return LPP(myActiveSummandsValue); }


  void RedCog::GeobucketFieldImpl::myMoveToNextLM()
  {
    const SparsePolyRing P = owner(myTmpLM);

    MoveLMToFront(myTmpLM, myActiveSummandsValue);
    P->myAppendClear(raw(myIgnoredPPsValue), raw(myTmpLM)); // myTmpLM is 0
  }


  bool RedCog::GeobucketFieldImpl::IamActiveZero() const
  {  return IsZero(myActiveSummandsValue); }


  void RedCog::GeobucketFieldImpl::myReduce(ConstRefRingElem g, long gLen)
  {
    CoCoA_ASSERT( !IamActiveZero() );
    const SparsePolyRing P = owner(g);
    CoCoA_ASSERT( IsField(CoeffRing(P)) );
    CoCoA_ASSERT( P == owner(myIgnoredPPsValue) );
    CoCoA::ReductionStep(myActiveSummandsValue, g, gLen);
  }

  std::ostream& RedCog::GeobucketFieldImpl::myOutput(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "(" << myIgnoredPPsValue
               << ") + (" << myActiveSummandsValue << ")";
  }




//--  GeobucketGCDImpl  ----------------------------------

//   RedCog::GeobucketGCDImpl::GeobucketGCDImpl(ConstRefRingElem f):
//     myActiveSummandsValue(f), myIgnoredPPsValue(owner(f))
//   {}


  RedCog::GeobucketGCDImpl::GeobucketGCDImpl(SparsePolyRing P):
    ReductionCogBase(),
    myActiveSummandsValue(P),
    myIgnoredPPsValue(P),
    myIgnoredPPsScaleValue(one(CoeffRing(P))),
    myTmpLM(P),
    myTmpScaleValue(one(CoeffRing(P)))
  { myReductionCount = 0; }


  void RedCog::GeobucketGCDImpl::myAssignReset(RingElem& f)
  {
    myAssignReset(f, NumTerms(f));
  }


  void RedCog::GeobucketGCDImpl::myAssignReset(RingElem& f, long fLen)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    CoCoA_ASSERT( P->myIsValid(raw(f)) );
    //    P->myAssignZero(raw(myActiveSummandsValue));  // to be added to geobucket!!!!
    P->myAssignZero(raw(myIgnoredPPsValue));
    myActiveSummandsValue.myAddClear(f, fLen);
    myIgnoredPPsScaleValue = 1;
    //    myTmpScaleValue = 1; // does not matter
    myReductionCount = 0;
  }


  void RedCog::GeobucketGCDImpl::myRelease(RingElem& f)
  {
    const SparsePolyRing P = owner(myIgnoredPPsValue);
    CoCoA_ASSERT( P == owner(f) );
    P->myMulByCoeff(raw(myIgnoredPPsValue), raw(myIgnoredPPsScaleValue));
    //myIgnoredPPsScaleValue = 1;
    AddClear(myIgnoredPPsValue, myActiveSummandsValue);
    P->myAssignZero(raw(f));
    if ( !IsZero(myIgnoredPPsValue) )
      P->myRemoveBigContent(raw(myIgnoredPPsValue));
    swap(f, myIgnoredPPsValue);
    myReductionCount = 0;
  }


  ConstRefPPMonoidElem RedCog::GeobucketGCDImpl::myActiveLPP() const
  { return LPP(myActiveSummandsValue); }


  void RedCog::GeobucketGCDImpl::myMoveToNextLM()
  {
    const SparsePolyRing P = owner(myTmpLM);
    RingElem cnt(content(myActiveSummandsValue));

    if ( !IsZero(myIgnoredPPsValue) )
    {
      // cnt = gcd(cnt, myIgnoredPPsScaleValue*content(myIgnoredPPsValue));
      cnt = gcd(cnt, myIgnoredPPsScaleValue);
      myIgnoredPPsScaleValue /= cnt;
      P->myMulByCoeff(raw(myIgnoredPPsValue), raw(myIgnoredPPsScaleValue));
      myIgnoredPPsScaleValue = 1;
    }
    myActiveSummandsValue.myDivByCoeff(cnt);
    MoveLMToFront(myTmpLM, myActiveSummandsValue);
    P->myAppendClear(raw(myIgnoredPPsValue), raw(myTmpLM)); // myTmpLM is 0
  }


  bool RedCog::GeobucketGCDImpl::IamActiveZero() const
  {  return IsZero(myActiveSummandsValue); }


  void RedCog::GeobucketGCDImpl::myReduce(ConstRefRingElem g, long gLen)
  {
    CoCoA_ASSERT( !IamActiveZero() );
    const SparsePolyRing P = owner(g);
    CoCoA_ASSERT( IsTrueGCDDomain(CoeffRing(P)) );
    CoCoA_ASSERT( P == owner(myIgnoredPPsValue) );

    ++myReductionCount;
    CoCoA::ReductionStepGCD(myActiveSummandsValue, g, myTmpScaleValue, gLen);
    if ( IamActiveZero() ) return;
    if ( (!IsZero(myIgnoredPPsValue)) )
    {
      if ( !IsOne(myTmpScaleValue) )
        myIgnoredPPsScaleValue *= myTmpScaleValue;
    }
    else
      if ( myReductionCount==50 )
      {
        RemoveBigContent(myActiveSummandsValue);
        myReductionCount = 0;
      }
  }


  std::ostream& RedCog::GeobucketGCDImpl::myOutput(std::ostream& out) const
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "(" << myIgnoredPPsValue
               << ") + (" << myActiveSummandsValue << ")";
  }



  //----------------------------------------------------------------------
  // Here are the pseudo-constructors:

  ReductionCog NewRedCogPolyField(const SparsePolyRing& P)
  { return ReductionCog(new RedCog::PolyFieldImpl(P)); }

  ReductionCog NewRedCogPolyGCD(const SparsePolyRing& P)
  { return ReductionCog(new RedCog::PolyGCDImpl(P)); }

  ReductionCog NewRedCogGeobucketField(const SparsePolyRing& P)
  { return ReductionCog(new RedCog::GeobucketFieldImpl(P)); }

  ReductionCog NewRedCogGeobucketGCD(const SparsePolyRing& P)
  { return ReductionCog(new RedCog::GeobucketGCDImpl(P)); }


} // namespace CoCoA



//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ReductionCog.C,v 1.18 2022/02/18 14:11:57 abbott Exp $
// $Log: ReductionCog.C,v $
// Revision 1.18  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2021/10/30 16:48:28  abbott
// Summary: Used keywords override and delete (redmine 1625 & 1627)
//
// Revision 1.16  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.15  2018/04/10 14:30:49  bigatti
// -- fixed includes
//
// Revision 1.14  2018/04/09 16:34:26  bigatti
// --
// back to original name
//
// Revision 1.2  2018/04/06 15:35:53  bigatti
// -- fixed includes
//
// Revision 1.1  2018/04/06 15:14:50  bigatti
// renamed ReductionCog.C
//
// Revision 1.12  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.11  2015/04/27 14:41:45  bigatti
// -- PolyFieldImpl::myMoveToNextLM is now one line
//
// Revision 1.10  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.9  2014/07/07 13:45:24  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.8  2014/04/30 16:10:59  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.7  2012/10/16 09:48:31  abbott
// Replaced  RefRingElem  by  RingElem&  (several times)
//
// Revision 1.6  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.5  2009/10/27 13:40:34  bigatti
// -- added check IamActiveZero
//
// Revision 1.4  2009/10/26 16:21:24  bigatti
// -- fixed rare bug in GeobucketGCDImpl::myReduce
//
// Revision 1.3  2009/10/02 13:48:12  bigatti
// -- aesthetics
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/07 17:23:46  cocoa
// -- for compilation with _Wextra: commented out names of unused arguments
//
// Revision 1.5  2006/11/24 17:53:14  cocoa
// -- using SmartPtrIRC for ref counting (was auto_ptr)
//
// Revision 1.4  2006/10/06 15:17:08  cocoa
// -- myReduce now returns void instead of RingElem
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/07/19 07:15:40  cocoa
// -- added CoCoA_ASSERT for input validity in myAssignReset
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/04/11 16:05:56  cocoa
// -- changed: the polynomial is made monic in GPoly::myReduce
//    instead of ReductionCog::myRelease
//
// Revision 1.5  2006/04/11 14:16:29  cocoa
// -- reorganization of fns between reduce,SparsePolyRing,GPoly
// -- added optional "len" argument to myAssignReset in ReductionCog
//
// Revision 1.4  2006/04/10 13:38:06  cocoa
// -- added myTmpLM field to avoid repeated creations/deletions
//
// Revision 1.3  2006/03/17 18:17:16  cocoa
// -- changed: use of ReductionCog for reduction (major cleanup)
//
// Revision 1.2  2006/03/16 15:21:53  cocoa
// -- rearranged
//
// Revision 1.1  2006/03/16 14:28:29  cocoa
// -- first import
//
