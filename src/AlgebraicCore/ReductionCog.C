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
