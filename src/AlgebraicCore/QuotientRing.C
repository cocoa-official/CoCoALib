//   Copyright (c)  2003-2009,2011,2014,2021  John Abbott and Anna M. Bigatti

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


// Source code for class GeneralQuotientRingImpl

#include "CoCoA/QuotientRing.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingElemInput.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingFpDouble.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"
//#include "CoCoA/SparsePolyRing.H"

#include <iostream>
using std::ostream;
#include <memory>
using std::unique_ptr;
#include <vector>
using std::vector;


namespace CoCoA
{

  const QuotientRingBase* QuotientRingPtr(const ring& R, const ErrorContext& ErrCtx)
  {
    const QuotientRingBase* ptr = QuotientRingPtr(R);
    if (ptr == nullptr)  CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::NotQuotientRing, ErrCtx);
    return ptr;
  }


  class QuotientingHomImpl: public RingHomBase
  {
  public:
    QuotientingHomImpl(const QuotientRing& RmodI);
    virtual void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const;
    virtual bool IamPartial() const { return false; }
  private:
    virtual void myOutputSelfDetails(std::ostream& out) const;
  };


  RingHom QuotientRingBase::myQuotientingHomCtor() const
  { return RingHom(new QuotientingHomImpl(QuotientRing(this))); }


  RingElem QuotientRingBase::mySymbolValue(const symbol& sym) const
  { return myQuotientingHomCtor()(myReprRing->mySymbolValue(sym)); }


  class GeneralQuotientRingImpl: public QuotientRingBase
  {
  private: // data members
    unique_ptr<RingElemAlias> myZeroPtr;
    unique_ptr<RingElemAlias> myOnePtr;

  public:
    GeneralQuotientRingImpl(const ring& R, const ideal& I);
    // Default copy ctor works fine.
    // Assignment disabled
    ~GeneralQuotientRingImpl() {};
  public: // disable assignment
    GeneralQuotientRingImpl& operator=(const GeneralQuotientRingImpl&) = delete;

  public: // functions that every ring must implement
    BigInt myCharacteristic() const override;
    long myLogCardinality() const override;
    bool IamCommutative() const override  {return IsCommutative(myBaseRingValue);} // BUG??? this could fail to recognize commutativity
    bool3 IamIntegralDomain3(bool QuickMode) const override  { if (QuickMode) return uncertain3; return bool3(IsPrime(myReducingIdeal)); /*return IsPrime3(myReducingIdeal);*/ }
    bool IamTrueGCDDomain() const override  {return false;} // we don't know how to compute GCDs
    bool IamField() const override  {return IsMaximal(myReducingIdeal);}
    bool IamFiniteField() const override;
    bool IamExact() const override  {return IsExact(myReprRing);}
    ConstRefRingElem myZero() const override  {return *myZeroPtr;}
    ConstRefRingElem myOne() const override  {return *myOnePtr;}
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& N) const override;
    RingElemRawPtr myNew(ConstRawPtr rawt) const override;
    void myDelete(RawPtr rawx) const override;                      // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;                           // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;                  // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;               // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;                   // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                                // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;                  // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                           // true iff x is invertible
    bool myIsZeroDivisor(ConstRawPtr rawx) const override;
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = gcd(x,y) if TrueGCDDomain; o/w error.
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;   // lhs = x^n, n>1, x not -1,0,1
    void mySymbols(std::vector<symbol>& SymList) const override  {myReprRing->mySymbols(SymList);}
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;              // out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;                            ///< x^n may be printed without parentheses
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;                     ///< first character of x printed is a minus sign
    void myOutputSelf(std::ostream& out) const override;                            // out << R
    void myOutputSelfLong(std::ostream& out) const override;                        // out << R (descr)
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                        // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;          // OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                                 // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                  // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                             // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                   // true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                  // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                    // false iff x overflows
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;              // x == y
  //???  void convert(string&, ConstRawPtr) const override;

    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))
    bool myImageLiesInSubfield(const RingHom& phi) const override;

    void myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const override;  ///< lhs = x^N (non-triv, N large); default gives error

  public: // functions that every QuotientRing must implement
    RingElem myCanonicalRepr(ConstRawPtr rawx) const override;
    void myReduction(RawPtr rawimage, ConstRawPtr rawarg) const override;
    RingHom myInducedHomCtor(const RingHom& InducingHom) const override; // NB domain(InducingHom) == myReprRing, and ker(InducingHom) contains myReducingIdeal.

  private:
    void myReduction(RawPtr rawx) const; // implementation detail


  private:
    class InducedHomImpl: public RingHomInducedBase
    {
    public:
      InducedHomImpl(const QuotientRing& RmodI, const RingHom& InducingHom);
      // Default copy ctor & assignment disabled in RingHomBase.
      // Default dtor works fine.
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return IsPartial(myInducingHom); }
    };

  private:
    class IdealImpl: public IdealBase
    {
    public:
      friend class GeneralQuotientRingImpl;
      IdealImpl(const QuotientRing& R, const vector<RingElem>& gens);
      ~IdealImpl();
      virtual IdealImpl* myClone() const override;
//???    virtual void swap(ideal& other);

    public: // functions every ideal must have
      const QuotientRing& myRing() const override  { return myR; } // covariant return type!
      bool IamZero() const override;
      bool IamOne() const override;
      bool myAssignMaximalFlag(bool b) const override;  // virtual + default impl
      bool myAssignPrimeFlag(bool b) const override;    // virtual + default impl

      void myReduceMod(RingElemRawPtr rawx) const override; // r elem of R, where I is ideal in R
      bool IhaveElem(RingElemConstRawPtr rawx) const override;
      void myAdd(const ideal&) override;
      void myMul(const ideal&) override;
      void myIntersect(const ideal&) override;
      void myColon(const ideal&) override;
      bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const override; // return true iff quot exists & is unique; if true then lhs = num/den modulo I

      const std::vector<RingElem>& myGens() const override  { return myGensValue; }
      const std::vector<RingElem>& myTidyGens(const CpuTimeLimit&) const override  { return myTidyGensValue; }
//??? use default??      virtual void myOutputSelf(OpenMathOutput&) const override;
    protected: // functions every ideal must implement
      void myTestIsMaximal() const override;
      void myTestIsPrimary() const override;
      void myTestIsPrime() const override;
      void myTestIsRadical() const override;

    private:
      QuotientRing myR;
      vector<RingElem> myGensValue;
      vector<RingElem> myTidyGensValue;
      ideal myPreimageIdeal;
      void SetTidyGensFromPreimageIdeal(); // implementation detail
      static const IdealImpl* ourGetPtr(const ideal& I);
    };
  };


  QuotientRingBase::QuotientRingBase(const ring& R, const ideal& I):
      myBaseRingValue(R),
      myDefiningIdeal(I),
      myReprRing(IsQuotientRing(R)?ReprRing(R):R),
      //      myReducingIdeal(IsQuotientRing(R)?(I+ReducingIdeal(R)):I)
      myReducingIdeal(IsQuotientRing(R)?(I+ideal(zero(R))):I)
  {
    CoCoA_ASSERT(RingOf(I) == R);
  }


  QuotientRing::QuotientRing(const QuotientRingBase* RingPtr):
      ring(RingPtr)
  {}


  QuotientRing::QuotientRing(const ring& R):
      ring(QuotientRingPtr(R, CoCoA_ERROR_CONTEXT))
  {}


  QuotientingHomImpl::QuotientingHomImpl(const QuotientRing& RmodI):
      RingHomBase(BaseRing(RmodI), RmodI)
  {}


  void QuotientingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  { QuotientRingPtr(myCodomain)->myReduction(rawimage, rawarg); }


  void QuotientingHomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << " canonical quotient";
  }


  GeneralQuotientRingImpl::GeneralQuotientRingImpl(const ring& R, const ideal& I):
      QuotientRingBase(R, I)
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    // NOTE: currently this impl does not make its own versions of 0 and 1:
    //       it simply aliases those of myReprRing, and relabels them as belonging to this ring.
    myZeroPtr.reset(new RingElemAlias(ring(this), raw(zero(myReprRing))));
    myOnePtr.reset(new RingElemAlias(ring(this), raw(one(myReprRing))));
    myRefCountZero();
  }


  BigInt GeneralQuotientRingImpl::myCharacteristic() const
  {
    // Handle the simple case of Z/I
    if (IsZZ(myReprRing))
    {
      if (IsZero(myReducingIdeal))  return BigInt(0);
      return ConvertTo<BigInt>(TidyGens(myReducingIdeal)[0]);
    }
    if (IsOne(myReducingIdeal))  { return BigInt(1); } // trivial ring has char=1
    if (IsPolyRing(myBaseRingValue) && IsField(CoeffRing(myBaseRingValue)))
    {
      return characteristic(myBaseRingValue);
    }
    CoCoA_THROW_ERROR1(ERR::NYI);
    return characteristic(myBaseRingValue); // just to keep compiler quiet
  }


  long GeneralQuotientRingImpl::myLogCardinality() const
  {
    if (!IamFiniteField())  return 0;
    if (IsZZ(myReprRing))  return 1;
    if (IsPolyRing(myBaseRingValue))
      return CoeffRing(myBaseRingValue)->myLogCardinality() * len(QuotientBasis(myReducingIdeal)); // SLUG SLUG SLUG!!!
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return -1;
  }


  bool GeneralQuotientRingImpl::IamFiniteField() const
  {
    if (IsZZ(myBaseRingValue))  { return IsProbPrime(myCharacteristic()); }  // mildly buggy ??? IsPrime ???

    // Currently we recognise only k[x]/I where k is finite field
    if (!IsPolyRing(myReprRing)) return false;
    if (!IsField(CoeffRing(myReprRing)))
      CoCoA_THROW_ERROR1(ERR::NYI);
    return IsFiniteField(CoeffRing(myReprRing)) &&
           IsZeroDim(myReducingIdeal) &&
           IsMaximal(myReducingIdeal);
  }


  RingElem GeneralQuotientRingImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    RingElem ans(myReprRing);
    myReprRing->myAssign(raw(ans), rawx);
    return ans;
  }


  void GeneralQuotientRingImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    myReprRing->myAssign(rawimage, rawarg);
    myReducingIdeal->myReduceMod(rawimage);
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew() const
  { return myReprRing->myNew(); }


  RingElemRawPtr GeneralQuotientRingImpl::myNew(const MachineInt& n) const
  {
    RingElemRawPtr rawans = myReprRing->myNew(n);
    myReduction(rawans);
    return rawans;
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew(const BigInt& N) const
  {
    RingElemRawPtr rawans = myReprRing->myNew(N);
    myReduction(rawans);
    return rawans;
  }


  RingElemRawPtr GeneralQuotientRingImpl::myNew(ConstRawPtr rawCopyMe) const
  { return myReprRing->myNew(rawCopyMe); }


  void GeneralQuotientRingImpl::myDelete(RawPtr rawx) const
  { myReprRing->myDelete(rawx); }


  void GeneralQuotientRingImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  { myReprRing->mySwap(rawx, rawy); }


  void GeneralQuotientRingImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  { myReprRing->myAssign(rawlhs, rawx); }


  void GeneralQuotientRingImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    myReprRing->myAssign(rawlhs, n);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    myReprRing->myAssign(rawlhs, N);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myAssignZero(RawPtr rawlhs) const
  { myReprRing->myAssignZero(rawlhs); }// assume zero is always reduced


  void GeneralQuotientRingImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myReprRing));
    myReprRing->myRecvTwinFloat(rawlhs, rawx);
  }


  void GeneralQuotientRingImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myReprRing->myNegate(rawlhs, rawx);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->myAdd(rawlhs, rawx, rawy);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->mySub(rawlhs, rawx, rawy);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->myMul(rawlhs, rawx, rawy);
    myReduction(rawlhs);
  }


  void GeneralQuotientRingImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    RingElem X = myCanonicalRepr(rawx);
    RingElem Y = myCanonicalRepr(rawy);
    RingElem ans(myReprRing);
    if (!myReducingIdeal->myDivMod(raw(ans), raw(X), raw(Y)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    myReduction(rawlhs, raw(ans));
  }


  bool GeneralQuotientRingImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsOne(rawy)) { myAssign(rawlhs, rawx); return true; }
    const RingElem X = myCanonicalRepr(rawx);
    const RingElem Y = myCanonicalRepr(rawy);
    RingElem ans(myReprRing);
    if (!myReducingIdeal->myDivMod(raw(ans), raw(X), raw(Y))) return false;
    myReduction(rawlhs, raw(ans));
    return true;
  }


  // We could simply use the default definition of this function.
  bool GeneralQuotientRingImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    if (myIsZero(rawx) ) return false;
    ring R(myReprRing);
    RingElem xR(R);
    R->myAssign(raw(xR), rawx);
    bool IsInv = IsOne(myReducingIdeal + ideal(xR));
    if ( !IsInv ) myReducingIdeal->myAssignMaximalFlag(false);
    return IsInv;
  }


  bool GeneralQuotientRingImpl::myIsZeroDivisor(ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) return true;
    if (IsTrue3(IsPrime3(myReducingIdeal))) return false;
    bool IsZD = !IsZero(colon(ideal(myZero()), ideal(RingElemAlias(ring(this),rawx))));
    if ( IsZD ) myReducingIdeal->myAssignPrimeFlag(false);
    return IsZD;
  }


  //  This seems hard to define in full generality; e.g. there are some GCD domains which are not euclidean.
  void GeneralQuotientRingImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
 {
   if (!IamTrueGCDDomain())
     CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
   if (IamField())  { myGcdInField(rawlhs, rawx, rawy); return; }

   //??? BUG INCOMPLETE BUG ???
   CoCoA_THROW_ERROR2(ERR::NYI, "gcd in generic quotient ring which is not a field");
 }


  void GeneralQuotientRingImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    myBinaryPower(rawlhs, rawx, n); //??? BUG/SLUG this is not always the best choice ???
  }


  void GeneralQuotientRingImpl::myOutput(std::ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams

    // THIS PRINT FORMAT IS UGLY!  BUG???
    RingElemAlias x(myReprRing, rawx); // alias for value to be printed
    if (!IsZZ(myReprRing))
    {
      out << "(" << x << ")";
      return;
    }
    if (DefaultResidueRepr()==GlobalSettings::ResidueRepr::NonNegative)
    {
      out << x;
      return;
    }
    // Special case of Z/(nZ) and using symmetric residues...
    const BigInt n = myCharacteristic();
    if (x > n/2)
      out << x-n;
    else
      out << x;
  }


  bool GeneralQuotientRingImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    if (myIsPrintedWithMinus(rawx))  return false;
    return true;  // either positive int or between parentheses
  }


  bool GeneralQuotientRingImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    if (IsZZ(myReprRing) && DefaultResidueRepr()==GlobalSettings::ResidueRepr::symmetric )
    {
      RingElemAlias x(myReprRing, rawx); // alias for value to be printed
      return x > myCharacteristic()/2;
    }
    return false;    
  }
  

  void GeneralQuotientRingImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    //    out << "QuotientRing(" << myBaseRingValue << ", " << myDefiningIdeal << ")";
    out << "RingWithID(" << myID << ", \"";
    myBaseRingValue->myOutputSelfShort(out);
    out << "/"<<myDefiningIdeal<<"\")";
  }


  void GeneralQuotientRingImpl::myOutputSelfLong(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    myOutputSelf(out);
    out <<"\n  with BaseRing  ";
    myBaseRingValue->myOutputSelfLong(out);
  }


  void GeneralQuotientRingImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "QuotientRing");
    OMOut << myBaseRingValue;     // superfluous since the ideal already contains the ring info...
    OMOut << myDefiningIdeal;
    OMOut->mySendApplyEnd();
  }


  void GeneralQuotientRingImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  { myReprRing->myOutput_OM(OMOut, rawx); }


  bool GeneralQuotientRingImpl::myIsZero(ConstRawPtr rawx) const
  { return myReprRing->myIsZero(rawx); }


  bool GeneralQuotientRingImpl::myIsOne(ConstRawPtr rawx) const
  { return myReprRing->myIsOne(rawx) || IsOne(myReducingIdeal); }


  bool GeneralQuotientRingImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    if (myReprRing->myIsMinusOne(rawx) || IsOne(myReducingIdeal)) return true;
    RingElem tmp(myReprRing);
    myReprRing->myAdd(raw(tmp), raw(myOne()), rawx);
    return IsElem(tmp, myReducingIdeal);
  }


  bool GeneralQuotientRingImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  { return myReprRing->myIsInteger(N, rawx); }


  bool GeneralQuotientRingImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  { return myReprRing->myIsRational(Q, rawx); }


  bool GeneralQuotientRingImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  { return myReprRing->myIsDouble(d, rawx); }


  bool GeneralQuotientRingImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  { return myReprRing->myIsEqual(rawx, rawy); }


  ideal GeneralQuotientRingImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  { return ideal(new IdealImpl(QuotientRing(this), gens)); }


  void GeneralQuotientRingImpl::myReduction(RawPtr rawx) const
  { myReducingIdeal->myReduceMod(rawx); }


  //---------------------------------------------------------------------------
  // Functions to do with homomorphisms


  RingHom GeneralQuotientRingImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    CoCoA_ASSERT(IsInKer(myDefiningIdeal, InducingHom));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  RingHom GeneralQuotientRingImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  { return myInducedHomCtor(phi(theta(myQuotientingHomCtor()))); }


  bool GeneralQuotientRingImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return IamField() || myReprRing->myImageLiesInSubfield(phi); // ***BUG BUG BUG*** wrong codomain!!!!
  }


  void GeneralQuotientRingImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  { myBinaryPower(rawlhs, rawx, N); }


  GeneralQuotientRingImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& RmodI, const RingHom& InducingHom):
      RingHomInducedBase(RmodI, InducingHom)
  { /* Compatibility already checked in InducedHom in QuotientRing.C */ }


  void GeneralQuotientRingImpl::InducedHomImpl::myApply(RawPtr rawimage, ConstRawPtr rawarg) const
  { myInducingHom->myApply(rawimage, rawarg); }


  //---------------------------------------------------------------------------
  // Functions to do with GQRIdeals

  inline const GeneralQuotientRingImpl::IdealImpl* GeneralQuotientRingImpl::IdealImpl::ourGetPtr(const ideal& I)
  { return dynamic_cast<const IdealImpl*>(I.myIdealPtr()); }


  GeneralQuotientRingImpl::IdealImpl::IdealImpl(const QuotientRing& RmodI, const vector<RingElem>& gens):
      myR(RmodI),
      myGensValue(gens),
      myPreimageIdeal(ReducingIdeal(RmodI))
  {
    vector<RingElem> OverGens;
    const long ngens = len(gens);
    for (long i=0; i < ngens; ++i)
      OverGens.push_back(RmodI->myCanonicalRepr(raw(gens[i])));
    myPreimageIdeal += ideal(ReprRing(RmodI), OverGens);
    SetTidyGensFromPreimageIdeal();
    //    IamPrime3Flag = uncertain3;   // default value
    //    IamMaximal3Flag = uncertain3; // default value
  }


  GeneralQuotientRingImpl::IdealImpl::~IdealImpl()
  {}


  GeneralQuotientRingImpl::IdealImpl* GeneralQuotientRingImpl::IdealImpl::myClone() const
  { return new IdealImpl(*this); }


  // NB the forwarded call to myPreimageIdeal does all the consistency checking
  bool GeneralQuotientRingImpl::IdealImpl::myAssignPrimeFlag(bool b) const
  {
    myPreimageIdeal->myAssignPrimeFlag(b);
    IamPrime3Flag = b;
    return b;
  }


  // NB the forwarded call to myPreimageIdeal does all the consistency checking
  bool GeneralQuotientRingImpl::IdealImpl::myAssignMaximalFlag(bool b) const
  {
    myPreimageIdeal->myAssignMaximalFlag(b);
    IamMaximal3Flag = b;
    return b;
  }


  bool GeneralQuotientRingImpl::IdealImpl::IamZero() const
  { return myPreimageIdeal == ReducingIdeal(myR); }


  bool GeneralQuotientRingImpl::IdealImpl::IamOne() const
  { return IsOne(myPreimageIdeal); }


  void GeneralQuotientRingImpl::IdealImpl::myTestIsMaximal() const
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
  //???    myAssignMaximalFlag(IsMaximal(myPreimageIdeal));
  }


  void GeneralQuotientRingImpl::IdealImpl::myTestIsPrimary() const
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
  //???    myAssignMaximalFlag(IsMaximal(myPreimageIdeal));
  }


  void GeneralQuotientRingImpl::IdealImpl::myTestIsPrime() const
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
    //???    myAssignPrimeFlag(IsPrime(myPreimageIdeal));
  }


  void GeneralQuotientRingImpl::IdealImpl::myTestIsRadical() const
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
  //???    myAssignMaximalFlag(IsMaximal(myPreimageIdeal));
  }


  void GeneralQuotientRingImpl::IdealImpl::myReduceMod(RawPtr rawx) const
  {
    if (IamZero()) return;
    myPreimageIdeal->myReduceMod(rawx);
  }


  bool GeneralQuotientRingImpl::IdealImpl::IhaveElem(RingElemConstRawPtr rawx) const
  {
    return myPreimageIdeal->IhaveElem(rawx);
  }


  void GeneralQuotientRingImpl::IdealImpl::myAdd(const ideal& J)
  {
    // ANNA: implement clever insert (skip 0s, deal with 1s)
    myGensValue.insert(myGensValue.end(), gens(J).begin(), gens(J).end());
    myPreimageIdeal += ourGetPtr(J)->myPreimageIdeal;
    SetTidyGensFromPreimageIdeal();
    IamPrime3Flag = ourGetPtr(J)->IamPrime3Flag; ///??? uncertain3?
    IamMaximal3Flag = ourGetPtr(J)->IamMaximal3Flag; ///??? uncertain3?
  }


  void GeneralQuotientRingImpl::IdealImpl::myMul(const ideal& J)
  {
    if (IamZero()) return;
    //    CoCoA_THROW_ERROR1(ERR::NYI);
    vector<RingElem> tmpV; //???tmpV.reserve(len(myGens)*len(gens(J))); // BUG  OVERFLOW????
    for (const RingElem& f: myGens())
      for (const RingElem& g: gens(J))
        tmpV.push_back(f*g);
    swap(tmpV, myGensValue);
  }


  void GeneralQuotientRingImpl::IdealImpl::myIntersect(const ideal& J)
  {
    MakeUnique(myPreimageIdeal)->myIntersect(ourGetPtr(J)->myPreimageIdeal);
    SetTidyGensFromPreimageIdeal();
    myGensValue = myTidyGensValue;
    IamPrime3Flag = ourGetPtr(J)->IamPrime3Flag; ///??? uncertain3?
    IamMaximal3Flag = ourGetPtr(J)->IamMaximal3Flag; ///??? uncertain3?
  }


  void GeneralQuotientRingImpl::IdealImpl::myColon(const ideal& J)
  {
    MakeUnique(myPreimageIdeal)->myColon(ourGetPtr(J)->myPreimageIdeal);
    SetTidyGensFromPreimageIdeal();
    myGensValue = myTidyGensValue;
    IamPrime3Flag = ourGetPtr(J)->IamPrime3Flag; ///??? uncertain3?
    IamMaximal3Flag = ourGetPtr(J)->IamMaximal3Flag; ///??? uncertain3?
    //??? CANNOT WORK OUT HOW TO USE THE ALGORITHM transform OR for_each :-(
    //transform(myPreimageIdeal->gens().begin(), myPreimageIdeal->gens().end(), back_inserter(myGens), myR->myReduction);
  }


  bool GeneralQuotientRingImpl::IdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
    const QuotientRing& RmodI = myR;
    const ring& R = ReprRing(RmodI);
    const RingElem N = RmodI->myCanonicalRepr(rawnum);
    const RingElem D = RmodI->myCanonicalRepr(rawden);
    RingElem ans(R);
    if (!myPreimageIdeal->myDivMod(raw(ans), raw(N), raw(D))) return false;
    RmodI->myReduction(rawlhs, raw(ans));
    return true;
  }


  void GeneralQuotientRingImpl::IdealImpl::SetTidyGensFromPreimageIdeal()
  {
    vector<RingElem> NewTidyGens;
    const vector<RingElem>& OverGens = TidyGens(myPreimageIdeal);
    const long ngens = len(OverGens);
    for (long i=0; i < ngens; ++i)
    {
      RingElem g(myR);
      myR->myReduction(raw(g), raw(OverGens[i]));
      if (!IsZero(g))
        NewTidyGens.push_back(g);
    }
    myTidyGensValue = NewTidyGens;
  }


  //---------------------------------------------------------------------------

  QuotientRing NewQuotientRing(const ring& R, const ideal& I)
  {
    if (R != RingOf(I))  CoCoA_THROW_ERROR1(ERR::IdealNotInRing);
    if (IsOne(I))  CoCoA_THROW_ERROR2(ERR::BadQuotRing, "ideal is <1>");
//???    if (IsZero(I)) CoCoA_THROW_ERROR1(ERR::BadQuotRing); ???
    if (!IsZero(I) && IsZZ(R))
    {
      // Try to make a RingFp first as it is fastest.
      if (IsGoodForRingFp(I)) return NewRingFp(I);
      // Couldn't make a RingFp, so try to make a RingFpDouble.
      if (IsGoodForRingFpDouble(I)) return NewRingFpDouble(I);
      // Couldn't make RingFp or RingFpDouble (char not prime or too big),
      // so fall through to (relatively slow) generic implementation.
    }
    return QuotientRing(new GeneralQuotientRingImpl(R, I));
  }


  QuotientRing NewQuotientRing(const ring& R, const std::string& str)
  { return NewQuotientRing(R, ideal(RingElems(R, str))); }
  

  QuotientRing NewQuotientRing(const ring& R, const std::vector<std::string>& L)
  {
    std::vector<RingElem> gensI;
    for (long i=0; i<len(L); ++i)
      gensI.push_back(RingElem(R, L[i]));
    return NewQuotientRing(R, ideal(R, gensI));
  }
  

  QuotientRing NewZZmod(const MachineInt& n)
  { return NewZZmod(BigInt(n)); }


  QuotientRing NewZZmod(const BigInt& N)
  {
    if (N == 1 || N == -1)  CoCoA_THROW_ERROR1(ERR::BadQuotRing);
    if (IsZero(N))  CoCoA_THROW_ERROR1(ERR::BadQuotRing); //??? disallow N==0???
    return NewQuotientRing(RingZZ(),
                           ideal(RingElem(RingZZ(), N)));
  }


  RingHom QuotientingHom(const QuotientRing& RmodI)
  { return RmodI->myQuotientingHomCtor(); }


  RingHom InducedHom(const QuotientRing& RmodI, const RingHom& phi)
  {
    if (domain(phi) != BaseRing(RmodI))
      CoCoA_THROW_ERROR2(ERR::BadInducingHom, "InducedHom(QuotientRing,hom)");

    if (!IsInKer(DefiningIdeal(RmodI), phi))
      CoCoA_THROW_ERROR2(ERR::BadInducingHomKer, "InducedHom(QuotientRing,hom)");

    return RmodI->myInducedHomCtor(phi);
  }


  RingElem CanonicalRepr(ConstRefRingElem r)
  {
    if (!IsQuotientRing(owner(r)))  CoCoA_THROW_ERROR1(ERR::NotElemQuotientRing);
    return QuotientRingPtr(owner(r))->myCanonicalRepr(raw(r));
  }


} // end of namespace CoCoA
