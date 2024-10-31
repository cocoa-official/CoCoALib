#ifndef CoCoA_PPMonoid_H
#define CoCoA_PPMonoid_H

//   Copyright (c)  2001-2017  John Abbott and Anna M. Bigatti
//   Author:  2001-2010  John Abbott

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


#include "CoCoA/PPOrdering.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"
// #include "CoCoA/SmartPtrIRC.H" --- already included in PPOrdering.H

// #include <iosfwd> --- already included in PPOrdering.H
// using std::ostream;
#include <vector>
// using std::vector;


namespace CoCoA
{

  class BigInt;         // fwd decl -- defined in BigInt.H
  class DivMask;        // fwd decl -- defined in DivMask.H
  class DivMaskRule;    // fwd decl -- defined in DivMask.H
  class MachineInt;     // fwd decl -- defined in MachineInt.H
  class OpenMathOutput; // fwd decl -- defined in OpenMath.H
  class degree;         // fwd decl -- defined in degree.H
  class matrix;         // fwd decl -- defined in matrix.H

  class PPMonoidBase; // forward declaration for SmartPtrIRC

  class PPMonoid
  {
  public:
    explicit PPMonoid(const PPMonoidBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    const PPMonoidBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    bool operator==(const PPMonoid& PPM) const { return mySmartPtr==PPM.mySmartPtr; }
  private: // data members
    SmartPtrIRC<const PPMonoidBase> mySmartPtr;
  };


  // These two classes are simply wrapped void*

  class PPMonoidElemConstRawPtr
  {
  public:
    explicit PPMonoidElemConstRawPtr(void const* ptr): myPtr(const_cast<void*>(ptr)) {}
    // default copy ctor, assignment, and dtor are fine
    bool operator==(const PPMonoidElemConstRawPtr& pp) const;
    bool operator!=(const PPMonoidElemConstRawPtr& pp) const;
    void const* myRawPtr() const { return myPtr; }
  protected: // data member
    void* myPtr;
  };

  class PPMonoidElemRawPtr: public PPMonoidElemConstRawPtr
  {
    // NB data members are inherited from PPMonoidElemConstRawPtr
  public:
    explicit PPMonoidElemRawPtr(void* ptr): PPMonoidElemConstRawPtr(ptr) {}
    // default copy ctor, assignment, and dtor are fine
    void* myRawPtr() { return myPtr; }
  };


  //---------------------------------------------------------------------------
  // The next three classes mimic the trio of PPMonoidElem classes.
  class ConstRefPPMonoidElem
  {
  protected: // data members
    const PPMonoid myPPM;
    PPMonoidElemRawPtr const myPPPtr;  ///< VALUE NOT OWNED BY ME, deliberately NOT PPMonoidElemConstRawPtr

    // friend accessor functions (with non-member-fn syntax)
    friend PPMonoidElemConstRawPtr raw(const ConstRefPPMonoidElem& pp);
    friend const PPMonoid& owner(const ConstRefPPMonoidElem& pp);

  public:
    ConstRefPPMonoidElem(const PPMonoid& PPM, PPMonoidElemConstRawPtr rawpp);
    ConstRefPPMonoidElem(const ConstRefPPMonoidElem&) = default;
    // default dtor works OK -- deliberately NOT virtual (see documentation)
  public: // disable assignment
    ConstRefPPMonoidElem& operator=(const ConstRefPPMonoidElem& rhs) = delete;
  };


  class RefPPMonoidElem: public ConstRefPPMonoidElem
  {
  protected: // data members inherited from ConstRefPPMonoidElem
    // friend accessor functions (with non-member-fn syntax)
    friend PPMonoidElemRawPtr raw(RefPPMonoidElem& pp);

  public:
    RefPPMonoidElem(const PPMonoid& PPM, PPMonoidElemRawPtr rawpp);
    RefPPMonoidElem(const RefPPMonoidElem& t) = default;// to avoid compiler warning!
    // default dtor works OK  -- deliberately NOT virtual (see documentation)
    RefPPMonoidElem& operator=(const RefPPMonoidElem& rhs);
    RefPPMonoidElem& operator=(const ConstRefPPMonoidElem& rhs);
  };



//   class PPMonoidElem; // fwd decl for friend decl.
//   class PPMonoidElemCheckVecSize
//   {
//     friend class PPMonoidElem;
//     PPMonoidElemCheckVecSize() {};
//     PPMonoidElemCheckVecSize(long NumIndets, long VecSize) { if (VecSize != NumIndets) CoCoA_THROW_ERROR1(ERR::BadArraySize); };
//   };

  class PPMonoidElem: public RefPPMonoidElem
  {
    // NB data members inherited from ConstRefPPMonoidElem via RefPPMonoidElem
  public:
    explicit PPMonoidElem(const PPMonoid& Gamma);
    PPMonoidElem(const PPMonoid& Gamma, const std::vector<long>& v);
    PPMonoidElem(const PPMonoid& Gamma, const std::vector<BigInt>& v);
    PPMonoidElem(const PPMonoid& Gamma, PPMonoidElemRawPtr rawToBeOwned);
    PPMonoidElem(const PPMonoidElem& copy);
    PPMonoidElem(const ConstRefPPMonoidElem& copy);
    PPMonoidElem& operator=(const PPMonoidElem& rhs);
    PPMonoidElem& operator=(const ConstRefPPMonoidElem& rhs);
    ~PPMonoidElem();  ///< deliberately NOT virtual (see documentation)
  };


  //---------------------------------------------------------------------------
  //  class PPMonoidBase /* ABSTRACT CLASS */
  class PPMonoidBase: protected IntrusiveReferenceCount  /* ABSTRACT CLASS */
  {
    friend class SmartPtrIRC<const PPMonoidBase>; ///< morally "friend PPMonoid", so it can alter reference count

  protected:
    PPMonoidBase(const PPOrdering& ord, const std::vector<symbol>& IndetNames);
    virtual ~PPMonoidBase() {};
  public: // disable copy ctor and assignment
    PPMonoidBase(const PPMonoidBase& copy) = delete;
    PPMonoidBase& operator=(const PPMonoidBase&) = delete;

  public:
    typedef PPMonoidElemRawPtr RawPtr;           ///< just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; ///< just to save typing

    const PPOrdering& myOrdering() const { return myOrd; };
    virtual const std::vector<PPMonoidElem>& myIndets() const = 0;                ///< std::vector whose n-th entry is n-th indet as PPMonoidElem
    const std::vector<symbol>& mySymbols() const { return myIndetSymbols; }
    virtual const symbol& myIndetSymbol(long i) const;                            ///< the symbol of the i-th indet
    ConstRefPPMonoidElem mySymbolValue(const symbol& s) const;                    ///< returns the PPMonoidElem corresponding to sym

    virtual const PPMonoidElem& myOne() const = 0;
    virtual PPMonoidElemRawPtr myNew() const = 0;                                 ///< ctor from nothing
    virtual PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const = 0;    ///< ctor from another pp
    virtual PPMonoidElemRawPtr myNew(const std::vector<long>& v) const = 0;       ///< ctor from exp vector
    PPMonoidElemRawPtr myNewCheckVecSize(const std::vector<long>& v) const;       ///< ctor (with check) from exp vector
    virtual PPMonoidElemRawPtr myNew(const std::vector<BigInt>& v) const;         ///< ctor from exp vector (default goes via vector<long>)
    PPMonoidElemRawPtr myNewCheckVecSize(const std::vector<BigInt>& v) const;     ///< ctor (with check) from exp vector

    virtual void myDelete(RawPtr rawpp) const = 0;                                ///< dtor, frees pp
    virtual void mySwap(RawPtr rawpp1, RawPtr rawpp2) const = 0;                  ///< swap(pp1, pp2)
    virtual void myAssignOne(RawPtr rawpp) const = 0;                             ///< pp = 1
    virtual void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const = 0;            ///< pp = pp1
    virtual void myAssign(RawPtr rawpp, const std::vector<long>& v) const = 0;    ///< pp = v (assign from exp vector)
    virtual void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0; ///< pp = pp1*pp2
    virtual void myMulIndetPower(RawPtr rawpp, long var, long exp) const = 0;     ///< pp *= indet(var)^exp, assumes exp >= 0
    virtual void myMulIndetPower(RawPtr rawpp, long var, const BigInt& EXP) const;    ///< pp *= indet(var)^EXP, NOT PURE
    virtual void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0; ///< pp = pp1/pp2, assumes the quotient exists
    virtual void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0;///< pp = pp1/gcd(pp1,pp2)
    virtual void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0; ///< pp = gcd(pp1,pp2)
    virtual void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0; ///< pp = lcm(pp1,pp2)
    virtual void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const = 0;           ///< pp = radical(pp1)
    void myPower(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const;               ///< pp = pp1^exp, assumes exp >= 0
    void myPower(RawPtr rawpp, ConstRawPtr rawpp1, const BigInt& EXP) const;          ///< pp = pp1^EXP, assumes EXP >= 0
    virtual void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const = 0;// pp = pp1^exp (non-trivial), assumes exp >= 0
    virtual void myPowerBigExp(RawPtr rawpp, ConstRawPtr rawpp1, const BigInt& EXP) const; ///< pp = pp1^EXP (non-triv, EXP large); default gives error
    virtual void myPowerOverflowCheck(ConstRawPtr rawpp1, long exp) const = 0;    ///< throw if pp1^exp would overflow, assumes exp >= 0

    virtual bool myIsOne(ConstRawPtr rawpp) const = 0;                            ///< is pp = 1?
    virtual bool myIsIndet(long& index, ConstRawPtr rawpp) const = 0;             ///< true iff pp is an indet
    virtual bool myIsIndetPosPower(long& index, BigInt& EXP, ConstRawPtr rawpp) const; ///< true iff pp is a positive power of an indet
    virtual bool myIsIndetPosPower(long& index, long& pow, ConstRawPtr rawpp) const; ///< true iff pp is a positive power of an indet and the exponent fits into a long (no checks!)
    virtual bool myIsIndetPosPower(ConstRawPtr rawpp) const; ///< true iff pp is a positive power of an indet (no checks!)
    virtual bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0;   ///< are pp1, pp2 coprime?
    virtual bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const;         ///< is pp1 equal to pp2?
    virtual bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0; ///< is pp1 divisible by pp2?
    virtual bool myIsSqFree(ConstRawPtr rawpp) const = 0;                         ///< is pp equal to its radical?

    virtual int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0;          ///< <0, =0, >0 as pp1 < = > pp2
    //    virtual int myHomogCmp(ConstRawPtr rawt1, ConstRawPtr rawt2) const = 0;   ///< <0, =0, >0 as t1 < = > t2 assuming t1 and t2 have same multi-degree
    //    virtual int myHomogDegRevLex(ConstRawPtr rawt1, ConstRawPtr rawt2) const = 0; ///< <0, =0, >0 as t1 < = > t2 ??? degrevlex assuming t1 and t2 have same multi-degree TO BE REMOVED

    virtual long myStdDeg(ConstRawPtr rawpp) const = 0;                           ///< standard degree of pp
    virtual void myWDeg(degree& d, ConstRawPtr rawpp) const = 0;                  ///< d = grading(pp)
    virtual int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const = 0;      ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    virtual int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const = 0;      ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2) wrt the first n rows of weights
    virtual long myExponent(ConstRawPtr rawpp, long i) const = 0;                 ///< exponent of ith indet in pp
    virtual void myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long i) const = 0; ///< EXP = degree of ith indet in pp

    virtual void myExponents(std::vector<long>& v, ConstRawPtr rawpp) const = 0;  ///< get exponents, SHOULD BE vector<BigInt> ????
    virtual void myBigExponents(std::vector<BigInt>& v, ConstRawPtr rawpp) const = 0;  ///< get exponents, SHOULD BE myExponents ???
    virtual void myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const = 0;   ///< v[i] = true if i-th indet has exponent != 0
    virtual void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const = 0; ///< computes the DivMask for pp according to DivMaskImpl
    virtual void myOutputSelf(std::ostream& out) const = 0;                       ///< print value of PPMonoid
//???    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const = 0;              ///< print OpenMath value of PPMonoid
    virtual void myOutput(std::ostream& out, ConstRawPtr rawpp) const;            ///< NOT PURE!!
    virtual void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawpp) const;        ///< NOT PURE!!

  protected: // data members (common to all PPOrdering implementations)
    const PPOrdering myOrd;
    const std::vector<symbol> myIndetSymbols;
    const long myNumIndets;
    friend long NumIndets(const PPMonoid& M);
  };

  const std::vector<PPMonoidElem>& indets(const PPMonoid& PPM);      ///< std::vector whose n-th entry is n-th indet as PPMonoidElem
  const symbol& IndetSymbol(const PPMonoid& PPM, long idx);
  const std::vector<symbol>& symbols(const PPMonoid& PPM);                  ///< std::vector of the symbols in PPM
  ConstMatrixView OrdMat(const PPMonoid& M);
  ConstMatrixView GradingMat(const PPMonoid& M);

  long StdDeg(ConstRefPPMonoidElem pp);                              ///< standard degree of pp
  long deg(ConstRefPPMonoidElem pp);                                 ///< standard degree of pp
  degree wdeg(ConstRefPPMonoidElem pp);                              ///< degree according to grading
  int CmpWDeg(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);   ///< <0 =0 >0 according as wdeg(t1) < = > wdeg(t2)
  int CmpWDegPartial(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2, long n);   ///< <0 =0 >0 according as wdeg(t1) < = > wdeg(t2) wrt the first n rows of weights
  long exponent(ConstRefPPMonoidElem pp, long var);                  ///< exponent in pp of indet of index var
  BigInt BigExponent(ConstRefPPMonoidElem pp, long var);             ///< exponent in pp of indet of index var
  const std::vector<long>& exponents(std::vector<long>& v, ConstRefPPMonoidElem pp);  ///< rtn val is const ref to 1st arg
  std::vector<long> exponents(ConstRefPPMonoidElem pp);              ///< creates a new vector!!
  const std::vector<BigInt>& BigExponents(std::vector<BigInt>& v, ConstRefPPMonoidElem pp);  ///< SHOULD BE exponents? ///< rtn val is const ref to 1st arg
  std::vector<BigInt> BigExponents(ConstRefPPMonoidElem pp);         ///< creates a new vector!!
  void swap(RefPPMonoidElem pp1, RefPPMonoidElem pp2);               ///< swap(t1, t2);

  // comparisons
  bool IsOne(ConstRefPPMonoidElem pp);
  bool IsIndet(long& index, ConstRefPPMonoidElem pp);
  bool IsIndet(ConstRefPPMonoidElem pp);
  bool IsIndetPosPower(long& index, BigInt& exp, ConstRefPPMonoidElem pp);
  bool IsIndetPosPower(long& index, long& exp, ConstRefPPMonoidElem pp);
  bool IsIndetPosPower(ConstRefPPMonoidElem pp);
  int cmp(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);         ///< <0, =0, >0 as pp1 < = > pp2
  bool operator==(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2); ///< pp1 == pp2;
  bool operator!=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2); ///< pp1 != pp2;
  bool operator<(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);  ///< pp1 < pp2;
  bool operator<=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2); ///< pp1 <= pp2;
  bool operator>(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);  ///< pp1 > pp2;
  bool operator>=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2); ///< pp1 => pp2;

  std::ostream& operator<<(std::ostream& out, const PPMonoid& PPM);
  std::ostream& operator<<(std::ostream& out, ConstRefPPMonoidElem pp);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, ConstRefPPMonoidElem t);


  // arithmetic
  PPMonoidElem operator*(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2); ///< pp1*pp2;
  PPMonoidElem operator/(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2); ///< pp1/pp2;
  PPMonoidElem colon(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);     ///< pp1:pp2
///??? saturation???
  PPMonoidElem gcd(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);       ///< gcd(pp1, pp2)
  PPMonoidElem lcm(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);       ///< lcm(pp1, pp2)
  PPMonoidElem radical(ConstRefPPMonoidElem pp);                              ///< radical(pp)
  std::vector<long> IndetsIn(ConstRefPPMonoidElem pp);                        ///< list of indices of indets actually in pp
  PPMonoidElem power(ConstRefPPMonoidElem pp, long exp);                      ///< pp^exp
  PPMonoidElem power(ConstRefPPMonoidElem pp, const BigInt& EXP);             ///< pp^EXP
  void PowerOverflowCheck(ConstRefPPMonoidElem pp, long exp);                 ///< error if  pp^exp  would overflow the exponents
  PPMonoidElem root(ConstRefPPMonoidElem pp, const MachineInt& exp);          ///< pp^(1/exp)
  bool IsPower(ConstRefPPMonoidElem pp, const MachineInt& exp);               /// true iff pp = t^exp for some t
  bool IsCoprime(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);         ///< is gcd(pp1, pp2)=1?
  bool IsDivisible(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2);       ///< is pp1 divisible by pp2?
  bool IsSqFree(ConstRefPPMonoidElem pp);                                     ///< is pp = radical(pp)?

  bool IsFactorClosed(const std::vector<PPMonoidElem>& S);

  void AssignOne(RefPPMonoidElem dest);
  RefPPMonoidElem operator*=(RefPPMonoidElem dest, ConstRefPPMonoidElem t);
  RefPPMonoidElem operator/=(RefPPMonoidElem dest, ConstRefPPMonoidElem t);

  const PPMonoidElem& indet(const PPMonoid& M, long i); ///< x[i]
  const PPMonoidElem& indet(const PPMonoid& M, const BigInt& i); ///< x[i]
  PPMonoidElem IndetPower(const PPMonoid& M, long var, long exp);
  PPMonoidElem IndetPower(const PPMonoid& M, long var, const BigInt& EXP);

  PPMonoid NewPPMonoid(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  PPMonoid NewPPMonoid(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);


  //---------------------------------------------------------------------------
  // INLINE FUNCTIONS

  inline bool PPMonoidElemConstRawPtr::operator==(const PPMonoidElemConstRawPtr& pp) const
  { return myPtr == pp.myPtr; }


  inline bool PPMonoidElemConstRawPtr::operator!=(const PPMonoidElemConstRawPtr& pp) const
  { return !operator==(pp); } // weird syntax


  //---------------------------------------------------------------------------
  // inline functions on PPMonoids -- the order of appearance might be important!

  inline bool operator!=(const PPMonoid& M1, const PPMonoid& M2)
  { return !(M1 == M2); }


  //----------------------------------------------------------------------

  inline ConstRefPPMonoidElem one(const PPMonoid& PPM)
  { return PPM->myOne(); }


  inline long NumIndets(const PPMonoid& M)
  { return M->myNumIndets; }


  inline const PPOrdering& ordering(const PPMonoid& M)
  { return M->myOrdering(); }


  inline long GradingDim(const PPMonoid& M)
  { return GradingDim(ordering(M)); }


  inline bool IsStdGraded(const PPMonoid& M)
  { return IsStdGraded(ordering(M)); }


  inline const std::vector<PPMonoidElem>& indets(const PPMonoid& PPM)
  { return PPM->myIndets(); }


  inline const symbol& PPMonoidBase::myIndetSymbol(long var) const
  {
    CoCoA_ASSERT(0 <= var && var < len(myIndetSymbols));
    return myIndetSymbols[var];
  }


  inline ConstRefPPMonoidElem::ConstRefPPMonoidElem(const PPMonoid& PPM, PPMonoidElemConstRawPtr rawpp):
      myPPM(PPM), myPPPtr(const_cast<void*>(rawpp.myRawPtr()))
  {}


  inline RefPPMonoidElem::RefPPMonoidElem(const PPMonoid& PPM, PPMonoidElemRawPtr rawpp):
      ConstRefPPMonoidElem(PPM, rawpp)
  {}


  inline const PPMonoid& owner(const ConstRefPPMonoidElem& rawpp)
  { return rawpp.myPPM; }


  inline PPMonoidElemRawPtr raw(RefPPMonoidElem& rawpp)
  { return rawpp.myPPPtr; }


  inline PPMonoidElemConstRawPtr raw(const ConstRefPPMonoidElem& rawpp)
  { return rawpp.myPPPtr; }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PPMonoid.H,v 1.43 2024/07/31 07:52:46 abbott Exp $
// $Log: PPMonoid.H,v $
// Revision 1.43  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.42  2023/12/01 20:26:39  abbott
// Summary: Improved comment
//
// Revision 1.41  2023/03/09 22:33:44  abbott
// Summary: Added explicit default copy ctor (to keep compiler quiet)
//
// Revision 1.40  2022/09/14 19:15:04  abbott
// Summary: Revised interface to exponents (redmine 1694)
//
// Revision 1.39  2022/03/07 13:44:51  abbott
// Summary: Added explicit default copy ctor (to avoid compiler warning, g++ v.9)
//
// Revision 1.38  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.37  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.36  2021/10/30 17:11:00  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.35  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.34  2020/02/14 12:20:43  abbott
// Summary: Renamed indets(PP) to IndetsIn(PP)
//
// Revision 1.33  2017/12/01 17:13:38  bigatti
// -- myOrdering inline code is now within the class definition
//
// Revision 1.32  2016/11/03 12:25:24  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.31  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.30  2015/06/30 12:52:13  abbott
// Summary: Added new fns myPowerOverflowCheck and myIndetsIn
// Author: JAA
//
// Revision 1.29  2015/04/13 14:42:07  abbott
// Summary: Added myPowerOverflowCheck (1st version)
// Author: JAA
//
// Revision 1.28  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.27  2014/07/03 15:36:35  abbott
// Summary: Cleaned up impl of PPMonoids: moved myIndetSymbols & myNumIndets to base class
// Author: JAA
//
// Revision 1.26  2014/02/25 16:27:57  abbott
// Summary: Added new fn IsFactorClosed
// Author: JAA
//
// Revision 1.25  2013/07/30 14:57:19  bigatti
// -- reduced empty lines
// -- added IsStdGraded
//
// Revision 1.24  2013/05/31 09:15:43  abbott
// Improved a comment.
//
// Revision 1.23  2012/04/17 19:55:40  abbott
// Added 2 new fns: root and IsPower.
//
// Revision 1.22  2012/02/10 17:07:28  abbott
// Added pseudo-ctor mem fn myNew(vector<BigInt>), and some related fns.
// Added new fn indets.
//
// Revision 1.21  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.20  2011/06/23 16:05:55  abbott
// Improved a comment.
//
// Revision 1.19  2011/03/22 15:26:46  bigatti
// -- added myIsIndetPosPower/IsIndetPosPower
//
// Revision 1.18  2011/03/16 15:18:48  bigatti
// -- added IsIndet(pp), IsIndetPosPower(pp)
//
// Revision 1.17  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.16  2010/11/30 11:18:11  bigatti
// -- renamed IndetName --> IndetSymbol
//
// Revision 1.15  2010/11/25 09:31:01  bigatti
// -- forward declaration --> fwd decl
//
// Revision 1.14  2010/11/05 16:21:08  bigatti
// -- added ZZExponents
//
// Revision 1.13  2010/11/02 16:03:12  bigatti
// -- added indet(*, ZZ) [for CoCoA-5]
// -- indet(PPM, i) now returns "const PPMonoidElem&"
//
// Revision 1.12  2010/10/01 15:20:33  bigatti
// -- added mySymbolValue
// -- added RingElem(R, sym)
//
// Revision 1.11  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.10  2010/02/03 12:11:08  bigatti
// -- added myIsIndetPosPower with exponent variable of type long
//    (no checks for overflow)
//
// Revision 1.9  2010/02/02 16:44:31  abbott
// Added radical & IsRadical (via mem fns myRadical & myIsRadical)
// for PPMonoidElems.
//
// Revision 1.8  2009/11/26 16:12:29  bigatti
// -- commented out some includes
//
// Revision 1.7  2009/09/22 14:01:33  bigatti
// -- added myCmpWDegPartial (ugly name, I know....)
// -- cleaned up and realigned code in PPMonoid*.C files
//
// Revision 1.6  2008/10/07 12:18:29  abbott
// Added deg function, same as StdDeg.
//
// Revision 1.5  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.4  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.3  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/09/24 14:15:16  abbott
// Renamed IsIndetPower to IsIndetPosPower.
// Changed layout of some comments.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.14  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.13  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.12  2007/01/11 14:06:15  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.11  2006/12/06 17:32:24  cocoa
// -- removed #include "config.H"
//
// Revision 1.10  2006/11/24 17:26:22  cocoa
// -- reorganized includes of header files
// -- doxygen style comments
//
// Revision 1.9  2006/11/23 17:35:48  cocoa
// -- changed: PPMonoid is now a class (instead of typedef)
//
// Revision 1.8  2006/11/22 14:53:46  cocoa
// -- changed: PPMonoid defined as class (instead of typedef)
//
// Revision 1.7  2006/11/03 15:37:47  cocoa
// -- cleaned up code after testing on usage of SmartPtrIRC
//
// Revision 1.6  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.5  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/10/06 10:08:58  cocoa
// Corrected a comment.
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.10  2006/04/21 14:59:22  cocoa
// Changed return type of indet(PPM, i) to ConstRefPPMonoidElem
// instead of PPMonoidElem (which made a copy).
//
// Revision 1.9  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.8  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.7  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.6  2006/03/07 10:08:48  cocoa
// -- fixed: PPMonoidElem(const ConstRefPPMonoidElem& copy);  [added const]
//
// Revision 1.5  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.4  2006/02/01 16:56:13  cocoa
// Added some missing assignment operators for (Ref)PPMonoidElems.
//
// Revision 1.3  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.2  2005/11/17 15:43:04  cocoa
// -- added ctor PPMonoidElem(ConstRefPPMonoidElem& copy);
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.7  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.6  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.5  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.4  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.3  2005/06/23 15:42:41  cocoa
// Fixed typo in GNU fdl -- all doc/*.txt files affected.
// Minor corrections to PPMonoid (discovered while writing doc).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.14  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.13  2004/11/11 13:41:14  cocoa
// -- moved CoCoA includes to the top
//
// Revision 1.12  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
//
// Revision 1.11  2004/11/02 14:56:33  cocoa
// -- changed *Print* into *Output* (myPrint --> myOutput)
// -- changed *Var* into *Indet* (myPrintVarName --> myOutputIndetName)
// -- removed suffix "IgnoreDivMask"
// -- added myComputeDivMask
// -- improved storing of IndetNames
// -- changed ExpvElem into SmallExponent_t
//
// Revision 1.10  2004/10/29 16:11:10  cocoa
// -- new type   PPMonoidEvOv_entry;
// -- function IsDivisible had wrong semantics --> swapped arguments everywhere
//
// Revision 1.9  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namespace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.8  2004/06/30 16:46:06  cocoa
// End of day check-in: mostly tidying up, and correcting some
// overly lax access permissions.
//
// Revision 1.7  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.6  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.5  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.4  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.3  2004/01/28 15:20:36  cocoa
// Added CmpDeg function and myCmpDeg member function.
// Aligned name of data member of PPMonoid with the convention used
// for the data member of ring.
//
// Revision 1.2  2003/10/01 10:35:31  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.22  2003/09/23 16:59:33  bigatti
// - added PPMonoidDiv_entry* in RawPP
//
// Revision 1.21  2003/06/23 16:56:59  abbott
// Minor cleaning prior to public release.
//
// Revision 1.20  2003/05/14 16:38:54  abbott
// Changed name of assign1 to AssignOne.
// Added new member functions PPMonoidBase::IsOne, IsOneIgnoreDivMask,
// IsEqual, and IsEqualIgnoreDivMask (following the pattern of
// IsDivisible).
// Added PPMonoidBase::ZeroRefCount.
// Added inline definitions of IsOne and IsEqual.
//
// Revision 1.19  2003/04/29 13:59:29  abbott
// Uh oh.  Lots of changes.  Checking in prior to restructuring
// to align with new PPOrdering code.
//
// Revision 1.18  2003/01/15 15:58:39  abbott
// Checking in prior to major change.
// Minor changes:
//   coprime has become AreCoprime,
//   divides has become IsDivisible,
//   new member fn IsDivisibleIgnoreDivMask,
//   signature of deg has changed (now produces values of type degree)
//
// Revision 1.17  2002/11/15 16:24:50  abbott
// MAJOR CHANGE reflecting a new "philosophy".
// The new philosophy is that PPMonoid::elems will be "slow and safe" since
// they are now highly decoupled from the internal representation of
// polynomials (which use "order vectors").
// Eliminated HalfPP, alias and temp subclasses.
// Every RawPP now contains its own "divisibility mask" because Anna needs
// the DivMask check to be inline.
//
// Revision 1.16  2002/06/22 17:15:14  abbott
// Numerous changes to accommodate the new idea of PPOrdering.
//
// Revision 1.15  2002/03/06 11:22:38  abbott
// Added stubs for two new comparison functions (implemented via
// member functions) to allow comparison between a RawPP and a
// HalfPP, and between two HalfPPs.
//
// Revision 1.14  2002/02/08 12:07:22  abbott
// MAJOR CHANGE: added two new proxy classes "alias" and "temp" (see PP.txt)
// Numerous consequential changes.  Changed order of some declarations.
//
// Revision 1.13  2002/01/30 14:18:19  abbott
// Replaced VarName function by PrintVarName; simplifies maintaining backward
// compatibility with gcc 2.95, and may even be slightly easier to use.
// Added owner(...) and raw(..) functions for PPMonoid::elem objects.
//
// Revision 1.12  2002/01/10 17:04:39  abbott
// Revised names of two of the members of the RawPP union to
// concord with their implementations.  The new names also
// make more sense than the old ones.
//
// Revision 1.11  2002/01/09 14:32:40  abbott
// Added/corrected HalfPP type and operations.
//
// Revision 1.10  2001/12/12 18:32:23  bigatti
// - HalfMul and OtherHalfMul
//
// Revision 1.9  2001/12/07 18:17:46  abbott
// Added field cpctexp in RawPP for PP_compact values.
// The field needs a better name.
//
// Revision 1.8  2001/12/05 15:18:40  abbott
// Act first, think later...  Fixed another copy of the name of struct PPMonoid_DenseUniq.
//
// Revision 1.7  2001/12/05 15:17:46  abbott
// Corrected name of struct PPmonoid_DenseUniq_entry.
//
// Revision 1.6  2001/12/04 19:53:08  abbott
// Names changed in accordance with new coding standards.
// (public) Function NumVars is now a non-member (but calls a private member
// function to allow inheritance to work as desired).
//
// Revision 1.5  2001/11/19 15:27:24  abbott
// Added support for "compact_uniq" PP implementation:
//   one forward declaration of a struct,
//   a new field to the PPmonoid::raw_elem union.
//
// Revision 1.4  2001/11/16 18:40:00  bigatti
// added: <iostream> and std::
// for compatibility with gcc-3
//
// Revision 1.3  2001/11/05 18:13:29  abbott
// Added homogeneous PP comparison functions: i.e. assuming two PPs have the
// same (multi)degree it compares the PPs using the monoid's ordering.
//
// Anna also wanted an explicit homog_degrevlex; JAA hopes this function will
// eventually disappear, though it may be needed for efficiency reasons.
//
// The names homog_cmp and homog_degrevlex are too cumbersome.
//
// Revision 1.2  2001/10/31 20:25:56  abbott
// MAIN CHANGE:
//   PPmonoid::raw_elem is now a union of the different types it can hold
//   rather than a void*.  JAA is not sure if it is really cleaner to use
//   a union than a void* (with casting), but the compiler gave enough
//   problems to make the change worthwhile.  Furthermore, a union does
//   give some extra checks/restrictions on what can be done to raw_elem.
//   One significant difficulty was producing an l-value when using a
//   (reinterpret) cast hidden inside a file static inline function.
//   Also it is easier to debug with unions since the debugger will allow
//   field selection in unions -- reinterpret_casts were a headache.
//
// The "copy constructor" for a PPmonoid::elem was altered to require a const
// copy (previously the const was missing).
//
// Revision 1.1  2001/10/25 17:53:53  abbott
// Initial revision
//

#endif
