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

#endif
