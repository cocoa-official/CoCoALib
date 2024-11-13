#ifndef CoCoA_ring_H
#define CoCoA_ring_H

//   Copyright (c)  2005-2012,2014  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/SmartPtrIRC.H"
#include "CoCoA/assert.H"
#include "CoCoA/bool3.H"

#include "gmp.h"

#include <iosfwd>
// using std::ostream;
#include <utility>
// using std::swap;
#include <vector>
// using std::vector;


namespace CoCoA
{

  class bool3;          // fwd decl -- defined in bool3.H
  class OpenMathOutput; // fwd decl -- defined in OpenMath.H
  class OpenMathInput;  // fwd decl -- defined in OpenMath.H
  class symbol;         // fwd decl -- defined in symbol.H
  class MachineInt;     // fwd decl -- defined in MachineInt.H
  class ideal;          // fwd decl -- defined in ideal.H
  class RingHom;        // fwd decl -- defined in RingHom.H
  class BigInt;         // fwd decl -- defined in BigInt.H
  class BigRat;         // fwd decl -- defined in BigRat.H

  class RingBase; // fwd decl for SmartPtrIRC

  class ring
  {
  public:
    ring();  // default ctor sets ring to RingZZ()
    explicit ring(const RingBase* RingPtr): mySmartPtr(RingPtr) {}
    // Default copy ctor works fine.
    // Default assignment works fine.
    const RingBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    const RingBase* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions IsRingFp, AsRingFp, etc.
  private: // data members
    SmartPtrIRC<const RingBase> mySmartPtr;
  };


  class RingElemConstRawPtr
  {
  public:
    explicit RingElemConstRawPtr(void const* ptr): myPtr(ptr) {}
    // default copy ctor, assignment, and dtor are fine
    void const* myRawPtr() const { return myPtr; }
  private: // data members
    void const* myPtr;
  };

  class RingElemRawPtr: public RingElemConstRawPtr
  {
    // data members are inherited from RingElemConstRawPtr
  public:
    explicit RingElemRawPtr(void* ptr): RingElemConstRawPtr(ptr) {}
    // default copy ctor, assignment, and dtor are fine
    void* myRawPtr() const { return const_cast<void*>(RingElemConstRawPtr::myRawPtr()); }
  };


  class RingElemAlias; // fwd decl, defined 3 lines below
  typedef const RingElemAlias& ConstRefRingElem;

  class RingElemAlias
  {
  protected:   // data members
    ring myR;
    RingElemRawPtr myValuePtr;  // VALUE NOT OWNED BY ME -- deliberately not RingElemConstRawPtr (see doc)

  public:
    RingElemAlias(const ring& R, RingElemConstRawPtr rawx);
    RingElemAlias(const RingElemAlias&) = default;
    // default dtor works OK  -- deliberately NOT virtual (see documentation)
  public: // disable assignment
    RingElemAlias& operator=(const RingElemAlias& rhs) = delete;
  protected:
    const ring& myOwner() const { return myR; }
    const RingElemConstRawPtr& myRawPtr() const { return myValuePtr; } // NB implicit conversion to RingElemConstRawPtr
    // friend accessor functions (with non-member-fn syntax)
    friend const RingElemConstRawPtr& raw(ConstRefRingElem x);
    friend const ring& owner(ConstRefRingElem x);
  };


  class RingElem: public RingElemAlias
  {
    // Data members inherited from ConstRefRingElem, as is accessor myOwner.
    // NB the value pointed to by myValuePtr is OWNED by the RingElem object.

  protected:
    const RingElemRawPtr& myRawPtr() const { return myValuePtr; }
    // friend accessor functions (with non-member-fn syntax)
    friend const RingElemRawPtr& raw(RingElem& x);

  public:
    RingElem(); // default ctor initializes to 0 in RingZZ()
    explicit RingElem(const ring& R);
    RingElem(const ring& R, const MachineInt& n);
    RingElem(const ring& R, const mpz_t N);
    RingElem(const ring& R, const BigInt& N);
    RingElem(const ring& R, const mpq_t N);
    RingElem(const ring& R, const BigRat& Q);
    RingElem(const ring& R, const RingElemRawPtr rawx);
    RingElem(const ring& R, const symbol& s);
    RingElem(const ring& R, const std::string& s);
    RingElem(const ring& R, ConstRefRingElem rhs); // maps rhs into R
    RingElem(const RingElem& copy); // have to define this, default version is unsuitable, cannot be explicit as it is needed in function return
    RingElem(RingElem&& copy) /*noexcept*/; // move ctor (cannot be noexcept because allocs)
    RingElem(ConstRefRingElem copy); // becomes rather a hassle if made explicit
    ~RingElem();  ///< -- deliberately NOT virtual (see documentation)
    RingElem& operator=(const RingElem& rhs);
    RingElem& operator=(RingElem&& rhs);
    RingElem& operator=(ConstRefRingElem rhs);
    RingElem& operator=(const MachineInt& n);
    RingElem& operator=(const BigInt& N);
    RingElem& operator=(const BigRat& Q);
    friend void swap(RingElem& x, RingElem& y) noexcept; // inline defn below
  };


  // This is a bit like an auto_ptr to a RingElem; it can surrender ownership.
  class AutoRingElem
  {
  private: // data members
    const ring myR;
    RingElemRawPtr myValuePtr;  // either nullptr or VALUE IS OWNED BY ME
  public:
    AutoRingElem(const ring& R, RingElemRawPtr rawx);
    ~AutoRingElem();//inline, see below
    friend RingElemRawPtr raw(AutoRingElem& x);
    friend RingElemRawPtr release(AutoRingElem& x);
  public: // disable assignment and copy ctor
    AutoRingElem(const AutoRingElem&) = delete;
    AutoRingElem& operator=(const AutoRingElem&) = delete;
  };


  class RingBase: protected IntrusiveReferenceCount   // abstract class
  {
    friend class SmartPtrIRC<const RingBase>; // Morally "friend ring", so it can alter reference count.

  protected: // data member
    friend long RingID(const ring& R);
    static long NewRingID();
    const long myID; // Every concrete ring has a fixed distinct (numerical) ID.

  protected:
    RingBase();             // See inline functions below.
    virtual ~RingBase();    // Destructor of a base class must be virtual.
  public: // Disable assignment and copy ctor.
    RingBase(const RingBase&) = delete;
    RingBase& operator=(const RingBase&) = delete;
  public:
    typedef RingElemConstRawPtr ConstRawPtr;  // for brevity
    typedef RingElemRawPtr RawPtr;            // for brevity

    // Functions which every ring must implement.
    //===== ring: ...
    virtual const ring& myBaseRing() const = 0;
    virtual BigInt myCharacteristic() const = 0;
    virtual long myLogCardinality() const;    ///< Default: 0
    virtual void mySymbols(std::vector<symbol>& SymList) const; ///< append symbols in ring to SymList
    //===== ring: queries
    virtual bool IamCommutative() const = 0;
    virtual bool3 IamIntegralDomain3(bool QuickMode) const = 0;
    virtual bool IamTrueGCDDomain() const;    ///< Default: !IamField()
    virtual bool IamOrderedDomain() const;    ///< Default: false
    virtual bool IamField() const = 0;
    virtual bool IamFiniteField() const = 0;
    virtual bool IamExact() const = 0;
    //===== ring: printing
    virtual void myOutputSelf(std::ostream& out) const = 0;  ///< out << R
    virtual void myOutputSelfShort(std::ostream& out) const; ///< default impl -- used for finite fields when coeffring of polyring
    virtual void myOutputSelfLong(std::ostream& out) const;  ///< describe
    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const = 0; ///< OMOut << R
    //========== member functions for RingElem ==========
    //===== RingElem: memory
    virtual ConstRefRingElem myZero() const = 0;
    virtual ConstRefRingElem myOne() const = 0;
    virtual RingElem myFieldGen() const; /*default throws error -- def'd only for finite fields */
    virtual RingElemRawPtr myNew() const = 0;
    virtual RingElemRawPtr myNew(const MachineInt& n) const = 0;
    virtual RingElemRawPtr myNew(const BigInt& N) const = 0;
    virtual RingElemRawPtr myNew(const BigRat& Q) const;               ///< Default: map num & den then divide
    virtual RingElemRawPtr myNew(const symbol& s) const;
    virtual RingElemRawPtr myNew(ConstRefRingElem rhs) const;          ///< maps rhs into R
    virtual RingElemRawPtr myNew(ConstRawPtr rawx) const = 0;
    virtual void myDelete(RawPtr rawx) const = 0;                      ///< destroys x (incl all resources)
    //===== RingElem: assignments
    virtual void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const = 0;  ///< lhs = x
    virtual void myAssign(RawPtr rawlhs, const MachineInt& n) const = 0; ///< lhs = n
    virtual void myAssign(RawPtr rawlhs, const BigInt& N) const = 0;   ///< lhs = N
    virtual void myAssign(RawPtr rawlhs, const BigRat& Q) const;       ///< lhs = Q, by default map num & den then divide
    virtual void myAssignZero(RawPtr rawlhs) const = 0;                ///< lhs = 0
    virtual void mySwap(RawPtr rawx, RawPtr rawy) const /*?noexcept?*/ = 0;  ///< swap(x, y)
    virtual void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const = 0; ///< Default: throws error
    //===== RingElem: operations
    virtual void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const = 0;  ///< lhs = -x
    virtual void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const = 0; ///< lhs = x+y
    virtual void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const = 0; ///< lhs = x-y
    virtual void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const = 0; ///< lhs = x*y
    virtual void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const = 0; ///< lhs = x/y
    virtual void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;     ///< lhs = gcd(x,y) if TrueGCDDomain;
    virtual void myLcm(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;     ///< lhs = lcm(x,y) if TrueGCDDomain;
    virtual void myGcdQuot(RawPtr rawlhs, RawPtr rawxquot, RawPtr rawyquot, ConstRawPtr rawx, ConstRawPtr rawy) const; ///< lhs = gcd(x,y)  xquot = x/gcd, yquot = y/gcd  if TrueGCDDomain;
////    virtual void myExgcd(RawPtr rawlhs, RawPtr rawxcofac, RawPtr rawycofac, ConstRawPtr rawx, ConstRawPtr rawy) const; ///< lhs = gcd(x,y) = xcofac*x + ycofac*y  if TrueGCDDomain;
    virtual void myNormalizeFrac(RawPtr num, RawPtr den) const;
    virtual void myNormalizeFracNoGcd(RawPtr num, RawPtr den) const;

    void myPower(RawPtr rawlhs, ConstRawPtr rawx, long n) const;                     ///< lhs = x^n, n>=0
    void myPower(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const;            ///< lhs = x^N, N>=0
    virtual void mySquare(RawPtr rawlhs, ConstRawPtr rawx) const;                    ///< lhs = x^2, default impl
    virtual RingElem mySymbolValue(const symbol& sym) const = 0;                     ///< returns the RingElem corresponding to sym
    //===== RingElem: printing
    virtual void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const = 0;        ///< OMOut << x
    virtual void myOutput(std::ostream& out, ConstRawPtr rawx) const = 0;            ///< out << x
    virtual bool myIsPrintAtom(ConstRawPtr rawx) const;                              ///< x^n may be printed without parentheses
    virtual bool myIsPrintedWithMinus(ConstRawPtr /*rawx*/) const { return false; }  ///< first character of x printed is a minus sign. Default false
    //===== RingElem: queries
    virtual bool myIsZero(ConstRawPtr rawx) const = 0;                               ///< x == 0
    virtual bool myIsOne(ConstRawPtr rawx) const = 0;                                ///< x == 1
    virtual bool myIsMinusOne(ConstRawPtr rawx) const;                               ///< x == -1
    virtual bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const = 0;            ///< x == y
    virtual bool myIsInteger(BigInt& N, ConstRawPtr rawx) const = 0;                 ///< true iff x is integer
    virtual bool myIsRational(BigRat& Q, ConstRawPtr rawx) const = 0;                    ///< true iff x is rational
    virtual bool myIsDouble(double& d, ConstRawPtr rawx) const;                      ///< false iff x overflows
    virtual bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const = 0;///< lhs = x/y, if divisible
    virtual bool myIsZeroDivisor(ConstRawPtr rawx) const;
    virtual bool myIsInvertible(ConstRawPtr rawx) const;                             ///< true iff x is invertible
    virtual bool myIsIrred(ConstRawPtr rawx) const;                                  ///< true iff x is irreducible (if defined)
    virtual bool myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const;   ///< lhs += y*z, result says whether lhs == 0.
    virtual bool myIsZeroAddMul(RawPtr rawlhs, RawPtr rawtmp, ConstRawPtr rawy, ConstRawPtr rawz) const;   ///< lhs += y*z, result says whether lhs == 0.
    //===== only for arithmetically ordered rings
    virtual int myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const;                     ///< result is <0, =0, >0 according as x<y, x=y, x>y
    virtual int myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const;                  ///< equiv to myCmp(abs(x),abs(y))
    virtual int mySign(ConstRawPtr rawx) const;                                      ///< -1,0,+1 according as x <0,=0,>0
    virtual BigInt myFloor(ConstRawPtr rawx) const;                                  ///< return largest integer N such that N <= x
    virtual BigInt myCeil(ConstRawPtr rawx) const;                                   ///< return smallest integer N such that N >= x
    virtual BigInt myNearestInt(ConstRawPtr rawx) const;                             ///< return integer N such that abs(N-x) <= 1/2
    //===== ideal
    virtual ideal myIdealCtor(const std::vector<RingElem>& gens) const = 0;
    //===== RingHom
    virtual RingHom myCompose(const RingHom& phi, const RingHom& theta) const = 0;   ///< phi(theta(...))
    virtual bool myImageLiesInSubfield(const RingHom& phi) const = 0;                ///< true iff image set lies in subfield

  protected: // These are almost an implementation detail
    virtual void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const = 0; ///< lhs = x^n, n>1, x not -1,0,1
    virtual void myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const;///< lhs = x^N, N big, x not -1,0,1; default gives error
    void mySequentialPower(RawPtr rawlhs, ConstRawPtr rawx, long n) const;           ///< lhs = x^n (assumes n >= 0)
    void myBinaryPower(RawPtr rawlhs, ConstRawPtr rawx, long n) const;               ///< lhs = x^n (assumes n > 0)
    void myBinaryPower(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const;      ///< lhs = x^N (assumes N > 0)
    void myGcdInField(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const;      ///< lhs = gcd(x, y) in a field, result is 0 or 1
  private:   // This is an implementation detail
    void myBinaryPowerLoop(RawPtr rawlhs, ConstRawPtr rawx, long n) const;
  };


  BigInt characteristic(const ring& R);
  long LogCardinality(const ring& R);
  std::vector<symbol> symbols(const ring& R);
  bool IsPID(const ring& R);
  bool3 IsPID3(const ring& R);

  inline void swap(RingElem& x, RingElem& y)  noexcept
  { std::swap(x.myR, y.myR);  std::swap(x.myValuePtr, y.myValuePtr); }

  // Some syntactic sugar for arithmetic on RingElems
  bool operator==(ConstRefRingElem, ConstRefRingElem);
  bool operator!=(ConstRefRingElem, ConstRefRingElem);
  RingElem operator-(ConstRefRingElem);
  RingElem operator+(ConstRefRingElem, ConstRefRingElem);
  RingElem operator-(ConstRefRingElem, ConstRefRingElem);
  RingElem operator*(ConstRefRingElem, ConstRefRingElem);
  RingElem operator/(ConstRefRingElem, ConstRefRingElem);
  RingElem gcd(ConstRefRingElem, ConstRefRingElem);
  RingElem lcm(ConstRefRingElem, ConstRefRingElem);
  bool IsCoprime(ConstRefRingElem, ConstRefRingElem);
////  RingElem exgcd(RingElem& xcoeff, RingElem& ycoeff, ConstRefRingElem x, ConstRefRingElem y);
  void GcdQuot(RingElem& gcd, RingElem& xquot, RingElem& yquot, ConstRefRingElem x, ConstRefRingElem y);
  RingElem& operator+=(RingElem&, ConstRefRingElem);
  RingElem& operator-=(RingElem&, ConstRefRingElem);
  RingElem& operator*=(RingElem&, ConstRefRingElem);
  RingElem& operator/=(RingElem&, ConstRefRingElem);
  std::ostream& operator<<(std::ostream&, ConstRefRingElem);
//??? see io.H  std::ostream& operator<<(std::ostream&, const std::vector<RingElem>&);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, ConstRefRingElem);
  bool IsZero(ConstRefRingElem);
  bool IsOne(ConstRefRingElem);
  bool IsMinusOne(ConstRefRingElem);
  bool IsInteger(BigInt& N, ConstRefRingElem x);
  bool IsRational(BigRat& Q, ConstRefRingElem x);
  bool IsDouble(double& d, ConstRefRingElem x);
  bool IsInvertible(ConstRefRingElem);
  bool IsZeroDivisor(ConstRefRingElem);
  bool IsIrred(ConstRefRingElem);
  bool IsDivisible(ConstRefRingElem, ConstRefRingElem); // gives error if testing in a field
  bool IsDivisible_AllowFields(ConstRefRingElem, ConstRefRingElem); // fields are allowed
  bool IsDivisible(RingElem& lhs, ConstRefRingElem, ConstRefRingElem);
  bool IsDivisible_AllowFields(RingElem& lhs, ConstRefRingElem, ConstRefRingElem);
  RingElem power(ConstRefRingElem x, const MachineInt& n); ///< NB exponent may be negative
  RingElem power(ConstRefRingElem x, const BigInt& N);         ///< NB exponent may be negative
  ConstRefRingElem zero(const ring& R);
  ConstRefRingElem one(const ring& R);
  std::ostream& operator<<(std::ostream& out, const ring& R);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ring& R);

  RingElem FieldGen(const ring& Fq);    ///< mult generator, only for finite fields
  bool IsPthPower(ConstRefRingElem x);  ///< only in Fp and Fp[x,y,...]
  RingElem PthRoot(ConstRefRingElem x); ///< only in Fp and Fp[x,y,...]

  RingElem radical(ConstRefRingElem x); /// squarefree "factor support"
  
  RingElem binomial(RingElem x, const MachineInt& n);
  RingElem binomial(RingElem x, const BigInt& N);

  // Syntactic sugar for operations on RingElems in an ordered domain.
  int sign(ConstRefRingElem x);
  RingElem abs(ConstRefRingElem x);
  BigInt floor(ConstRefRingElem x);
  BigInt ceil(ConstRefRingElem x);
  BigInt NearestInt(ConstRefRingElem x);
  int CmpDouble(ConstRefRingElem x, double z);
  int cmp(ConstRefRingElem x, ConstRefRingElem y);
//???  int cmp(const RingElem& x, const RingElem& y);
  int CmpAbs(ConstRefRingElem x, ConstRefRingElem y);
  bool operator<(ConstRefRingElem x, ConstRefRingElem y);
  bool operator<=(ConstRefRingElem x, ConstRefRingElem y);
  bool operator>(ConstRefRingElem x, ConstRefRingElem y);
  bool operator>=(ConstRefRingElem x, ConstRefRingElem y);

  // More syntactic sugar: arithmetic between RingElems and MachineInts
  bool operator==(ConstRefRingElem, const MachineInt&);
  bool operator!=(ConstRefRingElem, const MachineInt&);
  RingElem operator+(ConstRefRingElem, const MachineInt&);
  RingElem operator-(ConstRefRingElem, const MachineInt&);
  RingElem operator*(ConstRefRingElem, const MachineInt&);
  RingElem operator/(ConstRefRingElem, const MachineInt&);
  RingElem gcd(ConstRefRingElem, const MachineInt&);
  RingElem lcm(ConstRefRingElem, const MachineInt&);
  bool IsCoprime(ConstRefRingElem, const MachineInt&);
  bool operator==(const MachineInt&, ConstRefRingElem);
  bool operator!=(const MachineInt&, ConstRefRingElem);
  RingElem operator+(const MachineInt&, ConstRefRingElem);
  RingElem operator-(const MachineInt&, ConstRefRingElem);
  RingElem operator*(const MachineInt&, ConstRefRingElem);
  RingElem operator/(const MachineInt&, ConstRefRingElem);
  RingElem gcd(const MachineInt&, ConstRefRingElem);
  RingElem lcm(const MachineInt&, ConstRefRingElem);
  bool IsCoprime(const MachineInt&, ConstRefRingElem);
  RingElem& operator+=(RingElem&, const MachineInt&);
  RingElem& operator-=(RingElem&, const MachineInt&);
  RingElem& operator*=(RingElem&, const MachineInt&);
  RingElem& operator/=(RingElem&, const MachineInt&);
  bool IsDivisible(const MachineInt&, ConstRefRingElem);
  bool IsDivisible_AllowFields(const MachineInt&, ConstRefRingElem);
  bool IsDivisible(ConstRefRingElem, const MachineInt&);
  bool IsDivisible_AllowFields(ConstRefRingElem, const MachineInt&);
  bool IsDivisible(RingElem& lhs, const MachineInt&, ConstRefRingElem);
  bool IsDivisible_AllowFields(RingElem& lhs, const MachineInt&, ConstRefRingElem);
  bool IsDivisible(RingElem& lhs, ConstRefRingElem, const MachineInt&);
  bool IsDivisible_AllowFields(RingElem& lhs, ConstRefRingElem, const MachineInt&);

  // Ordered domain operations between RingElems and MachineInts
  int cmp(ConstRefRingElem x, const MachineInt& n);
  int cmp(const MachineInt& n, ConstRefRingElem y);
  bool operator<(ConstRefRingElem x, const MachineInt& y);
  bool operator<=(ConstRefRingElem x, const MachineInt& y);
  bool operator>(ConstRefRingElem x, const MachineInt& y);
  bool operator>=(ConstRefRingElem x, const MachineInt& y);
  bool operator<(const MachineInt& x, ConstRefRingElem y);
  bool operator<=(const MachineInt& x, ConstRefRingElem y);
  bool operator>(const MachineInt& x, ConstRefRingElem y);
  bool operator>=(const MachineInt& x, ConstRefRingElem y);

  // More syntactic sugar: arithmetic between RingElems and BigInts
  bool operator==(ConstRefRingElem, const BigInt&);
  bool operator!=(ConstRefRingElem, const BigInt&);
  RingElem operator+(ConstRefRingElem, const BigInt&);
  RingElem operator-(ConstRefRingElem, const BigInt&);
  RingElem operator*(ConstRefRingElem, const BigInt&);
  RingElem operator/(ConstRefRingElem, const BigInt&);
  RingElem gcd(ConstRefRingElem, const BigInt&);
  RingElem lcm(ConstRefRingElem, const BigInt&);
  bool IsCoprime(ConstRefRingElem, const BigInt&);

  bool operator==(const BigInt&, ConstRefRingElem);
  bool operator!=(const BigInt&, ConstRefRingElem);
  RingElem operator+(const BigInt&, ConstRefRingElem);
  RingElem operator-(const BigInt&, ConstRefRingElem);
  RingElem operator*(const BigInt&, ConstRefRingElem);
  RingElem operator/(const BigInt&, ConstRefRingElem);
  RingElem gcd(const BigInt&, ConstRefRingElem);
  RingElem lcm(const BigInt&, ConstRefRingElem);
  bool IsCoprime(const BigInt&, ConstRefRingElem);
  RingElem& operator+=(RingElem&, const BigInt&);
  RingElem& operator-=(RingElem&, const BigInt&);
  RingElem& operator*=(RingElem&, const BigInt&);
  RingElem& operator/=(RingElem&, const BigInt&);
  bool IsDivisible(const BigInt& N, ConstRefRingElem);
  bool IsDivisible_AllowFields(const BigInt& N, ConstRefRingElem);
  bool IsDivisible(ConstRefRingElem, const BigInt& N);
  bool IsDivisible_AllowFields(ConstRefRingElem, const BigInt& N);
  bool IsDivisible(RingElem& lhs, const BigInt& N, ConstRefRingElem);
  bool IsDivisible_AllowFields(RingElem& lhs, const BigInt& N, ConstRefRingElem);
  bool IsDivisible(RingElem& lhs, ConstRefRingElem, const BigInt& N);
  bool IsDivisible_AllowFields(RingElem& lhs, ConstRefRingElem, const BigInt& N);

  // Ordered domain operations on RingElems and BigInts
  int cmp(ConstRefRingElem x, const BigInt& N);
  int cmp(const BigInt& N, ConstRefRingElem y);
  bool operator<(ConstRefRingElem x, const BigInt& y);
  bool operator<=(ConstRefRingElem x, const BigInt& y);
  bool operator>(ConstRefRingElem x, const BigInt& y);
  bool operator>=(ConstRefRingElem x, const BigInt& y);
  bool operator<(const BigInt& x, ConstRefRingElem y);
  bool operator<=(const BigInt& x, ConstRefRingElem y);
  bool operator>(const BigInt& x, ConstRefRingElem y);
  bool operator>=(const BigInt& x, ConstRefRingElem y);

  // More syntactic sugar: arithmetic between RingElems and BigRats
  bool operator==(ConstRefRingElem, const BigRat&);
  bool operator!=(ConstRefRingElem, const BigRat&);
  RingElem operator+(ConstRefRingElem, const BigRat&);
  RingElem operator-(ConstRefRingElem, const BigRat&);
  RingElem operator*(ConstRefRingElem, const BigRat&);
  RingElem operator/(ConstRefRingElem, const BigRat&);
///???  RingElem gcd(ConstRefRingElem, const BigRat&);
  bool operator==(const BigRat&, ConstRefRingElem);
  bool operator!=(const BigRat&, ConstRefRingElem);
  RingElem operator+(const BigRat&, ConstRefRingElem);
  RingElem operator-(const BigRat&, ConstRefRingElem);
  RingElem operator*(const BigRat&, ConstRefRingElem);
  RingElem operator/(const BigRat&, ConstRefRingElem);
///???  RingElem gcd(const BigRat&, ConstRefRingElem);
  RingElem& operator+=(RingElem&, const BigRat&);
  RingElem& operator-=(RingElem&, const BigRat&);
  RingElem& operator*=(RingElem&, const BigRat&);
  RingElem& operator/=(RingElem&, const BigRat&);
//   bool IsDivisible(const BigRat& N, ConstRefRingElem);
//   bool IsDivisible(ConstRefRingElem, const BigRat& N);

  // Ordered domain operations on RingElems and BigRats
  int cmp(ConstRefRingElem x, const BigRat& N);
  int cmp(const BigRat& N, ConstRefRingElem y);
  bool operator<(ConstRefRingElem x, const BigRat& y);
  bool operator<=(ConstRefRingElem x, const BigRat& y);
  bool operator>(ConstRefRingElem x, const BigRat& y);
  bool operator>=(ConstRefRingElem x, const BigRat& y);
  bool operator<(const BigRat& x, ConstRefRingElem y);
  bool operator<=(const BigRat& x, ConstRefRingElem y);
  bool operator>(const BigRat& x, ConstRefRingElem y);
  bool operator>=(const BigRat& x, ConstRefRingElem y);


  //---------------------------------------------------------------------------
  // inline functions on rings -- the order of appearance might be important!


  inline bool operator==(const ring& R1, const ring& R2)
  { return R1.myRawPtr() == R2.myRawPtr(); }
//    return R1.mySmartPtr == R2.mySmartPtr;  // equivalent formulation

  inline bool operator!=(const ring& R1, const ring& R2)
  { return !(R1 == R2); }


  //----------------------------------------------------------------------
  // Inline functions.  The "owner" and "raw" accessor functions must come first.


  inline AutoRingElem::AutoRingElem(const ring& R, RingElemRawPtr rawx):
      myR(R),
      myValuePtr(rawx)
  {}


  inline AutoRingElem::~AutoRingElem()
  {
    if (myValuePtr.myRawPtr())
      myR->myDelete(myValuePtr);
  }

  inline RingElemRawPtr raw(AutoRingElem& x)
  {
    CoCoA_ASSERT(x.myValuePtr.myRawPtr() != nullptr);
    return x.myValuePtr;
  }

  inline RingElemRawPtr release(AutoRingElem& x)
  {
    CoCoA_ASSERT(x.myValuePtr.myRawPtr() != nullptr);
    RingElemRawPtr rawans = x.myValuePtr;
    x.myValuePtr = RingElemRawPtr(nullptr);
    return rawans;
  }


  inline RingElemAlias::RingElemAlias(const ring& R, RingElemConstRawPtr rawx):
      myR(R),
      myValuePtr(const_cast<void*>(rawx.myRawPtr()))
  {}


  inline const ring& BaseRing(const ring& R)
  { return R->myBaseRing(); }


//??? BEWARE reference to temporary???? e.g. const ring& R = owner(a+b);
  inline const ring& owner(ConstRefRingElem x)
  { return x.myOwner(); }

  inline const RingElemConstRawPtr& raw(ConstRefRingElem x)
  { return x.myRawPtr(); }

  inline const RingElemRawPtr& raw(RingElem& x)
  { return x.myRawPtr(); }

  inline long RingID(const ring& R)
  { return R->myID; }

  inline RingBase::RingBase():
    IntrusiveReferenceCount(),
    myID(NewRingID())
  {}


  inline RingBase::~RingBase()
  {}


  inline RingElem::RingElem(const ring& R):
      RingElemAlias(R, R->myNew())
  {}


  inline RingElem::RingElem(const ring& R, RingElemRawPtr rawx):
      RingElemAlias(R, rawx)
  {}


  inline RingElem::RingElem(const RingElem& copy):
      RingElemAlias(owner(copy), owner(copy)->myNew(raw(copy)))
  {}

  inline RingElem::RingElem(RingElem&& copy) /*noexcept*/ :
      RingElemAlias(owner(copy), owner(copy)->myNew())
  { std::swap(myValuePtr, copy.myValuePtr); }
//  { std::swap(myRawPtr(),copy.myRawPtr()); }
//  { owner(copy)->mySwap(myRawPtr(),copy.myRawPtr()); }


  inline RingElem::RingElem(ConstRefRingElem copy):
      RingElemAlias(owner(copy), owner(copy)->myNew(raw(copy)))
  {}


  inline ConstRefRingElem zero(const ring& R)
  { return R->myZero(); }

  inline ConstRefRingElem one(const ring& R)
  { return R->myOne(); }

  inline bool IsTrueGCDDomain(const ring& R)
  { return R->IamTrueGCDDomain(); }

  inline bool IsCommutative(const ring& R)
  { return R->IamCommutative(); }

  inline bool IsIntegralDomain(const ring& R)
  { return IsTrue3(R->IamIntegralDomain3(false)); }

  inline bool3 IsIntegralDomain3(const ring& R)
  { return R->IamIntegralDomain3(true); }

  inline bool IsField(const ring& R)
  { return R->IamField(); }

  inline bool IsFiniteField(const ring& R)
  { return R->IamFiniteField(); }

  inline bool IsExact(const ring& R)
  { return R->IamExact(); }

  inline bool IsOrderedDomain(const ring& R)
  { return R->IamOrderedDomain(); }

  //??? bool contains(const RingBase& big_ring, const RingBase& little_ring);

  // Define the not-equals functions here:

  inline bool operator!=(ConstRefRingElem x, ConstRefRingElem y)
  { return !(x == y); }

  inline bool operator!=(ConstRefRingElem r, const MachineInt& n)
  { return !(r == n); }

  inline bool operator!=(const MachineInt& n, ConstRefRingElem r)
  { return !(r == n); }

  inline bool operator!=(ConstRefRingElem r, const BigInt& N)
  { return !(r == N); }

  inline bool operator!=(const BigInt& N, ConstRefRingElem r)
  { return !(r == N); }

} // end of namespace CoCoA

#endif
