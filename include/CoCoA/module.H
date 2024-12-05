#ifndef CoCoA_module_H
#define CoCoA_module_H

//   Copyright (c)  2003-2009  John Abbott,  Anna M. Bigatti

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


// Header file for the classes module and ModuleBase, and ModuleElem


//#include "CoCoA/ring.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/SmartPtrIRC.H"

#include <iosfwd>
// using std::ostream;


namespace CoCoA
{

  class BigInt;              // fwd decl -- defined in BigInt.H
  class RingElemAlias;       // fwd decl -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;    // fwd decl -- defined in ring.H
  class OpenMathOutput;      // fwd decl -- defined in OpenMath.H
  class BigRat;              // fwd decl -- defined in BigRat.H
  class RingElemConstRawPtr; // fwd decl -- defined in ring.H
  class ring;                // fwd decl -- defined in ring.H

  union ModuleRawPtr
  {
    void* ptr;
  };

  class ModuleBase; // fwd decl for SmartPtrIRC

  //// class module --------------------------------------------------
  class module
  {
  public:
    explicit module(const ModuleBase* ModulePtr): mySmartPtr(ModulePtr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    const ModuleBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    const ModuleBase* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions Is(..), As(..), etc.
    //bool operator==(const module& M) const {return mySmartPtr==M.mySmartPtr;}
    bool operator==(const module& M) const;
  private: // data members
    SmartPtrIRC<const ModuleBase> mySmartPtr;
  };


  //// class ModuleElem --------------------------------------------------
  class ModuleElem
  {
  public:
    explicit ModuleElem(const module& M);
//???    ModuleElem(const module& M, ConstRawPtr rawvalue);
    ModuleElem(const ModuleElem& copy);
    ~ModuleElem();
    ModuleElem& operator=(const ModuleElem& rhs);
    ConstRefRingElem operator[](long pos) const;  ///< v[pos] (READ ONLY), defined only if element of a FreeModule

    friend ModuleRawPtr& raw(ModuleElem& v);
    friend const ModuleRawPtr& raw(const ModuleElem& v);
    friend const module& owner(const ModuleElem& v);

  private: // data members
    module myM;
    ModuleRawPtr myValue;
  };


  //// class ModuleBase --------------------------------------------------
  class ModuleBase: protected IntrusiveReferenceCount  // abstract class
  {
    friend class SmartPtrIRC<const ModuleBase>; // morally "friend module", so it can alter reference count

  protected:
    ModuleBase(); // just zeroes myRefCount, inline (see below)
    virtual ~ModuleBase(); // this class is abstract!!
  public:
    typedef ModuleRawPtr RawPtr;
    typedef const ModuleRawPtr& ConstRawPtr;

    virtual const ring& myRing() const = 0;
    virtual bool IamZero() const = 0;

    virtual const ModuleElem& myZero() const = 0;
    virtual void myNew(RawPtr& rawv) const = 0;
    virtual void myNew(RawPtr& rawv, ConstRawPtr rawv1) const = 0;
    virtual void myDelete(RawPtr& rawv) const = 0;                                    // destroys x (incl all resources)
    virtual void mySwap(RawPtr& rawv, RawPtr& raww) const = 0;                        // swap(x, y);
    virtual void myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const = 0;                // lhs = v;
    virtual void myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const = 0;                // lhs = -v
    virtual void myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const = 0; // lhs = v+w;
    virtual void mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const = 0; // lhs = v-w;

    virtual void myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const = 0; // lhs = r*v;
/////?????    virtual void myMul(RawPtr& rawlhs, ConstRawPtr rawv, RingElemConstRawPtr rawx) const = 0; // lhs = v*r;
    virtual void myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const = 0; // lhs = v/r;  FUNNY ARG ORDER!!!
    virtual void myOutput(std::ostream& out, ConstRawPtr rawv) const = 0;             // out << v
    virtual void myOutputSelf(std::ostream& out) const = 0;                           // out << M
    virtual void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const = 0;         // OMOut << v
    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const = 0;                       // OMOut << M
    virtual bool myIsZero(ConstRawPtr rawv) const = 0;                                // v == 0
//???    virtual bool IsZeroAddMul(RawPtr& lhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.
    virtual bool myIsEqual(ConstRawPtr rawv, ConstRawPtr raww) const = 0;
  ///  virtual void convert(string&, RawPtr) const = 0;
  };


  // Some syntactic sugar for arithmetic on ModuleElems
  ModuleElem operator-(const ModuleElem&);
  ModuleElem operator+(const ModuleElem&, const ModuleElem&);
  ModuleElem operator-(const ModuleElem&, const ModuleElem&);
  ModuleElem operator*(ConstRefRingElem, const ModuleElem&);
  ModuleElem operator*(const ModuleElem&, ConstRefRingElem);
  ModuleElem operator/(const ModuleElem&, ConstRefRingElem);

  ModuleElem& operator+=(ModuleElem&, const ModuleElem&);
  ModuleElem& operator-=(ModuleElem&, const ModuleElem&);
  ModuleElem& operator*=(ModuleElem&, ConstRefRingElem);
  ModuleElem& operator/=(ModuleElem&, ConstRefRingElem);

  // Arith between ModuleElems and MachineInts
  ModuleElem operator*(const MachineInt&, const ModuleElem&);
  ModuleElem operator*(const ModuleElem&, const MachineInt&);
  ModuleElem operator/(const ModuleElem&, const MachineInt&);

  ModuleElem& operator*=(ModuleElem&, const MachineInt&);
  ModuleElem& operator/=(ModuleElem&, const MachineInt&);

  // Arith between ModuleElems and BigInts
  ModuleElem operator*(const BigInt&, const ModuleElem&);
  ModuleElem operator*(const ModuleElem&, const BigInt&);
  ModuleElem operator/(const ModuleElem&, const BigInt&);

  ModuleElem& operator*=(ModuleElem&, const BigInt&);
  ModuleElem& operator/=(ModuleElem&, const BigInt&);

  // Arith between ModuleElems and BigRats
  ModuleElem operator*(const BigRat&, const ModuleElem&);
  ModuleElem operator*(const ModuleElem&, const BigRat&);
  ModuleElem operator/(const ModuleElem&, const BigRat&);

  ModuleElem& operator*=(ModuleElem&, const BigRat&);
  ModuleElem& operator/=(ModuleElem&, const BigRat&);


  std::ostream& operator<<(std::ostream&, const ModuleElem& v);
  std::ostream& operator<<(std::ostream&, const module& M);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ModuleElem& v);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const module& M);
  bool IsZero(const ModuleElem&);
  //  const ModuleElem& zero(const module& M);
  bool operator==(const ModuleElem& v1, const ModuleElem& v2);
  bool operator!=(const ModuleElem& v1, const ModuleElem& v2);


  /////////////////////////////////////////////////////////////////////////////
  // Inline functions


  inline const ring& RingOf(const module& M)
  {return M->myRing();}


  inline bool IsZero(const module& M)
  {return M->IamZero();}


  inline const ModuleElem& zero(const module& M)
  {return M->myZero();}


  inline const module& owner(const ModuleElem& v)
  {return v.myM;}


  inline ModuleRawPtr& raw(ModuleElem& v)
  {return v.myValue;}


  inline const ModuleRawPtr& raw(const ModuleElem& v)
  {return v.myValue;}


  inline bool operator!=(const module& M1, const module M2)
  {return !(M1 == M2);}


  inline ModuleBase::ModuleBase():
    IntrusiveReferenceCount()
  {
    myRefCountInc(); // NOT ZERO, otherwise FreeModule ctor dies ???
  }


  //??? bool contains(const module& big_module, const module& little_module);


} // end of namespace CoCoA

#endif
