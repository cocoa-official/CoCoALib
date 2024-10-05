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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/module.H,v 1.21 2022/02/18 14:11:52 abbott Exp $
// $Log: module.H,v $
// Revision 1.21  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.20  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.19  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.18  2020/06/22 15:41:02  abbott
// Summary: Added commented out mem fn for right-mult of ModuleElem by RingElem
//
// Revision 1.17  2014/07/30 14:02:17  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.16  2013/07/31 09:49:22  bigatti
// -- added IsZero(module)
//
// Revision 1.15  2013/06/06 05:44:29  bigatti
// -- minor style/spaces changes
//
// Revision 1.14  2013/06/03 14:02:18  bigatti
// -- changed operator== in now defined in .C
//
// Revision 1.13  2013/01/23 14:04:49  bigatti
// -- removed doxygen includes (obsolete)
//
// Revision 1.12  2012/10/24 12:09:34  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.11  2012/10/11 14:37:02  abbott
// Removed non-const operator[], thereby removing the need for RefRingElem.
//
// Revision 1.10  2011/11/09 13:51:29  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.9  2011/08/24 10:22:45  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.8  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.7  2011/03/14 10:29:28  abbott
// Removed superfluous #include directive.
//
// Revision 1.6  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.5  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.4  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.3  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/01/15 13:33:34  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.7  2006/12/06 17:18:40  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/24 17:23:23  cocoa
// -- removed OpenMathFwd.H
//
// Revision 1.5  2006/11/23 17:37:10  cocoa
// -- changed: module is now a class (instead of typedef)
//
// Revision 1.4  2006/11/14 17:46:20  cocoa
// -- changed: reference counting in modules now uses SmartPtrIRC
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.7  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.6  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2005/11/29 13:04:47  cocoa
// -- added "const" to myCompt argument
//
// Revision 1.2  2005/11/24 16:09:38  cocoa
// -- added operator[] for ModuleElem
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.7  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.6  2004/11/11 14:04:01  cocoa
// -- minor changes for doxygen
//
// Revision 1.5  2004/11/09 16:10:43  cocoa
// -- minor changes for doxygen
//
// Revision 1.4  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.3  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.2  2004/01/28 15:35:37  cocoa
// Major update -- was very much "old style" code, and didn't
// compiler under the new organization.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.1  2003/05/30 15:10:47  abbott
// Initial revision
//
//

#endif
