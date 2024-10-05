#ifndef CoCoA_OpenMath_H
#define CoCoA_OpenMath_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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

#include <iosfwd>
// using std::ostream;
#include <string>
// using std::string


namespace CoCoA
{

  class BigInt;      // fwd decl, defined in BigInt.H
  class MachineInt;  // fwd decl, defined in MachineInt.H

  class OpenMathSymbol
  {
  public:
    OpenMathSymbol(); // fills symbol with a default "unset" value
    OpenMathSymbol(const char* const cd, const char* const name);
    OpenMathSymbol(const std::string& cd, const std::string& name);
    // Default copy ctor, assignment and dtor are fine.

    void myOutputSelf(std::ostream& out) const;

  private: // data members
    std::string myCD;
    std::string myName;
    friend const std::string& CD(const OpenMathSymbol& s);   // accessor
    friend const std::string& name(const OpenMathSymbol& s); // accessor
  };

  std::ostream& operator<<(std::ostream& out, const OpenMathSymbol& oms);

  //---------------------------------------------------------------------------
  // Output functions etc.

  class OpenMathOutputBase; // fwd decl, defined below.
  class OpenMathOutput
  {
  public:
    explicit OpenMathOutput(OpenMathOutputBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    OpenMathOutputBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
  private: // data members
    SmartPtrIRC<OpenMathOutputBase> mySmartPtr;
  };


//  Removed next line to permit auto conversion from bool to bool3
//  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const MachineInt& n);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, int n);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, unsigned int n);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, long n);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, unsigned long n);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const OpenMathSymbol& s);


  class OpenMathOutputBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<OpenMathOutputBase>; // morally "friend OpenMathOutput", so it can alter reference count
  protected:
    OpenMathOutputBase(); // just sets ref count to zero.
    virtual ~OpenMathOutputBase();
  public:
    virtual void mySend(const MachineInt& n) = 0;
    virtual void mySend(const BigInt& N) = 0;
    virtual void mySend(const OpenMathSymbol& s) = 0;
    virtual void mySendSymbol(const char* const cd, const char* const name) { mySend(OpenMathSymbol(cd,name)); };
    virtual void mySendApplyStart() = 0;
    virtual void mySendApplyEnd() = 0;
    virtual void mySendObjectStart() = 0;
    virtual void mySendObjectEnd() = 0;

    //void OutputAttribute() = 0;//??????
  };


  //---------------------------------------------------------------------------
  // Input functions etc.

  enum class OpenMathTag { Obj, Apply, Int, Sym, Attr, ReachedEOF };

  class OpenMathInputBase;  // fwd decl, defined below.
  class OpenMathInput
  {
  public:
    explicit OpenMathInput(OpenMathInputBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    OpenMathInputBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
  private: // data members
    SmartPtrIRC<OpenMathInputBase> mySmartPtr;
  };

  class OpenMathInputBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<OpenMathInputBase>; // morally "friend OpenMathInput", so it can alter reference count
  protected:
    OpenMathInputBase(); // just sets ref count to zero.
    virtual ~OpenMathInputBase();
  public:
    virtual void advance() = 0;
    virtual OpenMathTag myCurrTag() = 0;
    virtual long NumDescendants() const = 0;
    virtual bool myRecv(long & n) = 0;
//????    virtual bool myRecv(unsigned long & n) = 0; // UNSIGNED SHORT INTEGER
    virtual bool myRecv(BigInt& N) = 0;
    virtual bool myRecv(OpenMathSymbol& s) = 0; // result is false if node is of wrong type

//       virtual void InputApplyStart() = 0;
//       virtual void InputApplyEnd() = 0;
//       virtual void InputObjectStart() = 0;
//       virtual void InputObjectEnd() = 0;

    //void InputAttribute() = 0;//??????
  };

  OpenMathInput& operator>>(OpenMathInput& OMIn, long& n);
  OpenMathInput& operator>>(OpenMathInput& OMIn, OpenMathSymbol& s); //????

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/OpenMath.H,v 1.13 2022/02/18 14:11:47 abbott Exp $
// $Log: OpenMath.H,v $
// Revision 1.13  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.12  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.11  2021/01/07 14:18:58  abbott
// Summary: Corrected copyright
//
// Revision 1.10  2012/05/04 16:59:43  abbott
// For future compatibility with auto conversion from bool to 3-way bool,
// removed sending of MachineInt; instead added operators for sending
// (unsigned) int and (unsigned) long [but no other integral types for now].
//
// Revision 1.9  2012/02/02 15:55:48  abbott
// Replaced include of MachineInt by fwd decl.
//
// Revision 1.8  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.7  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2011/08/12 15:48:31  abbott
// Added virtual mem fns for BigInt type.
//
// Revision 1.5  2011/03/11 14:49:08  abbott
// Changed size_t into long.
//
// Revision 1.4  2010/03/18 13:53:48  abbott
// Minor rationalization to OpenMath implementation: moved op<<(OMOut,ZZ) to ZZ files.
//
// Revision 1.3  2008/12/16 21:10:32  abbott
// Replaced the various output fns for different sort of machine integers by a
// single one for MachineInt.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/27 13:04:35  cocoa
// Added explicit OpenMath output operators for all integer types.
//
// Revision 1.4  2006/11/23 17:10:52  cocoa
// -- changed: OpenMathOutput and OpenMathInput are now a class (instead of typedef)
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


#endif
