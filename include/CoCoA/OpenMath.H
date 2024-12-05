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

#endif
