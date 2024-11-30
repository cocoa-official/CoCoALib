#ifndef CoCoA_PPOrdering_H
#define CoCoA_PPOrdering_H

//   Copyright (c)  2001-2017  John Abbott and Anna M. Bigatti
//   Author:  2004-2015  John Abbott

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


#include "CoCoA/MachineInt.H"
#include "CoCoA/SmartPtrIRC.H"

#include <iosfwd>
// using std::ostream;


namespace CoCoA
{

  class OpenMathOutput; // forward declaration -- defined in OpenMath.H
  class matrix;         // forward declaration -- defined in matrix.H
  class MatrixView;     // forward declaration -- defined in matrix.H
  class ConstMatrixView;// forward declaration -- defined in matrix.H

  class PPOrderingBase; // forward declaration for SmartPtrIRC
  //---------------------------------------------------------------------------
  class PPOrdering
  {
  public:
    explicit PPOrdering(const PPOrderingBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    const PPOrderingBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    const PPOrderingBase* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions Is(..), As(..), etc.
    bool operator==(const PPOrdering& PPO) const { return mySmartPtr==PPO.mySmartPtr; }
  private: // data members
    SmartPtrIRC<const PPOrderingBase> mySmartPtr;
  };


  std::ostream& operator<<(std::ostream& out, const PPOrdering& PPO);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const PPOrdering& PPO);

  //---------------------------------------------------------------------------
  class PPOrderingBase: private IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<const PPOrderingBase>; // morally "friend PPOrdering", so it can alter ref count.

  protected:
    PPOrderingBase(long NumIndets, long GradingDim);
    virtual ~PPOrderingBase();
  public: // disable copy ctor & assignment
    PPOrderingBase(const PPOrderingBase& copy) = delete;
    PPOrderingBase& operator=(const PPOrderingBase&) = delete;
  public:
    virtual void myOutputSelf(std::ostream& out) const = 0;
    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const = 0;
    virtual const ConstMatrixView& myOrdMat() const = 0;
    virtual bool IamStdGraded() const = 0;

  protected:
    friend long NumIndets(const PPOrdering& PPO);
    friend long GradingDim(const PPOrdering& PPO);
    friend bool IsStdGraded(const PPOrdering& PPO);

    long myNumIndets;        ///< Number of indeterminates
    long myGradingDim;       ///< zero if not graded, otherwise grading is in Z^myGradingDim
  };


  //---------------------------------------------------------------------------
  // Inline non-member functions

  inline long NumIndets(const PPOrdering& PPO)
  { return PPO->myNumIndets; }


  inline long GradingDim(const PPOrdering& PPO)
  { return PPO->myGradingDim; }


  inline bool IsStdGraded(const PPOrdering& PPO)
  { return PPO->IamStdGraded(); }

  //---------------------------------------------------------------------------
  // Pseudo-constructor for PPOrdering defined by a matrix
  PPOrdering NewMatrixOrdering(const ConstMatrixView& OrderMatrix, const MachineInt& GradingDim);


  // Functions for recognizing the underlying implementation of a PP ordering.
  bool IsLex(const PPOrdering& PPO);
  bool IsXel(const PPOrdering& PPO);
  bool IsStdDegLex(const PPOrdering& PPO);
  bool IsStdDegRevLex(const PPOrdering& PPO);
  bool IsMatrixOrdering(const PPOrdering& PPO);

  bool IsTermOrdering(const PPOrdering& PPO);

  ConstMatrixView OrdMat(const PPOrdering& PPO);
  ConstMatrixView GradingMat(const PPOrdering& PPO);

  
  // The next few classes are to allow one to use "lex", "StdDegLex" or "StdDegRevLex"
  // as the ordering (without having to create the PPOrdering object).  This is just
  // syntactic sugar.
  class PPOrderingCtor
  {
  public:
    virtual ~PPOrderingCtor() {}
    virtual PPOrdering operator()(const MachineInt& NumIndets) const = 0;
  };

  struct LexCtor: public PPOrderingCtor
  {
    PPOrdering operator()(const MachineInt& NumIndets) const;
  };

  struct XelCtor: public PPOrderingCtor
  {
    PPOrdering operator()(const MachineInt& NumIndets) const;
  };

  struct StdDegLexCtor: public PPOrderingCtor
  {
    PPOrdering operator()(const MachineInt& NumIndets) const;
  };

  struct StdDegRevLexCtor: public PPOrderingCtor
  {
    PPOrdering operator()(const MachineInt& NumIndets) const;
  };

  extern LexCtor lex;
  extern XelCtor xel;
  extern StdDegLexCtor StdDegLex;
  extern StdDegRevLexCtor StdDegRevLex;


} // end of namespace CoCoA

#endif
