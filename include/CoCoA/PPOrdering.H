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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PPOrdering.H,v 1.21 2022/02/18 14:11:47 abbott Exp $
// $Log: PPOrdering.H,v $
// Revision 1.21  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.20  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.19  2021/10/30 17:11:35  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.18  2021/02/22 21:09:25  abbott
// Summary: Implemented xel (redmine 1536)
//
// Revision 1.17  2017/12/01 17:14:46  bigatti
// -- just header with Copyright
//
// Revision 1.16  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.15  2015/12/01 13:11:00  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.14  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.13  2014/12/10 11:39:23  bigatti
// -- added NewMatrixOrdering with 2 args (remove with 3 args?)
//
// Revision 1.12  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.11  2014/01/28 16:44:40  abbott
// Added new fns:  IsMatrixOrdering  and  IsTermOrdering
//
// Revision 1.10  2013/07/30 15:01:07  bigatti
// -- added IamStdGraded
//
// Revision 1.9  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/09/13 15:48:36  abbott
// Changed "struct" into "class" PPOrdering (bug found by clang).
//
// Revision 1.7  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.6  2010/02/03 16:54:32  bigatti
// -- added virtual destructor
//
// Revision 1.5  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.4  2009/11/26 17:00:01  bigatti
// -- enabled operator== (used in an assert)
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.7  2006/11/27 13:45:56  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/24 17:33:26  cocoa
// -- reorganized includes of header files
//
// Revision 1.5  2006/11/22 14:50:33  cocoa
// -- changed: PPOrdering defined as class (instead of typedef)
//
// Revision 1.4  2006/11/14 17:16:20  cocoa
// -- removed commented code for reference counting
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
// Revision 1.5  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.4  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.3  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.13  2004/11/29 16:22:34  cocoa
// -- added function for computing adjoint and inverse for DenseMatrix
//    (so adjoint/inverse matrix is computed by OrdvArith and is no
//    longer needed by PPOrdering)
//
// Revision 1.12  2004/11/11 14:41:16  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.11  2004/11/05 16:44:20  cocoa
// -- deleted MatrixOrderingMod32749Impl (implemented in OrdvArith)
// -- changed C++ matrices into "matrix" over RingZ
//
// Revision 1.10  2004/11/03 17:54:44  cocoa
// -- added implementation of GetMatrix (OrdMat)
// -- added some functions for order matrices modulo 32749:
//    they will be deleted soon
//
// Revision 1.9  2004/11/02 14:49:03  cocoa
// -- new code for matrix orderings
//
// Revision 1.8  2004/10/29 16:09:21  cocoa
// -- added MatrixOrderingMod32749Impl (not tested)
// -- fixed assignment of myAdjointMatrix for MatrixOrderingImpl
//
// Revision 1.7  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//

#endif
