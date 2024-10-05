#ifndef CoCoA_OrdvArith_H
#define CoCoA_OrdvArith_H

//   Copyright (c)  2005-2008,2014  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/assert.H"
#include "CoCoA/long32or64.H"
#include "CoCoA/SmartPtrIRC.H"

#include <vector>
// using std::vector;
#include <iosfwd>
// using std::ostream;

namespace CoCoA
{

  class PPOrdering;       // forward declaration -- defined in PPOrdering.H
  class degree;           // forward declaration -- defined in degree.H
  class ConstMatrixView;  // forward declaration -- defined in matrix.H

  namespace OrdvArith
  {

    typedef unsigned long OrdvElem; // NB must be at least as big as SmallExponent_t!!
                                    //    Checked by CoCoA_ASSERT(myPackingDensity >= 1) in ctor.

    class base; // forward reference
    class reference
    {
    public:
      explicit reference(const base* ptr): mySmartPtr(ptr) {}
      // assignment disabled because SmartPtrIRC has no assignment
      const base* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
      const base* myRawPtr() const { return mySmartPtr.myRawPtr(); }
      bool operator==(const reference& ref) const { return mySmartPtr==ref.mySmartPtr; }
    private: // data members
      SmartPtrIRC<const base> mySmartPtr;
    };

    //    class base
    class base: protected IntrusiveReferenceCount
    {
      friend class SmartPtrIRC<const base>;  // morally "friend base", so it can alter reference count

    protected:
      base(long NumIndets, long OrdvWords, long NumOrdvEntries);
      virtual ~base();

    public:
      virtual void myAssignZero(OrdvElem* ordv) const;
      virtual void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const = 0;
      void myAssign(OrdvElem* dest, const OrdvElem* src) const;
      void mySwap(OrdvElem* ordv1, OrdvElem* ordv2) const;
      void myMul(OrdvElem* ordv, const OrdvElem* ordv1, const OrdvElem* ordv2) const;
      virtual void myMulIndetPower(OrdvElem* ordv, long var, long exp) const; ///< NOT PURE!!  Assumes exp >= 0
      void myDiv(OrdvElem* ordv, const OrdvElem* ordv1, const OrdvElem* ordv2) const;
      void myPower(OrdvElem* ordv, const OrdvElem* ordv1, long exp) const;  // assumes exp >= 0
      void myPowerOverflowCheck(const OrdvElem* ordv, long exp) const;  // assumes exp >= 0
      int myCmp(const OrdvElem* l_ordv, const OrdvElem* r_ordv) const noexcept;
      virtual void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const = 0;
      virtual long myExponent(const OrdvElem* ordv, long var) const = 0;
      virtual long myStdDeg(const OrdvElem* ordv) const;
      virtual void myWDeg(degree& d, const OrdvElem* ordv) const;
      int myCmpWDeg(const OrdvElem* ordv1, const OrdvElem* ordv2) const { return myCmpWDegPartial(ordv1, ordv2, myGradingDim); }
      virtual int myCmpWDegPartial(const OrdvElem* ordv1, const OrdvElem* ordv2, long GrDim) const;  // assumes GrDim >= 0
      virtual bool myIsZero(const OrdvElem* ordv) const;
      virtual bool myIsIndet(long& index, const OrdvElem* ordv) const;
      virtual void myOutputSelf(std::ostream& out) const = 0;

    protected:
      OrdvElem myOrdvGetNth(const OrdvElem* ordv, long n) const;
      void myCompress(OrdvElem* ordv, const std::vector<OrdvElem>& buffer) const;
      void myDecompress(std::vector<OrdvElem>& buffer, const OrdvElem* ordv, long NumCompts) const;

    protected: // data members
      OrdvElem myOrdvMask; // also used as max value!!
      long myPackingDensity;   // because C++ shift ops work strangely if shift length = word length
      long myBitsPerOrdvEntry;
      long myNumIndets;        ///< Number of indeterminates
      long myGradingDim;       ///< zero if not graded, otherwise grading is in Z^myGradingDim
      long myOrdvWords;        ///< byte-size-of-an-order-vector/sizeof(OrdvElem)
      long myOrdvWordsForCmp;  ///< might be smaller than myOrdvWords
      mutable std::vector<OrdvElem> myOrdvBuffer;  // only for single threaded  *****!!!!!NOT THREADSAFE!!!!!*****
      mutable std::vector<long> myExpvBuffer;      // only for single threaded  *****!!!!!NOT THREADSAFE!!!!!*****
  
      friend long OrdvWords(const reference& OA) noexcept;
    };


    class LexImpl: public base
    {
    public:
      LexImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    class XelImpl: public base
    {
    public:
      XelImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    class StdDegLexImpl: public base
    {
    public:
      StdDegLexImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myStdDeg(const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    class StdDegRevLexImpl: public base
    {
    public:
      StdDegRevLexImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myStdDeg(const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    class StdDegRevLexImpl2: public base
    {
    public:
      StdDegRevLexImpl2(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myStdDeg(const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    };


    class MatrixOrderingImpl: public base // !!!BUG!!! only partially implemented
    {
    public:
      MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    private:
      int myOrderMatrixDet;
      std::vector< std::vector<int> > myOrderMatrix;
      std::vector< std::vector<int> > myAdjointOrderMatrix;
    };


    class MatrixOrderingMod32749Impl: public base
    {
    public:
      MatrixOrderingMod32749Impl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    private:
      std::vector< std::vector<int> > myOrderMatrix;
      std::vector< std::vector<int> > myInverseOrderMatrix;
    private:
      static constexpr int ourMODULUS = 32749; // largest prime < 2^15
    };


    class MatrixOrdering32bitImpl: public base
    {
    public:
      MatrixOrdering32bitImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override;
      long myExponent(const OrdvElem* ordv, long var) const override;
      void myOutputSelf(std::ostream& out) const override;
    private:
      using  INT = long;   // int32 ???
#ifdef CoCoA_32BIT_LONG
      static constexpr INT ourMODULUS = 32749; // largest prime < 2^15
#else
      static constexpr INT ourMODULUS = 2147483587; // largest prime < 2^31
#endif
      std::vector< std::vector<INT> > myOrderMatrix;
      std::vector< std::vector<INT> > myInverseOrderMatrix;
    };



    //---------------------------------------------------------------------------
    // INLINE FUNCTIONS


    inline long OrdvWords(const reference& OA) noexcept
    {
      return OA->myOrdvWords;
    }


    inline void base::myMul(OrdvElem* ordv, const OrdvElem* ordv1, const OrdvElem* ordv2) const
    {
#ifdef CoCoA_DEBUG
      // if debugging is active, check for overflow
      std::vector<OrdvElem> expand1(myNumIndets);
      std::vector<OrdvElem> expand2(myNumIndets);
      myDecompress(expand1, ordv1, myNumIndets);
      myDecompress(expand2, ordv2, myNumIndets);
      for (long i=0; i < myNumIndets; ++i)
      {
        CoCoA_ASSERT("Exponent Overflow" && expand1[i] <= myOrdvMask-expand2[i]);
      }
#endif
      for (long i=0; i < myOrdvWords; ++i)
        ordv[i] = ordv1[i]+ordv2[i];
    }


    inline void base::myDiv(OrdvElem* ordv, const OrdvElem* ordv1, const OrdvElem* ordv2) const
    {
#ifdef CoCoA_DEBUG
      // if debugging is active, check for underflow
      std::vector<OrdvElem> expand1(myNumIndets);
      std::vector<OrdvElem> expand2(myNumIndets);
      myDecompress(expand1, ordv1, myNumIndets);
      myDecompress(expand2, ordv2, myNumIndets);
      for (long i=0; i < myNumIndets; ++i)
      {
        CoCoA_ASSERT("Exponent Underflow" && expand1[i] >= expand2[i]);
      }
#endif
      for (long i=0; i < myOrdvWords; ++i)
        ordv[i] = ordv1[i]-ordv2[i];
    }


    inline int base::myCmp(const OrdvElem* l_ordv, const OrdvElem* r_ordv) const noexcept
    {
      for (long i=0; i < myOrdvWordsForCmp; ++i)
        if (l_ordv[i] != r_ordv[i]) return (l_ordv[i] < r_ordv[i])? -1:1;
      return 0;
    }


  } // end of namespace OrdvArith


  OrdvArith::reference NewOrdvArith(const PPOrdering& PPO);
  std::ostream& operator<<(std::ostream& out, const OrdvArith::reference& OA);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/OrdvArith.H,v 1.27 2024/04/12 16:55:52 abbott Exp $
// $Log: OrdvArith.H,v $
// Revision 1.27  2024/04/12 16:55:52  abbott
// Summary: Changed name of CPP symbol which signifies that longs are 32 bits (redmine 1661/1804)
//
// Revision 1.26  2024/03/25 08:04:25  abbott
// Summary: Modifications to work  on MinGW (& 32 bitters) (redmine 1661)
//
// Revision 1.25  2024/03/08 20:02:25  abbott
// Summary: Used constexpr (redmone 1511)
//
// Revision 1.24  2023/03/27 13:53:42  abbott
// Summary: Added new MatrixOrdering32bitImpl (to replace 32749?)
//
// Revision 1.23  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.22  2021/10/30 17:10:30  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.21  2021/02/22 21:09:25  abbott
// Summary: Implemented xel (redmine 1536)
//
// Revision 1.20  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.19  2021/01/07 14:19:36  abbott
// Summary: Corrected copyright
//
// Revision 1.18  2017/12/12 14:13:48  abbott
// Summary: Corrected use of CoCoA_DEBUG CPP symbol
//
// Revision 1.17  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.16  2015/04/13 14:42:07  abbott
// Summary: Added myPowerOverflowCheck (1st version)
// Author: JAA
//
// Revision 1.15  2014/01/28 16:53:36  abbott
// Significant change: code is now threadsafe.
// Uses CPP flag CoCoA_THREADSAFE_HACK to select between single-threaded & multi-threaded
// code; single-threaded uses a "global" buffer.
//
// Revision 1.14  2013/03/27 11:37:04  abbott
// Added new (faster, cleaner) impl for StdDegRevLex.  Still no doc though.
//
// Revision 1.13  2013/03/15 10:56:39  abbott
// Changed OrdvElem into unsigned long.
// Corrected evil/subtle bug in CmpWDegPartial.
//
// Revision 1.12  2013/02/26 21:40:43  abbott
// Fixed evil subtle bug: shift operator<< and operator>> are UNDEFINED if
// shift amount is greater than or equal to wordsize.  NASTY!!!
// Added overflow check in power function.
// Some minor cleaning.
//
// Revision 1.11  2011/03/16 15:30:09  abbott
// Removed two "unsigned" (exponent args).
//
// Revision 1.10  2011/03/16 13:22:15  abbott
// Added comments (about GrDim) for myCmpWDegPartial.
//
// Revision 1.9  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.8  2009/09/22 14:01:33  bigatti
// -- added myCmpWDegPartial (ugly name, I know....)
// -- cleaned up and realigned code in PPMonoid*.C files
//
// Revision 1.7  2008/05/30 12:44:14  abbott
// Moved "ordering matrices" into their ownn special file.
//
// Revision 1.6  2008/03/26 16:52:04  abbott
// Added exponent overflow checks (also for ordvs) when CoCoA_DEBUG is active.
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
// Revision 1.2  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.11  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.10  2007/03/07 13:49:30  bigatti
// -- minor cleanup
//
// Revision 1.9  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.8  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
//
// Revision 1.7  2006/11/24 17:29:12  cocoa
// -- reorganized includes of header files
//
// Revision 1.6  2006/11/23 17:33:10  cocoa
// -- changed: OrdvArith::base is now a class (instead of typedef)
//
// Revision 1.5  2006/11/03 15:37:47  cocoa
// -- cleaned up code after testing on usage of SmartPtrIRC
//
// Revision 1.4  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.2  2006/01/20 17:44:20  cocoa
// Friday afternoon check-in: you have been warned!
// Fixed ref counting bug in SmartPtrIRC.
// Added sum/product functions in tmp.H.
// Minor tidying here and there.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.5  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
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
// Revision 1.9  2004/11/29 16:22:35  cocoa
// -- added function for computing adjoint and inverse for DenseMatrix
//    (so adjoint/inverse matrix is computed by OrdvArith and is no
//    longer needed by PPOrdering)
//
// Revision 1.8  2004/11/12 16:11:05  cocoa
// Some missing semicolons... oops.
// Now it compiles OK.
//
// Revision 1.7  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
//
// Revision 1.6  2004/11/05 16:37:27  cocoa
// -- change: inverse matrix computed modulo 32749 using RingHom
//
// Revision 1.5  2004/11/03 17:50:31  cocoa
// -- added  mySetMatrix and mySetInverseMatrixTmp  for MatrixOrderingMod32749Impl
//
// Revision 1.4  2004/11/02 18:21:21  cocoa
// -- changed: myGetExpvBuffer --> myExpvBufferRef
//
// Revision 1.3  2004/11/02 15:05:30  cocoa
// -- new: base::myGetExpvBuffer()
// -- new: base::myComputeExpvBuffer
// -- fixed: reference count in destructor
// -- new field: myExpvBuffer
// -- changed: class base is now protected
//
// Revision 1.2  2004/10/29 15:29:55  cocoa
// -- added MatrixOrderingMod32749Impl (not tested)
// -- added myOAPtr
//
// Revision 1.1  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//

#endif
