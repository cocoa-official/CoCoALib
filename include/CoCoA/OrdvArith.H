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
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
    };


    class XelImpl: public base
    {
    public:
      XelImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
    };


    class StdDegLexImpl: public base
    {
    public:
      StdDegLexImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myStdDeg(const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
    };


    class StdDegRevLexImpl: public base
    {
    public:
      StdDegRevLexImpl(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myStdDeg(const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
    };


    class StdDegRevLexImpl2: public base
    {
    public:
      StdDegRevLexImpl2(long NumIndets);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myStdDeg(const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
    };


    class MatrixOrderingImpl: public base // !!!BUG!!! only partially implemented
    {
    public:
      MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
    private:
      int myOrderMatrixDet;
      std::vector< std::vector<int> > myOrderMatrix;
      std::vector< std::vector<int> > myAdjointOrderMatrix;
    };


    class MatrixOrderingMod32749Impl: public base
    {
    public:
      MatrixOrderingMod32749Impl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
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
      void myAssignFromExpv(OrdvElem* ordv, const std::vector<long>& expv) const override final;
      void myComputeExpv(std::vector<long>& expv, const OrdvElem* ordv) const override final;
      long myExponent(const OrdvElem* ordv, long var) const override final;
      void myOutputSelf(std::ostream& out) const override final;
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

#endif
