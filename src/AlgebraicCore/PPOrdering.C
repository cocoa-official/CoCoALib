//   Copyright (c)  2004-2017,2021  John Abbott and Anna M. Bigatti
//   Author:  2004-2011,2015  John Abbott

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

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"  // for submat
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"  // for LongRange

#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;

namespace CoCoA
{

  std::ostream& operator<<(std::ostream& out, const PPOrdering& PPO)
  {
    if (!out) return out;  // short-cut for bad ostreams
    PPO->myOutputSelf(out);
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const PPOrdering& PPO)
  {
    PPO->myOutputSelf_OM(OMOut);
    return OMOut;
  }


  //---------------------------------------------------------------------------

  PPOrderingBase::PPOrderingBase(long NumIndets, long GradingDim):
    IntrusiveReferenceCount(),
    myNumIndets(NumIndets),
    myGradingDim(GradingDim)
  {
    CoCoA_ASSERT(NumIndets > 0);
    CoCoA_ASSERT(0 <= GradingDim && GradingDim <= NumIndets);
  }


  PPOrderingBase::~PPOrderingBase()
  {}



  //---------------------------------------------------------------------------

  namespace PPOrd
  {

    // created only via pseudo-ctor LexCtor::op()
    class LexImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering LexCtor::operator()(const MachineInt&) const;
      LexImpl(long NumIndets);
      LexImpl(const LexImpl&) = delete;
      LexImpl& operator=(const LexImpl&) = delete;
    public: // fns every PPOrdering must implement
      void myOutputSelf(std::ostream& out) const override;
      void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
      const ConstMatrixView& myOrdMat() const override;
      bool IamStdGraded() const override  { return (myNumIndets == 1); }
    private: // data member
      ConstMatrix myM;
    };


    // created only via pseudo-ctor LexCtor::op()
    class XelImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering XelCtor::operator()(const MachineInt&) const;
      XelImpl(long NumIndets);
      XelImpl(const LexImpl&) = delete;
      XelImpl& operator=(const LexImpl&) = delete;
    public: // fns every PPOrdering must implement
      void myOutputSelf(std::ostream& out) const override;
      void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
      const ConstMatrixView& myOrdMat() const override;
      bool IamStdGraded() const override  { return (myNumIndets == 1); }
    private: // data member
      ConstMatrix myM;
    };


    // created only via pseudo-ctor StdDegLexCtor::op()
    class StdDegLexImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering StdDegLexCtor::operator()(const MachineInt&) const;
      StdDegLexImpl(long NumIndets);
      StdDegLexImpl(const StdDegLexImpl&) = delete;
      StdDegLexImpl& operator=(const StdDegLexImpl&) = delete;
    public: // fns every PPOrdering must implement
      void myOutputSelf(std::ostream& out) const override;
      void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
      const ConstMatrixView& myOrdMat() const override;
      bool IamStdGraded() const override  { return true; }
    private: // data member
      ConstMatrix myM;
    };


    // created only via pseudo-ctor StdDegRevLexCtor::op()
    class StdDegRevLexImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering StdDegRevLexCtor::operator()(const MachineInt&) const;
      StdDegRevLexImpl(long NumIndets);
      StdDegRevLexImpl(const StdDegRevLexImpl&) = delete;
      StdDegRevLexImpl& operator=(const StdDegRevLexImpl&) = delete;
    public: // fns every PPOrdering must implement
      void myOutputSelf(std::ostream& out) const override;
      void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
      const ConstMatrixView& myOrdMat() const override;
      bool IamStdGraded() const override  { return true; }
    private: // data member
      ConstMatrix myM;
    };


    class MatrixOrderingImpl: public PPOrderingBase
    {
    private:
      friend PPOrdering CoCoA::NewMatrixOrdering(const ConstMatrixView& OrderMatrix, const MachineInt& GradingDim);
      MatrixOrderingImpl(/*long NumIndets,*/ const ConstMatrixView& OrderMatrix, long GradingDim);
      MatrixOrderingImpl(const MatrixOrderingImpl&) = delete;
      MatrixOrderingImpl& operator=(const MatrixOrderingImpl&) = delete;
    public: // fns every PPOrdering must implement
      void myOutputSelf(std::ostream& out) const override;
      void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
      const ConstMatrixView& myOrdMat() const override;
      bool IamStdGraded() const override;
    private: ///< data members (in addition to those inherited)
      matrix myDefiningMatrix;
    };


  } // end of namespace PPO


  //----------------------------------------------------------------------
  // Here is the pseudo-ctor for a matrix ordering:

  PPOrdering NewMatrixOrdering(const ConstMatrixView& OrderMatrix, const MachineInt& GradingDim)
  {
    if (IsNegative(GradingDim) || !IsSignedLong(GradingDim) || !IsInRange(0, GradingDim, NumCols(OrderMatrix)))
      CoCoA_THROW_ERROR2(ERR::OutOfRange, "GradingDim");
    //    std::cout << "--NewMatrixOrdering-called" << std::endl;
    const long GrDim = AsSignedLong(GradingDim);
    if (!IsTermOrdering(OrderMatrix))
    {
      if (rk(OrderMatrix) != NumRows(OrderMatrix))
        CoCoA_THROW_ERROR1(ERR::ReqFullRank);
      CoCoA_THROW_ERROR1(ERR::ReqTermOrdering);
    }
    return PPOrdering(new PPOrd::MatrixOrderingImpl(OrderMatrix, GrDim));
  }


  //------------------------------------------------------------------
  // Implementations

  namespace PPOrd
  {

    LexImpl::LexImpl(long NumIndets):
        PPOrderingBase(NumIndets, 0),
        myM(IdentityMat(RingZZ(),NumIndets))
    {}


    void LexImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "PPOrderingLex(" << myNumIndets << ")";
    }


    void LexImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingLex");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    const ConstMatrixView& LexImpl::myOrdMat() const
    {
      return myM;
    }


    //-------------------------------------------------------

    XelImpl::XelImpl(long NumIndets):
        PPOrderingBase(NumIndets, 0),
        myM(XelMat(NumIndets))
    {}


    void XelImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "PPOrderingXel(" << myNumIndets << ")";
    }


    void XelImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingXel");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    const ConstMatrixView& XelImpl::myOrdMat() const
    {
      return myM;
    }


    //-------------------------------------------------------

    StdDegLexImpl::StdDegLexImpl(long NumIndets):
        PPOrderingBase(NumIndets, 1),
        myM(StdDegLexMat(NumIndets))
    {}


    void StdDegLexImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "PPOrderingStdDegLex(" << myNumIndets << ")";
    }


    void StdDegLexImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingStdDegLex");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    const ConstMatrixView& StdDegLexImpl::myOrdMat() const
    {
      return myM;
    }


    //------------------------------------------------------------//


    StdDegRevLexImpl::StdDegRevLexImpl(long NumIndets):
        PPOrderingBase(NumIndets, 1),
        myM(StdDegRevLexMat(NumIndets))
    {}


    void StdDegRevLexImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "PPOrderingStdDegRevLex(" << myNumIndets << ")";
    }


    void StdDegRevLexImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingStdDegRevLex");
      OMOut << myNumIndets;
      OMOut->mySendApplyEnd();
    }


    const ConstMatrixView& StdDegRevLexImpl::myOrdMat() const
    {
      return myM;
    }


    //------------------------------------------------------------//


    MatrixOrderingImpl::MatrixOrderingImpl(/*long NumIndets,*/ const ConstMatrixView& OrderMatrix, long GradingDim):
        PPOrderingBase(NumCols(OrderMatrix), AsSignedLong(GradingDim)),
      myDefiningMatrix(NewDenseMat(RingZZ(), myNumIndets, myNumIndets))  // BUG!!! Nasty if NumIndets is large!!!
    {
      // This ctor ***ASSUMES*** the caller knows that OrderMatrix and GradingDim are good values.
      // Note: ctor for PPOrderingBase checks value of GradingDim and myNumIndets
      // AMB 2024-03: I simplified but not tested the following ASSERT call
      CoCoA_ASSERT(!HasNegEntry(FirstRows(OrderMatrix, myGradingDim)));
      CoCoA_ASSERT(IsTermOrdering(OrderMatrix));
      for (int i=0; i < myNumIndets; ++i)
        for (int j=0; j < myNumIndets; ++j)
        {
          SetEntry(myDefiningMatrix, i,j, ConvertTo<BigInt>(OrderMatrix(i,j)));
        }
    }


    void MatrixOrderingImpl::myOutputSelf(std::ostream& out) const
    {
      if (!out) return;  // short-cut for bad ostreams
      out << "PPOrderingMatrix(GrDim=" << myGradingDim << ", " << myDefiningMatrix << ")";
    }


    void MatrixOrderingImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("cocoa", "PPOrderingMatrix");
      OMOut << myNumIndets;
      OMOut << myGradingDim;
      OMOut->mySendApplyEnd();
    }


    const ConstMatrixView& MatrixOrderingImpl::myOrdMat() const
    {
      return myDefiningMatrix;
    }


    bool MatrixOrderingImpl::IamStdGraded() const
    {
      // ???Cache the result in a data member???
      if (myGradingDim != 1)  return false;
      if (IsZero(myDefiningMatrix(0,0)))  return false;
      for (long j=1; j < myNumIndets; ++j)
        if (myDefiningMatrix(0,j) != myDefiningMatrix(0,0))  return false;
      return true;
    }
    

    //------------------------------------------------------------//


  } // end of namespace PPO


  bool IsLex(const PPOrdering& PPO)
  {
    if (dynamic_cast<const PPOrd::LexImpl*>(PPO.myRawPtr())) return true;
    if (!dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())) return false;
    // must decide whether the matrix is Lex, possibly in disguise
    return false; // just to keep compiler quiet for the moment
  }


  bool IsStdDegLex(const PPOrdering& PPO)
  {
    if (dynamic_cast<const PPOrd::StdDegLexImpl*>(PPO.myRawPtr())) return true;
    if (!dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())) return false;
    // must decide whether the matrix is StdDegLex, possibly in disguise
    return false; // just to keep compiler quiet for the moment
  }


  bool IsStdDegRevLex(const PPOrdering& PPO)
  {
    if (dynamic_cast<const PPOrd::StdDegRevLexImpl*>(PPO.myRawPtr())) return true;
    if (!dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr())) return false;
    // must decide whether the matrix is StdDegRevLex, possibly in disguise
    return false; // just to keep compiler quiet for the moment
  }


  bool IsMatrixOrdering(const PPOrdering& PPO)
  {
    return (dynamic_cast<const PPOrd::MatrixOrderingImpl*>(PPO.myRawPtr()) != nullptr);
  }


  bool IsTermOrdering(const PPOrdering& PPO)
  {
    return IsLex(PPO) ||
           IsStdDegLex(PPO) ||
           IsStdDegRevLex(PPO) ||
           IsTermOrdering(OrdMat(PPO));
  }


  ConstMatrixView OrdMat(const PPOrdering& PPO)
  { return PPO->myOrdMat(); }


  ConstMatrixView GradingMat(const PPOrdering& PPO)
  { return FirstRows(OrdMat(PPO), GradingDim(PPO)); }


  //------------------------------------------------------------------
  // Common ordering pseudo-ctors
  
  long CheckIsPositive(const MachineInt& n, const char* const FnName)
  {
    if (IsNegative(n) || IsZero(n) || !IsSignedLong(n))
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    return AsSignedLong(n);
  }


  PPOrdering LexCtor::operator()(const MachineInt& NumIndets) const
  { return PPOrdering(new PPOrd::LexImpl(CheckIsPositive(NumIndets, "lex ordering"))); }


  PPOrdering XelCtor::operator()(const MachineInt& NumIndets) const
  { return PPOrdering(new PPOrd::XelImpl(CheckIsPositive(NumIndets, "xel ordering"))); }


  PPOrdering StdDegLexCtor::operator()(const MachineInt& NumIndets) const
  { return PPOrdering(new PPOrd::StdDegLexImpl(CheckIsPositive(NumIndets, "DegLex ordering"))); }


  PPOrdering StdDegRevLexCtor::operator()(const MachineInt& NumIndets) const
  { return PPOrdering(new PPOrd::StdDegRevLexImpl(CheckIsPositive(NumIndets, "DegRevLex ordering"))); }



  LexCtor lex;
  XelCtor xel;
  StdDegLexCtor StdDegLex;
  StdDegRevLexCtor StdDegRevLex;


} // end of namespace CoCoA
