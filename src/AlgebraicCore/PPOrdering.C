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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPOrdering.C,v 1.43 2024/08/28 15:26:39 bigatti Exp $
// $Log: PPOrdering.C,v $
// Revision 1.43  2024/08/28 15:26:39  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2
//
// Revision 1.42  2024/07/02 15:36:14  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.41  2024/03/20 15:34:11  bigatti
// Summary: simplified GradingMat (using FirstRows);
// simplified ASSERT in MatrixOrderingImpl ctor
//
// Revision 1.40  2024/03/12 14:29:32  bigatti
// Summary: NewMatrixOrdering now giving proper error is OrdMat is not full rank
//
// Revision 1.39  2023/10/05 17:54:11  abbott
// Summary: Improved IamStdGraded (redmine 1756)
//
// Revision 1.38  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.37  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.36  2021/10/30 11:53:48  abbott
// Summary: Used keyword override (redmine 1625); a little tidying too
//
// Revision 1.35  2021/02/22 21:09:25  abbott
// Summary: Implemented xel (redmine 1536)
//
// Revision 1.34  2020/06/17 19:00:33  abbott
// Summary: Updated comment
//
// Revision 1.33  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.32  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.31  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.30  2019/10/03 14:53:50  bigatti
// -- just spaces
//
// Revision 1.29  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.28  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.27  2018/05/17 15:38:49  bigatti
// -- renamed MatrixOperations --> MatrixOps
// -- sorted includes
//
// Revision 1.26  2017/12/01 17:29:21  bigatti
// // -- updated Copyright line
// // -- removed doxygen initial comment
// // -- some commented out debugging info
//
// Revision 1.25  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.24  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.23  2015/12/11 13:03:20  bigatti
// -- added IsUpperTriangular, IhaveNegEntry
// -- removed useless checks
//
// Revision 1.22  2015/12/10 14:09:45  bigatti
// -- removed IsPositiveGrading check (still working for proper fix)
//
// Revision 1.21  2015/12/08 14:05:11  abbott
// Summary: Renamed NewMatCompleteOrd to MakeTermOrd
//
// Revision 1.20  2015/12/04 15:40:32  abbott
// Summary: Removed from NewMatrixOrderig check that grading is positive
//
// Revision 1.19  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.18  2015/04/13 15:42:43  abbott
// Summary: Changed "rank" --> "rk"
// Author: JAA
//
// Revision 1.17  2014/12/10 11:39:23  bigatti
// -- added NewMatrixOrdering with 2 args (remove with 3 args?)
//
// Revision 1.16  2014/07/31 13:10:45  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.15  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.14  2014/01/28 16:44:40  abbott
// Added new fns:  IsMatrixOrdering  and  IsTermOrdering
//
// Revision 1.13  2013/07/30 14:57:45  bigatti
// -- added IamStdGraded
//
// Revision 1.12  2012/03/30 17:29:56  bigatti
// -- accepting and returning matrices over QQ
//
// Revision 1.11  2012/02/10 10:28:08  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.10  2012/02/08 16:13:41  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.9  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.8  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.7  2011/03/04 16:37:55  bigatti
// -- changed: matrix members functions args of type long instead of size_t
//
// Revision 1.6  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.5  2009/07/30 15:44:53  bigatti
// -- added comment for pseudo-constructors
//
// Revision 1.4  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.11  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.10  2007/03/07 14:06:53  bigatti
// -- CoCoA_ASSERT converted into CoCoA_ERROR (now CoCoA_THROW_ERROR)
//
// Revision 1.9  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.8  2007/01/18 17:05:22  bigatti
// -- changed namespace PPO into PPOrd (as for ModuleTermOrd)
//
// Revision 1.7  2006/11/27 14:25:15  cocoa
// -- reorganised #include files
//
// Revision 1.6  2006/11/22 14:50:33  cocoa
// -- changed: PPOrdering defined as class (instead of typedef)
//
// Revision 1.5  2006/11/16 11:27:20  cocoa
// -- reinserted myRefCountZero(): sometimes really necessary, in general safe
//
// Revision 1.4  2006/11/14 17:38:47  cocoa
// -- deleted commented out code about reference counting
// -- commented out myRefCountZero() (not necessary?)
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.5  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.2  2005/12/16 11:24:21  cocoa
// -- fixed  StdDegRevLexImpl::myOrdMatCopy  [1 --> -1]
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.2  2005/05/04 16:51:47  cocoa
// -- changed: check on matrix ordering computes "rank" (was "det")
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
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
// Revision 1.14  2004/11/29 16:22:35  cocoa
// -- added function for computing adjoint and inverse for DenseMatrix
//    (so adjoint/inverse matrix is computed by OrdvArith and is no
//    longer needed by PPOrdering)
//
// Revision 1.13  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.12  2004/11/11 14:42:01  cocoa
// -- change: cout --> GlobalLogput()
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
