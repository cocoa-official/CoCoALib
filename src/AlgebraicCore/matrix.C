//   Copyright (c)  2005,2008,2015,2020  John Abbott & Anna M. Bigatti

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


// Source code for classes matrix, MatrixView and MatrixBase, MatrixViewBase

#include "CoCoA/matrix.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/ring-AutomaticConversion.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;

namespace CoCoA
{


  // BUG: placeholder for Laura Torrente
  bool IsZero(const std::vector<RingElem>& v)
  {
    const long n = len(v);
    for (long i=0; i < n; ++i)
      if (!IsZero(v[i])) return false;
    return true;
  }


  RingElemAlias ConstMatrixView::operator()(const MachineInt& i, const MachineInt& j) const
  {  
////    const char* const FnName = "M(i,j)";
    const long I = mySmartPtr->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    const long J = mySmartPtr->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    return mySmartPtr->myEntry(I, J);
  }



  // // Check that the row index is in range (i.e. non-neg and not too large).
  // // Throws ERR::BadRowIndex if not.
  // void ConstMatrixViewBase::myCheckRowIndex(long i, const char* where) const
  // {
  //   //    if (IsNegative(i) || AsUnsignedLong(i) >= myNumRows())
  //   if (i < 0 || i >= myNumRows())
  //     CoCoA_THROW_ERROR(ERR::BadRowIndex, where);
  // }

  // long ConstMatrixViewBase::myCheckRowIndex(const MachineInt& i, const char* where) const
  // {
  //   if (IsNegative(i) || !IsSignedLong(i) || AsSignedLong(i) >= myNumRows())
  //     CoCoA_THROW_ERROR(ERR::BadRowIndex, where);
  //   return AsSignedLong(i);
  // }

    long ConstMatrixViewBase::myCheckRowIndex(const MachineInt& i, const ErrorContext& ErrCtx) const
  {
    if (IsNegative(i) || !IsSignedLong(i) || AsSignedLong(i) >= myNumRows())
//      CoCoA_THROW_ERROR_WITH_CONTEXT(ERR::BadRowIndex, ErrCtx);
      CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadIndex, "row index", ErrCtx);
    return AsSignedLong(i);
  }

  // // Check that the col index is in range (i.e. non-neg and not too large).
  // // Throws ERR::BadColIndex if not.
  // void ConstMatrixViewBase::myCheckColIndex(long j, const char* where) const
  // {
  //   //    if (IsNegative(j) || AsUnsignedLong(j) >= myNumCols())
  //   if (j < 0 || j >= myNumCols())
  //     CoCoA_THROW_ERROR(ERR::BadColIndex, where);
  // }

  // long ConstMatrixViewBase::myCheckColIndex(const MachineInt& j, const char* where) const
  // {
  //   if (IsNegative(j) || !IsSignedLong(j) || AsSignedLong(j) >= myNumCols())
  //     CoCoA_THROW_ERROR(ERR::BadColIndex, where);
  //   return AsSignedLong(j);
  // }

  long ConstMatrixViewBase::myCheckColIndex(const MachineInt& j, const ErrorContext& ErrCtx) const
  {
    if (IsNegative(j) || !IsSignedLong(j) || AsSignedLong(j) >= myNumCols())
//      CoCoA_THROW_ERROR_WITH_CONTEXT(ERR::BadColIndex, ErrCtx);
      CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadIndex, "column index", ErrCtx);
    return AsSignedLong(j);
  }


  bool ConstMatrixViewBase::IamEqual(ConstMatrixView M) const
  {
    // Naive, fully general version
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
      {
        if (myEntry(i,j) != M(i,j)) return false;
      }
    
    return true;
  }
  

  bool ConstMatrixViewBase::IamSymmetric() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      for (long j=i+1; j < myNumCols(); ++j)
        if (myEntry(i,j) != myEntry(j,i)) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IamAntiSymmetric() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myEntry(i,i))) return false;
    for (long i=0; i < myNumRows(); ++i)
      for (long j=i+1; j < myNumCols(); ++j)
        if (myEntry(i,j) != -myEntry(j,i)) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IamDiagonal() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        if (i!=j && !IsZero(myEntry(i,j))) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IamUpperTriangular() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < i; ++j)
        if (!IsZero(myEntry(i,j))) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IamLowerTriangular() const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      for (long j=i+1; j < myNumCols(); ++j)
        if (!IsZero(myEntry(i,j))) return false;
    return true;
  }
  

  bool ConstMatrixViewBase::IhaveNegEntry() const
  {
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < i; ++j)
        if (myEntry(i,j)<0) return true;
    return false;
  }
  

  bool ConstMatrixViewBase::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    for (long j=0; j < myNumCols(); ++j)
      if (!IsZero(myEntry(i,j))) return false;
    return true;
  }


  bool ConstMatrixViewBase::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myEntry(i,j))) return false;
    return true;
  }


  void ConstMatrixViewBase::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    // Put result into ans to make code exception clean.
    vector<RingElem> ans(myNumCols(), zero(myRing()));
    for (long j=0; j < myNumCols(); ++j)
      for (long i=0; i < myNumRows(); ++i)
        ans[j] += v[i]*myEntry(i, j);
    // We have successfully computed the answer, now swap it into lhs.
//???    swap(lhs, ans); // Would this invalidate iterators on lhs???
    for (long j=0; j < myNumCols(); ++j)
      swap(lhs[j], ans[j]);
  }


  void ConstMatrixViewBase::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    // Put result into ans to make code exception clean.
    vector<RingElem> ans(myNumRows(), zero(myRing()));
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        ans[i] += myEntry(i, j)*v[j];

    // We have successfully computed the answer, now swap it into lhs.
//???    swap(lhs, ans); // Would this invalidate iterators on lhs???
    for (long i=0; i < myNumRows(); ++i)
      swap(lhs[i], ans[i]);
  }


  // General case
  void ConstMatrixViewBase::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(myNumRows() == myNumCols());
    if (myNumRows() <= 5) { d = DetOfSmallMat(ConstMatrixView(this)); return; }
    if (IsPolyRing(myRing()) && NumIndets(myRing()) > 1)
    {
      d = DetDirect(ConstMatrixView(this));
      return;
    }
    if (IsIntegralDomain(myRing()))
    {
      if (IsField(myRing()))
        d = DetByGauss(ConstMatrixView(this));
      else
        d = DetByBareiss(ConstMatrixView(this));
      return;
    }
    d = DetDirect(ConstMatrixView(this));  // SLUG??? better to expand by minors?
    return;
  }


  long ConstMatrixViewBase::myRank() const
  {
    vector<long> discard;
    return RankByGauss(discard, ConstMatrixView(this));
  }


  void ConstMatrixViewBase::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    //    out << "matrix(" << myRing() << ")(" << myNumRows() << ", " << myNumCols() << ")\n[\n";
    if (IsZZ(myRing())) out << "matrix(ZZ,";
    else
    {
      if (IsQQ(myRing())) out << "matrix(QQ,";
      else
      {
        if (IsPolyRing(myRing()) && 
            myRing() == RingQQt(NumIndets(myRing())))
          out << "matrix( RingQQt(" << NumIndets(myRing()) << "),";
        else out << "matrix( /*" << myRing() << "*/";
      }
    }
    out << "\n [[";
    for (long i=0; i < myNumRows(); ++i)
    {
      if (i >0)  out << "],\n  [";
      for (long j=0; j < myNumCols(); ++j)
      {
        if (j > 0) out << ", ";
        out << myEntry(i, j);
      }
    }
    out << "]])";
  }


  void ConstMatrixViewBase::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("linalg2", "matrix");
    OMOut << myRing() << myNumRows() << myNumCols();
    for (long i=0; i < myNumRows(); ++i)
    {
      OMOut->mySendApplyStart();
      OMOut << OpenMathSymbol("linalg2", "matrixrow"); //??? << myNumCols; ???
      for (long j=0; j < myNumCols(); ++j)
      {
        OMOut << myEntry(i, j); //??? print only raw value???
      }
      OMOut->mySendApplyEnd();
    }
    OMOut->mySendApplyEnd();
  }


  /***************************************************************************/

  ConstMatrix::ConstMatrix(const ConstMatrix& M):
      ConstMatrixView(M->myClone())
  {}


  matrix::matrix(const matrix& M):
      MatrixView(M->myClone())
  {}

  /***************************************************************************/

  std::ostream& operator<<(std::ostream& out, const ConstMatrixView& M)
  {
    if (!out) return out;  // short-cut for bad ostreams
    M->myOutputSelf(out);
    return out;
  }


  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, ConstRefRingElem r)
  {
    const char* const FnName = "SetEntry(M,i,j,r)";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,FnName);
    const long I = M->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    const long J = M->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    if (!M->myIsWritable(I, J))
      CoCoA_THROW_ERROR(ERR::BadMatrixSetEntry, FnName);
    if (owner(r) == RingOf(M)) { M->mySetEntry(I, J, r); return; }
    const RingHom promote = AutomaticConversionHom(owner(r),RingOf(M),ErrMixed); // throws ErrMixed if not permitted 
    if (codomain(promote) == owner(r))
      CoCoA_THROW_ERROR("Cannot automatically map ring elem into smaller ring", FnName);
    M->mySetEntry(I,J, promote(r));
  }


  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, const MachineInt& n)
  {
    const char* const FnName = "SetEntry(M,i,j,n)";
    const long I = M->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    const long J = M->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    if (!M->myIsWritable(I, J))
      CoCoA_THROW_ERROR(ERR::BadMatrixSetEntry, FnName);
    M->mySetEntry(I, J, n);
  }


  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, const BigInt& N)
  {
    const char* const FnName = "SetEntry(M,i,j,N)";
    const long I = M->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    const long J = M->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    if (!M->myIsWritable(I, J))
      CoCoA_THROW_ERROR(ERR::BadMatrixSetEntry, FnName);
    M->mySetEntry(I, J, N);
  }


  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, const BigRat& Q)
  {
    const char* const FnName = "SetEntry(M,i,j,Q)";
    const long I = M->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    const long J = M->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    if (!M->myIsWritable(I, J))
      CoCoA_THROW_ERROR(ERR::BadMatrixSetEntry, FnName);
    M->mySetEntry(I, J, Q);
  }


  void SwapRows(matrix& M, const MachineInt& i1, const MachineInt& i2)
  {
////    const char* const FnName = "SwapRows";
    const long I1 = M->myCheckRowIndex(i1, CoCoA_ERROR_CONTEXT);
    const long I2 = M->myCheckRowIndex(i2, CoCoA_ERROR_CONTEXT);
    if (I1 == I2) return;
    return M->mySwapRows(I1,I2);
  }


  void SwapCols(matrix& M, const MachineInt& j1, const MachineInt& j2)
  {
////    const char* const FnName = "SwapCols";
    const long J1 = M->myCheckColIndex(j1, CoCoA_ERROR_CONTEXT);
    const long J2 = M->myCheckColIndex(j2, CoCoA_ERROR_CONTEXT);
    if (J1 == J2) return;
    return M->mySwapCols(J1, J2);
  }


  void swap(matrix& M1, matrix& M2)
  {
    M1.mySmartPtr.mySwap(M2.mySmartPtr);
  }


  void DeleteRow(matrix& M, const MachineInt& i)
  {
    const long I = M->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    const long nrows = NumRows(M);
    for (long row=I+1; row<nrows; ++row)
      M->mySwapRows(row, row-1);
    M->myResize(NumRows(M)-1, NumCols(M));
  }


  void DeleteCol(matrix& M, const MachineInt& j)
  {
    const long J = M->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    for (long col=J+1; col<NumCols(M); ++col)
      M->mySwapCols(col,col-1);
    M->myResize(NumRows(M), NumCols(M)-1);
  }


  void AddRowMul(matrix& M, long i, long j, ConstRefRingElem c)
  { M->myAddRowMul(i, j, c); }

  void AddColMul(matrix& M, long i, long j, ConstRefRingElem c)
  { M->myAddColMul(i, j, c); }


  /***************************************************************************/
  // ZeroMatImpl -- zero matrices
  
  class ZeroMatImpl: public ConstMatrixBase
  {
  private:
    friend ConstMatrix ZeroMat(const ring&, const MachineInt&, const MachineInt&); // pseudo-ctor
    ZeroMatImpl(const ring& R, long nrows, long ncols);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ZeroMatImpl(const ZeroMatImpl&) = delete;
    ZeroMatImpl& operator=(const ZeroMatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
    bool IamEqual(ConstMatrixView M) const override;
    bool IamSymmetric() const override  {return true;}
    bool IamAntiSymmetric() const override  {return true;}
    bool IamDiagonal() const override  {return true;}
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    virtual ConstMatrixBase* myClone() const override;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
  };


  /******************************************************************/
  // Functions for ZeroMatImpl

  bool IsZeroMatImpl(const ConstMatrixView& M)  
  {return dynamic_cast<ZeroMatImpl*>(M.mySmartPtr.myRawPtr()) != nullptr;}


  ZeroMatImpl::ZeroMatImpl(const ring& R, long NumRows, long NumCols):
      ConstMatrixBase(),
      myR(R),
      myNumRowsValue(NumRows),
      myNumColsValue(NumCols)
  {}


  const ring& ZeroMatImpl::myRing() const
  {
    return myR;
  }


  long ZeroMatImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long ZeroMatImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  RingElemAlias ZeroMatImpl::myEntry(long i, long j) const
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return zero(myR);
  }


  // This is exception safe only if "lhs[i]=0;" cannot throw.
  void ZeroMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    (void)(v); // to avoid compiler warning about unused parameter
// BUG: should check all elems of v belong to the right ring!
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    for (long i=0; i < myNumCols(); ++i)
      lhs[i] = 0;
  }


  // This is exception safe only if "lhs[i]=0;" cannot throw.
  void ZeroMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    (void)(v); // to avoid compiler warning about unused parameter
// BUG: should check all elems of v belong to the right ring!
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      lhs[i] = 0;
  }


  bool ZeroMatImpl::IamEqual(ConstMatrixView M) const
  {
    //    if (myRing() != RingOf(M)) return false;
    if (myRing() != RingOf(M))
      CoCoA_THROW_ERROR(ERR::MixedRings, "IamEqual");
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (myNumCols() == 0 || myNumRows() == 0) return true;
    if (IsZeroMatImpl(M))     return true;
    if (IsIdentityMatImpl(M)) return false;
    if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "ZeroMatImpl::IamEqual - diag" << std::endl;
      for (long i=0; i < myNumRows(); ++i) if (!IsZero(M(i,i))) return false;
      return true;
    }
    //    std::cout << "ZeroMatImpl::IamEqual" << std::endl;
    for (long i=0; i < myNumRows(); ++i)
      for (long j=0; j < myNumCols(); ++j)
        if (!IsZero(M(i,j))) return false;
    //    return ConstMatrixViewBase::IamEqual(M);
    return true;
  }


  bool ZeroMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= 0 && i < myNumRows());
    return true;
  }


  bool ZeroMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return true;
  }


  void ZeroMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    if (myNumRows() == 0)
      d = 1;
    else
      d = 0;
  }



  long ZeroMatImpl::myRank() const
  {
    return 0;
  }


  ConstMatrixBase* ZeroMatImpl::myClone() const
  {
    return new ZeroMatImpl(myR, myNumRowsValue, myNumColsValue);
  }


  ConstMatrix ZeroMat(const ring& R, const MachineInt& NumRows, const MachineInt& NumCols)
  {
    if (IsNegative(NumRows) || !IsSignedLong(NumRows)) CoCoA_THROW_ERROR2(ERR::BadIndex, "num rows");
    if (IsNegative(NumCols) || !IsSignedLong(NumCols)) CoCoA_THROW_ERROR2(ERR::BadIndex, "num cols");
    return ConstMatrix(new ZeroMatImpl(R, AsSignedLong(NumRows), AsSignedLong(NumCols)));
  }


  /***************************************************************************/
  // IdentityMatImpl -- identity matrices (necessarily square)

  class IdentityMatImpl: public ConstMatrixBase
  {
  private:
    friend ConstMatrix IdentityMat(const ring& R, const MachineInt& dim); // pseudo-ctor
    IdentityMatImpl(const ring& R, long dim);
    // default dtor is fine
  public: // disable default copy ctor and assignment
    IdentityMatImpl(const IdentityMatImpl&) = delete;
    IdentityMatImpl& operator=(const IdentityMatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
    bool IamEqual(ConstMatrixView M) const override;
    bool IamSymmetric() const override  {return true;}
    bool IamAntiSymmetric() const override  {return myNumRows()==0;}
    bool IamDiagonal() const override  {return true;}
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    ConstMatrixBase* myClone() const override;

  private: // data members
    const ring myR;
    long myDim;
  };


  /***************************************************************************/
  // Functions for IdentityMatImpl

  bool IsIdentityMatImpl(const ConstMatrixView& M)
  {return dynamic_cast<const IdentityMatImpl*>(M.mySmartPtr.myRawPtr()) != nullptr;}


  IdentityMatImpl::IdentityMatImpl(const ring& R, long dim):
      ConstMatrixBase(),
      myR(R),
      myDim(dim)
  {}


  const ring& IdentityMatImpl::myRing() const
  {
    return myR;
  }


  long IdentityMatImpl::myNumRows() const
  {
    return myDim;
  }


  long IdentityMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias IdentityMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == j) return one(myR);
    return zero(myR);
  }


  //BUG: should use FreeModule elems!
  void IdentityMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
    lhs = v;
  }


  void IdentityMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    lhs = v;
  }


  bool IdentityMatImpl::IamEqual(ConstMatrixView M) const
  {
    if (myRing() != RingOf(M)) return false;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (IsZeroMatImpl(M)) return NumCols(M) == 0;
    if (IsIdentityMatImpl(M)) return true;
    if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "IdentityMatImpl::IamEqual - diag" << std::endl;
      for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
      return true;
    }
    return ConstMatrixViewBase::IamEqual(M);
  }


  bool IdentityMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= i && i < myNumRows());
    return false;
  }


  bool IdentityMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return false;
  }


  void IdentityMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    d = 1; // NB correct even for 0x0 matrix.
  }


  long IdentityMatImpl::myRank() const
  {
    return myNumRows();
  }


  ConstMatrixBase* IdentityMatImpl::myClone() const
  {
    return new IdentityMatImpl(myR, myDim);
  }


  ConstMatrix IdentityMat(const ring& R, const MachineInt& dim)
  {
    if (IsNegative(dim) || !IsSignedLong(dim))  CoCoA_THROW_ERROR2(ERR::BadIndex, "dim");
    return ConstMatrix(new IdentityMatImpl(R, AsSignedLong(dim)));
  }


} // end of namespace CoCoA
