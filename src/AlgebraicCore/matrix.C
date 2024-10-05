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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/matrix.C,v 1.61 2024/07/22 15:43:58 abbott Exp $
// $Log: matrix.C,v $
// Revision 1.61  2024/07/22 15:43:58  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.60  2024/07/15 17:12:42  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.59  2024/07/14 10:29:08  abbott
// Summary: First attempt to use new throw-error macro
//
// Revision 1.58  2024/03/18 08:48:18  abbott
// Summary: Removed incorrect assertion (redmine 1085)
//
// Revision 1.57  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.56  2024/03/15 18:57:33  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.55  2022/08/05 20:48:33  abbott
// Summary: Took NumRows/NumCols out of for-loop exit conds
//
// Revision 1.54  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.53  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.52  2021/10/30 19:36:35  abbott
// Summary: Added in some missing "override" (according to clang)
//
// Revision 1.51  2021/10/30 16:48:28  abbott
// Summary: Used keywords override and delete (redmine 1625 & 1627)
//
// Revision 1.50  2020/10/05 19:27:32  abbott
// Summary: Changed SetEntry to do auto ringhom mapping
//
// Revision 1.49  2020/09/28 11:19:56  abbott
// Summary: Improved dispatch fn myDet
//
// Revision 1.48  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.47  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.46  2020/02/11 16:13:17  abbott
// Summary: Added check for bad ostream (see redmine 969)
//
// Revision 1.45  2020/01/27 19:55:46  abbott
// Summary: Added check for trivial case in SwapRows and SwapCols
//
// Revision 1.44  2019/10/15 11:54:09  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.43  2019/09/11 09:53:38  abbott
// Summary: Added call to DetOfSmallMat if size <= 5
//
// Revision 1.42  2018/05/17 15:57:00  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.41  2018/02/15 14:38:58  abbott
// Summary: For anna
//
// Revision 1.40  2016/02/09 15:04:11  bigatti
// -- added AddRowMul, AddColMul
//
// Revision 1.39  2015/12/11 15:51:23  bigatti
// -- added IsLowerTriangular
//
// Revision 1.38  2015/12/11 13:03:20  bigatti
// -- added IsUpperTriangular, IhaveNegEntry
// -- removed useless checks
//
// Revision 1.37  2015/12/01 12:28:31  abbott
// Summary: Change arg type from long to MachineInt: see issue 830
//
// Revision 1.36  2015/12/01 10:01:31  abbott
// Summary: Removed AssignZero (& myAssignZero) for matrices
//
// Revision 1.35  2015/11/30 21:53:56  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.34  2015/04/13 16:28:25  abbott
// Summary: Changed ConstMatrixViewBase::myDet: now chooses betw DetByGauss & DetDirect
// Author: JAA
//
// Revision 1.33  2014/08/16 13:42:26  abbott
// Summary: Added copy ctor for matrix, and myClone for MatrixBase
// Author: JAA
//
// Revision 1.32  2014/07/30 14:13:05  abbott
// Summary: Changed BaseRing into RingOf; myBaseRing --> myRing
// Author: JAA
//
// Revision 1.31  2014/07/07 13:27:29  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.30  2014/04/30 16:27:30  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.29  2014/04/17 14:08:05  bigatti
// -- moved Is... functions from matrix.H to here
//
// Revision 1.28  2014/04/11 15:44:28  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.27  2014/03/19 15:58:16  bigatti
// -- now printing the ring explicitely (when possible)
//
// Revision 1.26  2013/07/12 15:07:18  abbott
// Corrected a string (FnName passed to index checking fn).
//
// Revision 1.25  2013/06/17 12:56:57  bigatti
// -- added DeleteRow, DeleteCol
//
// Revision 1.24  2013/03/28 15:15:34  abbott
// Modified printing of matrices (removed two newlines); consequential changes to expected output of certain tests.
//
// Revision 1.23  2012/10/24 12:24:42  abbott
// Changed return type of matrix::operator().
//
// Revision 1.22  2012/10/16 09:50:40  abbott
// Changed  RefRingElem  into  RingElem&.
// Realigned some comments.
//
// Revision 1.21  2012/04/04 13:53:54  bigatti
// -- added SwapRows, SwapCols
//
// Revision 1.20  2012/03/27 11:53:29  bigatti
// -- changed printing style (and added ring)
//
// Revision 1.19  2011/11/09 17:08:53  bigatti
// -- changed output for matrix
//
// Revision 1.18  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.17  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.15  2011/03/10 17:57:57  bigatti
// -- fixed some CoCoA_ASSERT
//
// Revision 1.14  2011/03/09 09:11:20  bigatti
// -- removed call to AsSignedLong
//
// Revision 1.13  2011/03/08 17:27:01  bigatti
// -- changed: args for rows and cols are now  long  instead of  MachineInt
//
// Revision 1.12  2011/03/04 16:26:14  bigatti
// -- changed: functions args of type MachineInt instead of size_t
//             members functions args of type long instead of size_t
//
// Revision 1.11  2011/03/03 13:50:21  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.10  2011/02/22 14:23:59  bigatti
// -- fixed IamAntiSymmetric
//
// Revision 1.9  2011/02/22 13:15:10  bigatti
// -- added IsSymmetric, IsAntiSymmetric, IsDiagonal, and tests
//
// Revision 1.8  2011/02/15 10:04:05  bigatti
// -- added myIsEqual (default impl)
// -- added check for myIsWritable in SetEntry (4 defn)
//
// Revision 1.7  2011/01/31 14:10:30  bigatti
// -- added mySetEntry/SetEntry with BigRat entry
//
// Revision 1.6  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.5  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.4  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.3  2006/11/02 13:25:43  cocoa
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
// Revision 1.3  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.2  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.11  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.10  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.9  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.8  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.7  2005/03/31 16:59:42  cocoa
// Made special matrix ctors private, so a user has to pass via the
// pseudo-ctors (which do all the arg sanity checking).
//
// Revision 1.6  2005/03/30 17:15:14  cocoa
// Cleaned the SpecialMatrix code; a simple test compiles and
// runs fine.
//
// Revision 1.5  2005/03/29 17:19:59  cocoa
// -- added: rank
//
// Revision 1.4  2005/03/11 16:44:18  cocoa
// New abstract class structure for matrices.
// New types of special matrix.
//
// Revision 1.3  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.9  2004/11/29 16:19:02  cocoa
// -- changed syntax for function det
//
// Revision 1.8  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.7  2004/11/11 13:55:10  cocoa
// -- minor changes for doxygen
//
// Revision 1.6  2004/11/02 14:49:03  cocoa
// -- new code for matrix orderings
//
// Revision 1.5  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.4  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.3  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.2  2004/03/09 16:15:41  cocoa
// First version of matrices.  A simple exmaple compiles and runs,
// so I'm checking in.
//
// Revision 1.1  2004/01/28 16:29:17  cocoa
// First attempt at implementing matrices.
//
//
