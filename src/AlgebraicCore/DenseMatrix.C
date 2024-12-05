//   Copyright (c)  2005-2009,2011  John Abbott and Anna M. Bigatti

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


// Source code for classes DenseMatBase and DenseMatImpl

#include "CoCoA/DenseMatrix.H"

#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/MatrixFp.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H" // for len
#include "CoCoA/PolyRing.H" // for IsPolyRing, myNumIndets in myDet(..)


//#include <vector>
using std::vector;


namespace CoCoA
{

  class DenseMatBase: public MatrixBase
  {
  protected: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
  protected:
    DenseMatBase(const ring& R, long nrows, long ncols);
  };



  class DenseMatImpl: public DenseMatBase
  {
  public:
    DenseMatImpl(const ring& R, long nrows, long ncols); // creates zero matrix
    ~DenseMatImpl();
  private: // implementation detail
    void myDtorBody();

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    DenseMatImpl* myClone() const override;
    DenseMatImpl* myZeroClone(const ring& R, long NumRows, long NumCols) const override;

    bool myIsWritable(long /*i*/, long /*j*/) const override  { return true; }
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;
// Use generic versions of the following four functions.
//     void myMulByRow(vec& lhs, const vec& v) const override;
//     void myMulByCol(vec& lhs, const vec& v) const override;
//     bool myIsZeroRow(long i) const override; ///< tests whether row i is zero
//     bool myIsZeroCol(long j) const override; ///< tests whether col j is zero
    void myRowMul(long i, ConstRefRingElem c) override; ///< row(i) = c*row(i)
    void myColMul(long j, ConstRefRingElem c) override; ///< col(j) = c*col(j)
    void myAddRowMul(long i1, long i2, ConstRefRingElem c) override; ///< row(i1) += c*row(i2)
    void myAddColMul(long j1, long j2, ConstRefRingElem c) override; ///< col(j1) += c*col(j2)
    void mySwapRows(long i1, long i2) override;
    void mySwapCols(long j1, long j2) override;
    void myResize(long NumRows, long NumCols) override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

   private: // data members
    vector< vector< RingElemRawPtr > > myEntries;
  };



  /////////////////////////////////////////////////////////////////////////////
  // Inline functions

  inline DenseMatBase::DenseMatBase(const ring& R, long NumRows, long NumCols):
      myR(R),
      myNumRowsValue(NumRows),
      myNumColsValue(NumCols)
  {}


//----------------------------------------------------------------------

  DenseMatImpl::DenseMatImpl(const ring& R, long NumRows, long NumCols):
      DenseMatBase(R, NumRows, NumCols)
  {
    constexpr long MaxNumEntries = 100*1024L*1024L;  // arbitrary upper limit
    if (NumRows > 0 && NumCols > 0 && NumRows >= MaxNumEntries/NumCols)
      CoCoA_THROW_ERROR2(ERR::ArgTooBig, "dense mat constructor");
    // If the matrix entries were full objects we could do the following one-liner:
    // myEntries.swap(vector<vector<RingElem>>(NumRows, vector<RingElem>(NumCols,zero(R))));
    // INSTEAD: we do it this way to make it exception safe (see myDtorBody)
      myEntries.reserve(myNumRows()); // not necessary but might help
      for (long i=0; i < myNumRows(); ++i)
      {
        CheckForInterrupt("DenseMatImpl ctor");
        myEntries.push_back(vector<RingElemRawPtr>());
        myEntries[i].reserve(myNumCols()); // not necessary but might help
        for (long j=0; j < myNumCols(); ++j)
          myEntries[i].push_back(myR->myNew());
      }
  }


  DenseMatImpl::~DenseMatImpl()
  {
    myDtorBody();
  }


  // Ignore the values of myNumRows and myNumCols because this procedure
  // may be called on a partially constructed DenseMat if an exception
  // is thrown during the ctor.
  void DenseMatImpl::myDtorBody()
  {
    // Written this way, we destroy in reverse order of construction.
    // NB:  myEntries.size()-1 is **unsigned** so we do not use it.
    for (long i=len(myEntries)-1; i >= 0; --i)
      for (long j=len(myEntries[i])-1; j >= 0; --j)
        myR->myDelete(myEntries[i][j]);
  }


  const ring& DenseMatImpl::myRing() const
  {
    return myR;
  }


  long DenseMatImpl::myNumRows() const
  {
    return myNumRowsValue;
  }


  long DenseMatImpl::myNumCols() const
  {
    return myNumColsValue;
  }


  inline RingElemAlias DenseMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    return RingElemAlias(myR, myEntries[i][j]);
  }


  DenseMatImpl* DenseMatImpl::myClone() const
  {
    // Making this exception safe turns out to be a bit complicated.
    DenseMatImpl* ans = nullptr;
    try
    {
      ans = new DenseMatImpl(myR, myNumRows(), myNumCols());
      for (long i=0; i < myNumRows(); ++i)
      {
        for (long j=0; j < myNumCols(); ++j)
        {
          myR->myAssign(ans->myEntries[i][j],myEntries[i][j]); // could throw
        }
      }
    }
    catch (...)
    {
      ans->myDtorBody();
      throw;
    }
    return ans;
  }


  DenseMatImpl* DenseMatImpl::myZeroClone(const ring& R, long NumRows, long NumCols) const
  {
    return new DenseMatImpl(R, NumRows, NumCols);
  }


  RingElemRawPtr DenseMatImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    return myEntries[i][j];
  }


  void DenseMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], raw(r));
  }


  void DenseMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], n);
  }


  void DenseMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], N);
  }


  void DenseMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows());
    CoCoA_ASSERT(j < myNumCols());
    myR->myAssign(myEntries[i][j], Q);
  }


//   MOVED UP TO GENERIC IMPLEMENTATIONS IN ConstMatrixBase
//   void DenseMatImpl::myMulByRow(vec& lhs, const vec& v) const
//   {
//     CoCoA_ASSERT(lhs.size() == myNumCols() && v.size() == myNumRows());
//     vector<RingElem> ans;
//     ans.resize(myNumCols(), zero(myR));
//     for (long j=0; j < myNumCols(); ++j)
//       for (long i=0; i < myNumRows(); ++i)
//         ans[j] += v[i]*myEntry(i, j);
//     // We have successfully computed the answer, now swap it in.
//     swap(lhs, ans);
// //     for (long j=0; j < myNumCols(); ++j)
// //       swap(lhs[j], ans[j]);
//   }


//   void DenseMatImpl::myMulByCol(vec& lhs, const vec& v) const
//   {
//     CoCoA_ASSERT(lhs.size() == myNumRows() && v.size() == myNumCols());
//     vector<RingElem> ans;
//     ans.resize(myNumRows(), zero(myR));
//     for (long i=0; i < myNumRows(); ++i)
//       for (long j=0; j < myNumCols(); ++j)
//         ans[i] += myEntry(i, j)*v[j];

//     // We have successfully computed the answer, now swap it in.
//     swap(lhs, ans);
// //     for (long i=0; i < myNumRows(); ++i)
// //       swap(lhs[i], ans[i]);
//   }


//   bool DenseMatImpl::myIsZeroRow(long i) const
//   {
//     if (i >= myNumRows())
//       CoCoA_THROW_ERROR(ERR::BadRowIndex, "DenseMatImpl::myIsZeroRow");
//     for (long j=0; j < myNumCols(); ++j)
//       if (!myR->myIsZero(myEntries[i][j])) return false;
//     return true;
//   }


//   bool DenseMatImpl::myIsZeroCol(long j) const
//   {
//     if (j >= myNumCols())
//       CoCoA_THROW_ERROR(ERR::BadColIndex, "DenseMatImpl::myIsZeroCol");
//     for (long i=0; i < myNumRows(); ++i)
//       if (!myR->myIsZero(myEntries[i][j])) return false;
//     return true;
//   }


  void DenseMatImpl::myRowMul(long i, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myRowMul";
    myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    if (owner(c) != myR)
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (IsOne(c)) return;
    vector<RingElem> ans;
    ans.resize(myNumCols(), zero(myR));
    for (long j = 0; j < myNumCols(); ++j)
      myR->myMul(raw(ans[j]), raw(c), myEntries[i][j]);
    // Answer successfully computed in ans, swap it into the i-th row.
    for (long j = 0; j < myNumCols(); ++j)
      myR->mySwap(myEntries[i][j], raw(ans[j]));
  }


  void DenseMatImpl::myColMul(long j, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myColMul";
    myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    if (owner(c) != myR)
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (IsOne(c)) return;
    vector<RingElem> ans;
    ans.resize(myNumRows(), zero(myR));
    for (long i = 0; i < myNumRows(); ++i)
      myR->myMul(raw(ans[i]), raw(c), myEntries[i][j]);
    // Answer successfully computed in ans, swap it into the i-th row.
    for (long i = 0; i < myNumRows(); ++i)
      myR->mySwap(myEntries[i][j], raw(ans[i]));
  }


  void DenseMatImpl::myAddRowMul(long i1, long i2, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myAddRowMul";
    myCheckRowIndex(i1, CoCoA_ERROR_CONTEXT);
    myCheckRowIndex(i2, CoCoA_ERROR_CONTEXT);
    if (owner(c) != myR)
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (IsZero(c)) return;
    const long ncols = myNumCols();
    vector<RingElem> ans(ncols, zero(myR));
    for (long j = 0; j < ncols; ++j)
    {
      myR->myAssign(raw(ans[j]),myEntries[i1][j]);
      if (!myR->myIsZero(myEntries[i2][j]))
        myR->myIsZeroAddMul(raw(ans[j]), raw(c), myEntries[i2][j]);
    }
    // Answer successfully computed in ans, swap it into the i-th row
    for (long j = 0; j < ncols; ++j)
      myR->mySwap(raw(ans[j]), myEntries[i1][j]);
  }


  void DenseMatImpl::myAddColMul(long j1, long j2, ConstRefRingElem c)
  {
    const char* const FnName = "DenseMatImpl::myAddColMul";
    myCheckColIndex(j1, CoCoA_ERROR_CONTEXT);
    myCheckColIndex(j2, CoCoA_ERROR_CONTEXT);
    if (owner(c) != myR)
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (IsZero(c)) return;
    vector<RingElem> ans;
    ans.resize(myNumRows(), zero(myR));
    for (long i = 0; i < myNumRows(); ++i)
    {
      myR->myAssign(raw(ans[i]),myEntries[i][j1]);
      myR->myIsZeroAddMul(raw(ans[i]), raw(c), myEntries[i][j2]);
    }
    // Answer successfully computed in ans, swap it into the j-th col
    for (long i = 0; i < myNumRows(); ++i)
      myR->mySwap(raw(ans[i]), myEntries[i][j1]);
  }


  void DenseMatImpl::mySwapRows(long i1, long i2)
  {
///    const char* const FnName = "DenseMatImpl::mySwapRows";
    myCheckRowIndex(i1, CoCoA_ERROR_CONTEXT);
    myCheckRowIndex(i2, CoCoA_ERROR_CONTEXT);
    if (i1 == i2) return;
    std::swap(myEntries[i1], myEntries[i2]); // trouble with iterators???
//     for (long j=0; j < myNumCols(); ++j)
//       myR->mySwap(myEntries[i1][j], myEntries[i2][j]);
  }


  void DenseMatImpl::mySwapCols(long j1, long j2)
  {
///    const char* const FnName = "DenseMatImpl::mySwapCols";
    myCheckColIndex(j1, CoCoA_ERROR_CONTEXT);
    myCheckColIndex(j2, CoCoA_ERROR_CONTEXT);
    if (j1 == j2) return;
    // Must do this the "hard way".
    for (long i=0; i < myNumRows(); ++i)
      myR->mySwap(myEntries[i][j1], myEntries[i][j2]);
  }


  // NOT FULLY EXCEPTION CLEAN (do not see how without taking a performance hit, or making it even harder to read)
  void DenseMatImpl::myResize(long NumRows, long NumCols)
  {
    // First remove any excess rows
    if (NumRows < myNumRowsValue)
    {
      for (long i=NumRows; i < myNumRowsValue; ++i)
        for (long j=0; j < myNumColsValue; ++j)
          myR->myDelete(myEntries[i][j]);
      myEntries.resize(NumRows);
    }
    const long LenEntries = len(myEntries);
    // Now resize each remaining row.
    // Two cases: lengthen or shorten  (or no change)
    if (NumCols > myNumColsValue)
    {
      for (long i=0; i < LenEntries; ++i)
        for (long j=myNumColsValue; j < NumCols; ++j)
          myEntries[i].push_back(myR->myNew()); // could throw std::bad_alloc
    }
    if (NumCols < myNumColsValue)
    {
      RingElem useless=zero(myR);
      for (long i=0; i < LenEntries; ++i)
      {
        for (long j=NumCols; j < myNumColsValue; ++j)
          myR->myDelete(myEntries[i][j]);
        myEntries[i].resize(NumCols, raw(useless)); // second arg not used
      }
    }
    // Add new zero rows if number of rows is to increase.
    if (NumRows > myNumRowsValue)
    {
      for (long i=myNumRowsValue; i < NumRows; ++i)
      {
        myEntries.push_back(vector<RingElemRawPtr>());
        for (long j=0; j < NumCols; ++j)
          myEntries[i].push_back(myR->myNew()); // could throw std::bad_alloc
      }
    }
    myNumRowsValue = NumRows;
    myNumColsValue = NumCols;
  }

  void DenseMatImpl::myDet(RingElem& d) const
  {
    const long n = myNumRows(); // we know matrix is square
    if (n <= 5) { d = DetOfSmallMat(ConstMatrixView(this)); return; }
    if (IsRingFp(myR))
    { MatrixFp M(ConstMatrixView(this)); d = det(M); return; }
    if (IsZZ(myR))
    {
      // SLUG: Criterion should depend on both dimension and entry size
      if (n > 24) d = DetByCRT(ConstMatrixView(this));
      else d = DetByBareiss(ConstMatrixView(this));
      return;
    }
    if (IsQQ(myR))
    { d = DetOverQQ(ConstMatrixView(this)); return; }
    if (IsField(myR))
    { d = DetByGauss(ConstMatrixView(this)); return; }
    if (!IsIntegralDomain(myR))
    {
      d = DetDirect(ConstMatrixView(this));
      return;
    }
    if (IsPolyRing(myR) && n < 10 && NumIndets(myR) > 2)
    {
      bool UseDetDirect = true;
      for (int i=0; i < n && UseDetDirect; ++i)
        for (int j=0; j < n && UseDetDirect; ++j)
          UseDetDirect = IsMonomial(RingElemAlias(myR,myEntries[i][j]));
      if (UseDetDirect) { d = DetDirect(ConstMatrixView(this)); return; }
    }
    /*??!!useless check!!??*/    if (IsIntegralDomain(myR))
    { d = DetByBareiss(ConstMatrixView(this)); return; }
//BUG SHOULD NEVER GET HERE --> restructure code!
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "densematrix::myDet");
  }


  long DenseMatImpl::myRank() const
  {
    vector<long> discard;
    return RankByGauss(discard, ConstMatrixView(this));
  }


  matrix NewDenseMat(const ring& R, long NumRows, long NumCols)
  {
    if (NumRows < 0)  CoCoA_THROW_ERROR2(ERR::BadIndex, "row index");
    if (NumCols < 0)  CoCoA_THROW_ERROR2(ERR::BadIndex, "column index");
    return matrix(new DenseMatImpl(R, NumRows, NumCols));
  }


  // BUG??? The next 4+4 fns should be template fns?
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<long> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }

  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigInt> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }

  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigRat> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }

  matrix NewDenseMat(const ring& R, const std::vector< std::vector<RingElem> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMat()");
    const long NumRows = len(VV);
    if (NumRows == 0) return NewDenseMat(R, 0, 0);
    const long NumCols = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[i][j]);
    return ans;
  }


  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<long> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMatTranspose()");
    const long NumCols = len(VV);
    if (NumCols == 0) return NewDenseMat(R, 0, 0);
    const long NumRows = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[j][i]);
    return ans;
  }

  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<BigInt> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMatTranspose()");
    const long NumCols = len(VV);
    if (NumCols == 0) return NewDenseMat(R, 0, 0);
    const long NumRows = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[j][i]);
    return ans;
  }

  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<BigRat> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMatTranspose()");
    const long NumCols = len(VV);
    if (NumCols == 0) return NewDenseMat(R, 0, 0);
    const long NumRows = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[j][i]);
    return ans;
  }

  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<RingElem> >& VV)
  {
    if (!IsRectangular(VV))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, "NewDenseMatTranspose()");
    const long NumCols = len(VV);
    if (NumCols == 0) return NewDenseMat(R, 0, 0);
    const long NumRows = len(VV[0]);
    matrix ans(new DenseMatImpl(R, NumRows, NumCols));
    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(ans, i, j, VV[j][i]);
    return ans;
  }


  // I bet this makes lots of wasteful copies...
  matrix NewDenseMat(const ConstMatrixView& M)
  {
    const long r = NumRows(M);
    const long c = NumCols(M);
    matrix ans(new DenseMatImpl(RingOf(M), r, c));
    for (long i=0; i < r; ++i)
      for (long j=0; j < c; ++j)
      {
        CheckForInterrupt("NewDenseMat ctor");  // wasteful to check for every matrix entry???
        SetEntry(ans, i, j, M(i, j));
      }
    return ans;
  }

} // end of namespace CoCoA
