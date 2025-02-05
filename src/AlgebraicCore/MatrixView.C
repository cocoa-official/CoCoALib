//   Copyright (c)  2005,2008,2021  John Abbott and Anna M. Bigatti

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


#include "CoCoA/MatrixView.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H" // for len
#include "CoCoA/VectorOps.H" // for HasUniqueOwner


//#include <vector>
using std::vector;
#include <algorithm>
using std::min;
//#include <limits> // included in utils.H
using std::numeric_limits;
//#include <iostream> // for debugging only

namespace CoCoA
{

  /***************************************************************************/
  bool IsConstDiagMatImpl(ConstMatrixView M);
  bool IsDiagMatImpl(ConstMatrixView M);

  /***************************************************************************/
  // ConstTransposeMatImpl -- a transposed view of a const matrix

  class ConstTransposeMatImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView transpose(ConstMatrixView M); // pseudo-ctor
    ConstTransposeMatImpl(ConstMatrixView M);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ConstTransposeMatImpl(const ConstTransposeMatImpl&) = delete;
    ConstTransposeMatImpl& operator=(const ConstTransposeMatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

  private: // data members
    ConstMatrixView myM;
  };


  /***************************************************************************/
  // Functions for ConstTransposeMatImpl -- all are simple one-liners

  ConstTransposeMatImpl::ConstTransposeMatImpl(ConstMatrixView M):
      myM(M)
  {}


  const ring& ConstTransposeMatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long ConstTransposeMatImpl::myNumRows() const
  {
    return NumCols(myM);
  }


  long ConstTransposeMatImpl::myNumCols() const
  {
    return NumRows(myM);
  }


  RingElemAlias ConstTransposeMatImpl::myEntry(long i, long j) const
  {
    return myM(j, i);
  }


  void ConstTransposeMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    myM->myMulByCol(lhs, v);
  }


  void ConstTransposeMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    myM->myMulByRow(lhs, v);
  }



  bool ConstTransposeMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return myM->myIsZeroCol(i);
  }


  bool ConstTransposeMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myM->myIsZeroRow(j);
  }


  void ConstTransposeMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    myM->myDet(d);
  }


  long ConstTransposeMatImpl::myRank() const
  {
    return myM->myRank();
  }


  ConstMatrixView transpose(ConstMatrixView M)
  {
    return ConstMatrixView(new ConstTransposeMatImpl(M));
  }


  /***************************************************************************/
  // TransposeMatImpl -- a transposed view of a (non-const) matrix

  class TransposeMatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView transpose(MatrixView M); // pseudo-ctor
    TransposeMatImpl(MatrixView M);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    TransposeMatImpl(const TransposeMatImpl&) = delete;
    TransposeMatImpl& operator=(const TransposeMatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    // non-const member fns
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    MatrixView myM;
  };


  /***************************************************************************/
  // Functions for TransposeMatImpl -- all are simple one-liners

  TransposeMatImpl::TransposeMatImpl(MatrixView M):
      myM(M)
  {}


  const ring& TransposeMatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long TransposeMatImpl::myNumRows() const
  {
    return NumCols(myM);
  }


  long TransposeMatImpl::myNumCols() const
  {
    return NumRows(myM);
  }


  RingElemAlias TransposeMatImpl::myEntry(long i, long j) const
  {
    return myM(j, i);
  }


  void TransposeMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    myM->myMulByCol(lhs, v);
  }


  void TransposeMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    myM->myMulByRow(lhs, v);
  }



  bool TransposeMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return myM->myIsZeroCol(i);
  }


  bool TransposeMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myM->myIsZeroRow(j);
  }


  void TransposeMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    CoCoA_ASSERT(myNumRows() == myNumCols());
    myM->myDet(d);
  }


  long TransposeMatImpl::myRank() const
  {
    return myM->myRank();
  }


  bool TransposeMatImpl::myIsWritable(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM->myIsWritable(j,i);
  }


  RingElemRawPtr TransposeMatImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    return myM->myRawEntry(j,i);
  }


  void TransposeMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,r);
  }


  void TransposeMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,n);
  }


  void TransposeMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,N);
  }


  void TransposeMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myM->mySetEntry(j,i,Q);
  }


  MatrixView transpose(MatrixView M)
  {
    return MatrixView(new TransposeMatImpl(M));
  }


  /***************************************************************************/
  // SubmatImpl -- see only certain rows/cols of a matrix

  class ConstSubmatImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView submat(ConstMatrixView M, const std::vector<long>& RowIndices, const std::vector<long>& ColIndices); // pseudo-ctor
    ConstSubmatImpl(ConstMatrixView M, const std::vector<long>& rows, const std::vector<long>& cols);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ConstSubmatImpl(const ConstSubmatImpl&) = delete;
    ConstSubmatImpl& operator=(const ConstSubmatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use DEFAULT definitions of the following six functions.
//     void myMulByRow(vec& lhs, const vec& v) const override;
//     void myMulByCol(vec& lhs, const vec& v) const override;
//     bool myIsZeroRow(long i) const override;
//     bool myIsZeroCol(long j) const override;
//     void myDet(RingElem& d) const override;
//     long myRank() const override;

  private: // data members
    ConstMatrixView myM;
    std::vector<long> myRowTable;
    std::vector<long> myColTable;
  };


  /***************************************************************************/
  // Functions for ConstSubmatImpl

  ConstSubmatImpl::ConstSubmatImpl(ConstMatrixView M, const std::vector<long>& rows, const std::vector<long>& cols):
      myM(M),
      myRowTable(rows),
      myColTable(cols)
  {
    // Args already checked by pseudo-ctor submat, the only
    // function which calls this ctor
  }


  const ring& ConstSubmatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long ConstSubmatImpl::myNumRows() const
  {
    return len(myRowTable);
  }


  long ConstSubmatImpl::myNumCols() const
  {
    return len(myColTable);
  }


  RingElemAlias ConstSubmatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM(myRowTable[i], myColTable[j]);
  }


  ConstMatrixView submat(ConstMatrixView M, const std::vector<long>& rows, const vector<long>& cols)
  {
    const long LenRows = len(rows);
    const long LenCols = len(cols);
    for (long i=0; i < LenRows; ++i)
      M->myCheckRowIndex(rows[i], CoCoA_ERROR_CONTEXT);
    for (long j=0; j < LenCols; ++j)
      M->myCheckColIndex(cols[j], CoCoA_ERROR_CONTEXT);

    // Check all row indices are distinct.
    vector<bool> RowUsed(NumRows(M));
    for (long i=0; i < LenRows; ++i)
    {
      if (RowUsed[rows[i]])
        CoCoA_THROW_ERROR2(ERR::BadIndex, "repeated row indices");
      RowUsed[rows[i]] = true;
    }

    // Check all col indices are distinct.
    vector<bool> ColUsed(NumCols(M));
    for (long i=0; i < LenCols; ++i)
    {
      if (ColUsed[cols[i]])
        CoCoA_THROW_ERROR2(ERR::BadIndex, "repeated col indices");
      ColUsed[cols[i]] = true;
    }

    return ConstMatrixView(new ConstSubmatImpl(M, rows, cols));
  }


  /***************************************************************************/
  // SubmatImpl -- see only certain rows/cols of a matrix

  class SubmatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView submat(MatrixView M, const std::vector<long>& RowIndices, const std::vector<long>& ColIndices); // pseudo-ctor
    SubmatImpl(MatrixView M, const std::vector<long>& rows, const std::vector<long>& cols);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    SubmatImpl(const SubmatImpl&) = delete;
    SubmatImpl& operator=(const SubmatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use DEFAULT definitions of the following six functions.
//     void myMulByRow(vec& lhs, const vec& v) const override;
//     void myMulByCol(vec& lhs, const vec& v) const override;
//     bool myIsZeroRow(long i) const override;
//     bool myIsZeroCol(long j) const override;
//     void myDet(RingElem& d) const override;
//     long myRank() const override;

    // non-const member fns
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    MatrixView myM;
    std::vector<long> myRowTable;
    std::vector<long> myColTable;
  };


  /***************************************************************************/
  // Functions for SubmatImpl


  SubmatImpl::SubmatImpl(MatrixView M, const std::vector<long>& rows, const std::vector<long>& cols):
      myM(M),
      myRowTable(rows),
      myColTable(cols)
  {
    // Args already checked by pseudo-ctor submat, the only
    // function which calls this ctor
  }


  const ring& SubmatImpl::myRing() const
  {
    return RingOf(myM);
  }


  long SubmatImpl::myNumRows() const
  {
    return len(myRowTable);
  }


  long SubmatImpl::myNumCols() const
  {
    return len(myColTable);
  }


  RingElemAlias SubmatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM(myRowTable[i], myColTable[j]);
  }


  bool SubmatImpl::myIsWritable(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myM->myIsWritable(myRowTable[i], myColTable[j]);
  }


  RingElemRawPtr SubmatImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    return myM->myRawEntry(myRowTable[i], myColTable[j]);
  }


  void SubmatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], r);
  }

  void SubmatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], n);
  }

  void SubmatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], N);
  }


  void SubmatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i, j));
    myM->mySetEntry(myRowTable[i], myColTable[j], Q);
  }


  MatrixView submat(MatrixView M, const std::vector<long>& rows, const vector<long>& cols)
  {
    const long LenRows = len(rows);
    const long LenCols = len(cols);
    for (long i=0; i < LenRows; ++i)
      M->myCheckRowIndex(rows[i], CoCoA_ERROR_CONTEXT);
    for (long j=0; j < LenCols; ++j)
      M->myCheckColIndex(cols[j], CoCoA_ERROR_CONTEXT);

    // Check all row indices are distinct.
    for (long i=0; i < LenRows; ++i)
      for (long j=i+1; j < LenRows; ++j)
        if (rows[i] == rows[j])
          CoCoA_THROW_ERROR2(ERR::BadIndex, "repeated row indices");

    // Check all col indices are distinct.
    for (long i=0; i < LenCols; ++i)
      for (long j=i+1; j < LenCols; ++j)
        if (cols[i] == cols[j])
          CoCoA_THROW_ERROR2(ERR::BadIndex, "repeated col indices");

    return MatrixView(new SubmatImpl(M, rows, cols));
  }


  MatrixView ColMat(MatrixView M, long col)
  { return submat(M, LongRange(0, NumRows(M)-1), vector<long>(1, col)); }
  
  ConstMatrixView ColMat(ConstMatrixView M, long col)
  { return submat(M, LongRange(0, NumRows(M)-1), vector<long>(1, col)); }
  
  MatrixView RowMat(MatrixView M, long row)
  { return submat(M, vector<long>(1, row), LongRange(0, NumCols(M)-1) ); }
  
  ConstMatrixView RowMat(ConstMatrixView M, long row)
  { return submat(M, vector<long>(1, row), LongRange(0, NumCols(M)-1) ); }

  MatrixView FirstCols(MatrixView M, long ncols)
  { return submat(M, LongRange(0, NumRows(M)-1), LongRange(0, ncols-1) ); }
  
  ConstMatrixView FirstCols(ConstMatrixView M, long ncols)
  { return submat(M, LongRange(0, NumRows(M)-1), LongRange(0, ncols-1) ); }

  MatrixView FirstRows(MatrixView M, long nrows)
  { return submat(M, LongRange(0, nrows-1), LongRange(0, NumCols(M)-1) ); }
  
  ConstMatrixView FirstRows(ConstMatrixView M, long nrows)
  { return submat(M, LongRange(0, nrows-1), LongRange(0, NumCols(M)-1) ); }


  /****************************************************************************/

  class ConstColMatImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView ColMat(const std::vector<RingElem>& v); // pseudo-ctor
    ConstColMatImpl(const std::vector<RingElem>& v);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ConstColMatImpl(const ConstColMatImpl&) = delete;
    ConstColMatImpl& operator=(const ConstColMatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    const std::vector<RingElem>& myV; // NB this is a reference!!!
  };


  ConstColMatImpl::ConstColMatImpl(const std::vector<RingElem>& v):
      myR(owner(v[0])),
      myNumRowsValue(len(v)),
      myV(v)
  {
    CoCoA_ASSERT(!v.empty());
    CoCoA_ASSERT(HasUniqueOwner(v));
  }


  const ring& ConstColMatImpl::myRing() const
  { return myR; }


  long ConstColMatImpl::myNumRows() const
  { return myNumRowsValue; }


  long ConstColMatImpl::myNumCols() const
  { return 1; }


  RingElemAlias ConstColMatImpl::myEntry(long i, long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[i];
  }


  bool ConstColMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return IsZero(myV[i]);
  }

  bool ConstColMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myV[i])) return false;
    return true;
  }


  void ConstColMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(myNumRows() == 1);
    d = myV[0];
  }


  long ConstColMatImpl::myRank() const
  {
    if (myIsZeroCol(0)) return 0;
    return 1;
  }


  ConstMatrixView ColMat(const std::vector<RingElem>& v)
  {
    if (v.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return ConstMatrixView(new ConstColMatImpl(v));
  }

  /****************************************************************************/

  class ColMatImpl: public MatrixViewBase
  {
  private:
    friend MatrixView ColMat(std::vector<RingElem>& v); // pseudo-ctor
    ColMatImpl(std::vector<RingElem>& v);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ColMatImpl(const ColMatImpl&) = delete;
    ColMatImpl& operator=(const ColMatImpl&) = delete;

  public:  // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    // non-const member fns
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    std::vector<RingElem>& myV; // NB this is a reference!!!
  };


  ColMatImpl::ColMatImpl(std::vector<RingElem>& v):
      myR(owner(v[0])),
      myNumRowsValue(len(v)),
      myV(v)
  {
    CoCoA_ASSERT(!v.empty());
    CoCoA_ASSERT(HasUniqueOwner(v));
  }


  const ring& ColMatImpl::myRing() const
  { return myR; }


  long ColMatImpl::myNumRows() const
  { return myNumRowsValue; }


  long ColMatImpl::myNumCols() const
  { return 1; }


  RingElemAlias ColMatImpl::myEntry(long i, long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[i];
  }


  bool ColMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    return IsZero(myV[i]);
  }

  bool ColMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myV[i])) return false;
    return true;
  }


  void ColMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(myNumRows() == 1);
    d = myV[0];
  }


  long ColMatImpl::myRank() const
  {
    if (myIsZeroCol(0)) return 0;
    return 1;
  }


  bool ColMatImpl::myIsWritable(long i, long j) const
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr ColMatImpl::myRawEntry(long i, long j)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myV[i]);
  }


  void ColMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = r;
  }

  void ColMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = n;
  }

  void ColMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = N;
  }


  void ColMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[i] = Q;
  }


  MatrixView ColMat(std::vector<RingElem>& v)
  {
    if (v.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return MatrixView(new ColMatImpl(v));
  }


  /***************************************************************************/
  // Very simple definitions... but may not be very efficient.

  MatrixView RowMat(std::vector<RingElem>& v)
  { return transpose(ColMat(v)); }

  
  ConstMatrixView RowMat(const std::vector<RingElem>& v)
  { return transpose(ColMat(v)); }


  /***************************************************************************/
  // ConstDiagMatImpl -- diagonal matrices implemented in a pretty obvious manner

  class ConstDiagMatImpl: public ConstMatrixViewBase
  {
  public:
    ConstDiagMatImpl(const ring& R, const std::vector<RingElem>& DiagEntries);
    // Default dtor is OK

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
    bool IamAntiSymmetric() const override;
    bool IamDiagonal() const override  { return true; }
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

  private: // data members
    const ring myR;
    const long myDim;
    const std::vector<RingElem>& myDiagEntries; // NB reference to an external vector!!!
  };

  //----------------------------------------------------------------------

  bool IsConstDiagMatImpl(ConstMatrixView M)  
  {return dynamic_cast<ConstDiagMatImpl*>(M.mySmartPtr.myRawPtr()) != nullptr;}


  ConstDiagMatImpl::ConstDiagMatImpl(const ring& R, const std::vector<RingElem>& DiagEntries):
      myR(R),
      myDim(len(DiagEntries)),
      myDiagEntries(DiagEntries)
  {
#ifdef CoCoA_DEBUG
    for (long i=0; i < myDim; ++i)
      CoCoA_ASSERT(owner(myDiagEntries[i]) == myR);
#endif
  }


  const ring& ConstDiagMatImpl::myRing() const
  { return myR; }


  long ConstDiagMatImpl::myNumRows() const
  { return myDim; }


  long ConstDiagMatImpl::myNumCols() const
  { return myDim; }


  RingElemAlias ConstDiagMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == j)
      return myDiagEntries[i];
    else
      return zero(myR);
  }


  void ConstDiagMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    //BUG: Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? does this clobber iterators on lhs???
  }


  void ConstDiagMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    //BUG: Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? does this clobber iterators on lhs???
  }


  bool ConstDiagMatImpl::IamEqual(ConstMatrixView M) const
  {
    if (myRing() != RingOf(M)) return false;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (IsZeroMatImpl(M) || IsIdentityMatImpl(M) || IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "ConstDiagMatImpl::IamEqual" << std::endl;
      for (long i=0; i < myNumRows(); ++i)
        if (myDiagEntries[i] != M(i,i)) return false;
      return true;
    }
    return ConstMatrixViewBase::IamEqual(M);
  }


  bool ConstDiagMatImpl::IamAntiSymmetric() const
  {
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myDiagEntries[i])) return false;
    return true;
  }


  bool ConstDiagMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    return IsZero(myDiagEntries[i]);
  }


  bool ConstDiagMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return IsZero(myDiagEntries[j]);
  }


  void ConstDiagMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myR);
    // Correct also if myDim == 0
    RingElem ans = one(myR);
    for (long i=0; i < myDim; ++i)
      ans *= myDiagEntries[i];
    swap(d, ans);
  }


  // ??? better to use "count"?
  long ConstDiagMatImpl::myRank() const
  {
    CoCoA_ASSERT(IsIntegralDomain(myR));
    long rk=0;
    for (long i = 0; i < myDim; ++i)
      if (!IsZero(myDiagEntries[i]))
        ++rk;
    return rk;
  }


  ConstMatrixView DiagMat(const std::vector<RingElem>& DiagEntries)
  {
    if (DiagEntries.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(DiagEntries))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return ConstMatrixView(new ConstDiagMatImpl(owner(DiagEntries[0]), DiagEntries));
  }


  /***************************************************************************/
  // DiagMatImpl -- diagonal matrices implemented in a pretty obvious manner

  class DiagMatImpl: public MatrixViewBase
  {
  public:
    DiagMatImpl(const ring& R, std::vector<RingElem>& DiagEntries);
    // Default dtor is OK

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    bool IamEqual(ConstMatrixView M) const override;
    bool IamSymmetric() const override  {return true;}
    bool IamAntiSymmetric() const override;
    bool IamDiagonal() const override  {return true;}
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    // non-const functions
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const ring myR;
    const long myDim;
    std::vector<RingElem>& myDiagEntries; // NB reference to an external vector!!!
  };

  //----------------------------------------------------------------------

  bool IsDiagMatImpl(ConstMatrixView M)  
  { return dynamic_cast<DiagMatImpl*>(M.mySmartPtr.myRawPtr()) != nullptr; }


  DiagMatImpl::DiagMatImpl(const ring& R, std::vector<RingElem>& DiagEntries):
    myR(R),
    myDim(len(DiagEntries)),
    myDiagEntries(DiagEntries)
  {
#ifdef CoCoA_DEBUG
    for (long i=0; i < myDim; ++i)
      CoCoA_ASSERT(owner(myDiagEntries[i]) == myR);
#endif
  }


  const ring& DiagMatImpl::myRing() const
  { return myR; }


  long DiagMatImpl::myNumRows() const
  { return myDim; }


  long DiagMatImpl::myNumCols() const
  { return myDim; }


  RingElemAlias DiagMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == j)
      return myDiagEntries[i];
    else
      return zero(myR);
  }


  void DiagMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    // Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? does this clobber iterators on lhs???
  }


  void DiagMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    // Should also check owner of elements of v and lhs.
    vector<RingElem> ans(myDim, zero(myR));
    for (long i=0; i < myDim; ++i)
      ans[i] = v[i]*myDiagEntries[i];
    for (long i=0; i < myDim; ++i)
      swap(lhs[i], ans[i]);
//???    std::swap(lhs, ans); //??? would this clobber iterators on lhs???
  }


  bool DiagMatImpl::IamEqual(ConstMatrixView M) const
  {
    if (myRing() != RingOf(M)) return false;
    if (myNumRows() != NumRows(M)) return false;
    if (myNumCols() != NumCols(M)) return false;
    if (IsZeroMatImpl(M) || IsIdentityMatImpl(M) || IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
    {
      //      std::cout << "DiagMatImpl::IamEqual" << std::endl;
      for (long i=0; i < myNumRows(); ++i)
        if (myDiagEntries[i] != M(i,i)) return false;
      return true;
    }
    return ConstMatrixViewBase::IamEqual(M);
  }


  bool DiagMatImpl::IamAntiSymmetric() const
  {
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myDiagEntries[i])) return false;
    return true;
  }


  bool DiagMatImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    return IsZero(myDiagEntries[i]);
  }


  bool DiagMatImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return IsZero(myDiagEntries[j]);
  }


  void DiagMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myR);
    // Loop below is correct even when myDim == 0
    RingElem ans = one(myR);
    for (long i=0; i < myDim; ++i)
      ans *= myDiagEntries[i];
    swap(d, ans);
  }


  // ??? better to use "std::count"???
  long DiagMatImpl::myRank() const
  {
    CoCoA_ASSERT(IsIntegralDomain(myR));
    long rk=0;
    for (long i = 0; i < myDim; ++i)
      if (!IsZero(myDiagEntries[i]))
        ++rk;
    return rk;
  }


  bool DiagMatImpl::myIsWritable(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return i==j;
  }

  RingElemRawPtr DiagMatImpl::myRawEntry(long i, long j)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    return raw(myDiagEntries[i]);
  }


  void DiagMatImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = r;
  }


  void DiagMatImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = n;
  }


  void DiagMatImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = N;
  }


  void DiagMatImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(myIsWritable(i,j));
    myDiagEntries[i] = Q;
  }


  MatrixView DiagMat(std::vector<RingElem>& DiagEntries)
  {
    if (DiagEntries.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(DiagEntries))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return MatrixView(new DiagMatImpl(owner(DiagEntries[0]), DiagEntries));
  }

  /****************************************************************************/

  class MatByRowsImpl: public MatrixViewBase
  {
  private:
    friend MatrixView MatByRows(const MachineInt& nrows, const MachineInt& ncols, std::vector<RingElem>& v); // pseudo-ctor
    MatByRowsImpl(long nrows, long ncols, std::vector<RingElem>& v);
    static std::vector<RingElem>& CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    MatByRowsImpl(const MatByRowsImpl&) = delete;
    MatByRowsImpl& operator=(const MatByRowsImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
/// USE DEFAULT?    void myDet(RingElem& d) const override;
/// USE DEFAULT?    long myRank() const override;

    // non-const member fns
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const long myNumRowsValue;
    const long myNumColsValue;
    const ring myR;
    std::vector<RingElem>& myV; // NB this is a reference!!!
  private: // impl detail
    long myIndex(long i, long j) const { return i*myNumColsValue+j; }

  };


  std::vector<RingElem>& MatByRowsImpl::CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v)
  {
    if (nrows < 1 || ncols < 1)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    const long n = len(v);
    // condition just tests if n = nrows*ncols but avoids overflow
    if (n/nrows != ncols || n != nrows*ncols)  CoCoA_THROW_ERROR1(ERR::BadArraySize);
//???    if (v.empty()) CoCoA_THROW_ERROR2(ERR::Empty, "MatByRows ctor");
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return v;
  }


  MatByRowsImpl::MatByRowsImpl(long nrows, long ncols, std::vector<RingElem>& v):
    myNumRowsValue(nrows),
    myNumColsValue(ncols),
    myR(owner(CtorArgCheck(nrows, ncols,v)[0])),
    myV(v)
  {}


  const ring& MatByRowsImpl::myRing() const
  { return myR; } // return owner(myV[0]);


  long MatByRowsImpl::myNumRows() const
  { return myNumRowsValue; }


  long MatByRowsImpl::myNumCols() const
  { return myNumColsValue; }


  RingElemAlias MatByRowsImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[myIndex(i,j)];
  }


  bool MatByRowsImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    const long start = myIndex(i,0);
    const long end = myIndex(i,myNumColsValue-1);
    for (long j=start; j <= end; ++j)
      if (!IsZero(myV[j])) return false;
    return true;
  }

  bool MatByRowsImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myV[myIndex(i,j)])) return false;
    return true;
  }


  bool MatByRowsImpl::myIsWritable(long i, long j) const
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr MatByRowsImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myV[myIndex(i,j)]);
  }


  void MatByRowsImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = r;
  }

  void MatByRowsImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = n;
  }

  void MatByRowsImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = N;
  }


  void MatByRowsImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = Q;
  }


  MatrixView MatByRows(const MachineInt& nrows, const MachineInt& ncols, std::vector<RingElem>& v)
  {
    if (IsNegative(nrows) || !IsSignedLong(nrows)) CoCoA_THROW_ERROR2(ERR::OutOfRange, "num rows");
    if (IsNegative(ncols) || !IsSignedLong(ncols)) CoCoA_THROW_ERROR2(ERR::OutOfRange, "num cols");
    // We let ctor for MatByRowsImpl do the arg checking...
    return MatrixView(new MatByRowsImpl(AsSignedLong(nrows),AsSignedLong(ncols),v));
  }


  /****************************************************************************/

  class MatByRowsImpl2: public MatrixViewBase
  {
  private:
    friend MatrixView MatByRows(std::vector< std::vector<RingElem> >& vv); // pseudo-ctor
    MatByRowsImpl2(std::vector< std::vector<RingElem> >& vv);
    static std::vector< std::vector<RingElem> >& CtorArgCheck(std::vector< std::vector<RingElem> >& vv);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    MatByRowsImpl2(const MatByRowsImpl&) = delete;
    MatByRowsImpl2& operator=(const MatByRowsImpl2&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
/// USE DEFAULT?    void myDet(RingElem& d) const override;
/// USE DEFAULT?    long myRank() const override;

    // non-const member fns
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const ring myR;
    const long myNumRowsValue;
    const long myNumColsValue;
    std::vector< std::vector<RingElem> >& myVV; // NB this is a reference!!!
  };


  std::vector< std::vector<RingElem> >& MatByRowsImpl2::CtorArgCheck(std::vector< std::vector<RingElem> >& vv)
  {
    if (vv.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    const long nrows = len(vv);
    for (long i=0; i < nrows; ++i)
    {
      if (len(vv[i]) != len(vv[0]) )  CoCoA_THROW_ERROR1(ERR::IncompatDims);
      if (!HasUniqueOwner(vv[i]) || owner(vv[i][0]) != owner(vv[0][0]))
        CoCoA_THROW_ERROR1(ERR::MixedRings);
    }
    return vv;
  }


  MatByRowsImpl2::MatByRowsImpl2(std::vector< std::vector<RingElem> >& vv):
    myR(owner(CtorArgCheck(vv)[0][0])),
    myNumRowsValue(len(vv)),
    myNumColsValue(len(vv[0])),
    myVV(vv)
  {}
  

  const ring& MatByRowsImpl2::myRing() const
  { return myR; }


  long MatByRowsImpl2::myNumRows() const
  { return myNumRowsValue; }


  long MatByRowsImpl2::myNumCols() const
  { return myNumColsValue; }


  RingElemAlias MatByRowsImpl2::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myVV[i][j];
  }


  bool MatByRowsImpl2::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    for (long j=0; j < myNumColsValue; ++j)
      if (!IsZero(myVV[i][j])) return false;
    return true;
  }

  bool MatByRowsImpl2::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    for (long i=0; i < myNumRows(); ++i)
      if (!IsZero(myVV[i][j])) return false;
    return true;
  }


  bool MatByRowsImpl2::myIsWritable(long i, long j) const
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr MatByRowsImpl2::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myVV[i][j]);
  }


  void MatByRowsImpl2::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = r;
  }

  void MatByRowsImpl2::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = n;
  }

  void MatByRowsImpl2::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = N;
  }


  void MatByRowsImpl2::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myVV[i][j] = Q;
  }


  MatrixView MatByRows(std::vector< std::vector<RingElem> >& vv)
  {
    // We let ctor for MatByRowsImpl2 do the arg checking...
    return MatrixView(new MatByRowsImpl2(vv));
  }


  /****************************************************************************/

  class MatByColsImpl: public MatrixViewBase
  {
  private:
    friend MatrixView MatByCols(const MachineInt& nrows, const MachineInt& ncols, std::vector<RingElem>& v); // pseudo-ctor
    MatByColsImpl(long nrows, long ncols, std::vector<RingElem>& v);
    static std::vector<RingElem>& CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    MatByColsImpl(const MatByRowsImpl&) = delete;
    MatByColsImpl& operator=(const MatByRowsImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//     void myMulByRow(vec& lhs, const vec& v) const override;
//     void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
/// USE DEFAULT?    void myDet(RingElem& d) const override;
/// USE DEFAULT?    long myRank() const override;

    // non-const member fns
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const long myNumRowsValue;
    const long myNumColsValue;
    const ring myR;
    std::vector<RingElem>& myV; // NB this is a reference!!!
  private: // impl detail
    long myIndex(long i, long j) const { return j*myNumRowsValue+i; }

  };


  std::vector<RingElem>& MatByColsImpl::CtorArgCheck(long nrows, long ncols, std::vector<RingElem>& v)
  {
    if (nrows < 1 || ncols < 1)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    const long n = len(v);
    // Check that n = nrows*ncols (avoiding overflow)
    if (n/nrows != ncols || n != nrows*ncols)  CoCoA_THROW_ERROR1(ERR::IncompatDims);
//???    if (v.empty()) CoCoA_THROW_ERROR2(ERR::Empty, "MatByCols ctor");
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return v;
  }


  MatByColsImpl::MatByColsImpl(long nrows, long ncols, std::vector<RingElem>& v):
    myNumRowsValue(nrows),
    myNumColsValue(ncols),
    myR(owner(CtorArgCheck(nrows, ncols,v)[0])),
    myV(v)
  {}


  const ring& MatByColsImpl::myRing() const
  { return myR; }// return owner(myV[0]);


  long MatByColsImpl::myNumRows() const
  { return myNumRowsValue; }


  long MatByColsImpl::myNumCols() const
  { return myNumColsValue; }


  RingElemAlias MatByColsImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return myV[myIndex(i,j)];
  }


  bool MatByColsImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i >= 0 && i < myNumRows());
    for (long j=0; j < myNumCols(); ++j)
      if (!IsZero(myV[myIndex(i,j)])) return false;
    return true;
  }

  bool MatByColsImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    const long start = myIndex(0,j);
    const long end = myIndex(myNumRowsValue-1,j);
    for (long i=start; i <= end; ++i)
      if (!IsZero(myV[i])) return false;
    return true;
  }


  bool MatByColsImpl::myIsWritable(long i, long j) const
  {
    (void)(i); (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return true;
  }

  RingElemRawPtr MatByColsImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    return raw(myV[myIndex(i,j)]);
  }


  void MatByColsImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(owner(r) == myR);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = r;
  }

  void MatByColsImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = n;
  }

  void MatByColsImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = N;
  }


  void MatByColsImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    myV[myIndex(i,j)] = Q;
  }


  MatrixView MatByCols(const MachineInt& nrows, const MachineInt& ncols, std::vector<RingElem>& v)
  {
    if (IsNegative(nrows) || !IsSignedLong(nrows)) CoCoA_THROW_ERROR2(ERR::OutOfRange, "num rows");
    if (IsNegative(ncols) || !IsSignedLong(ncols)) CoCoA_THROW_ERROR2(ERR::OutOfRange, "num cols");
    // We let ctor for MatByColsImpl do the arg checking...
    return MatrixView(new MatByColsImpl(AsSignedLong(nrows),AsSignedLong(ncols),v));
  }


  /***************************************************************************/

  namespace // anonymous
  {
    
    class ZeroBlockImpl;
    class ZeroBlock: public MatrixView
    {
      // Data member inherited from ConstMatrixView via MatrixView
    public:
      explicit ZeroBlock(const ring& R);
      // Assignment disabled.
      // Default dtor is OK.
      ZeroBlockImpl* operator->() const;
    };

    class ZeroBlockImpl: public MatrixViewBase
    {
    public:
      ZeroBlockImpl(const ring& R): myR(R) {}
      // default copy ctor OK, default dtor OK; default assign not OK because a data field is a reference

    public: // functions every matrix type must implement
      const ring& myRing() const override  { return myR; }
      long myNumRows() const override  { return numeric_limits<long>::max(); } // called by (Const)MatrixViewBase::myCheckRowIndex
      long myNumCols() const override  { return numeric_limits<long>::max(); } // called by (Const)MatrixViewBase::myCheckColIndex
      RingElemAlias myEntry(long, long) const override  { return zero(myR); }
      bool myIsZeroRow(long) const override  { return true; }
      bool myIsZeroCol(long) const override  { return true; }

      bool myIsWritable(long, long) const override  { return false; }
      RingElemRawPtr myRawEntry(long, long) override  { CoCoA_THROW_ERROR1(ERR::ConstMatEntry); RingElem JUNK;return raw(JUNK); /*just to keep compiler quiet*/ }
      void mySetEntry(long, long, ConstRefRingElem) override  { CoCoA_ASSERT(0 && "ERR::ConstMatEntry"); }
      void mySetEntry(long, long, const MachineInt&) override  { CoCoA_ASSERT(0 && "ERR::ConstMatEntry"); }
      void mySetEntry(long, long, const BigInt&) override  { CoCoA_ASSERT(0 && "ERR::ConstMatEntry"); }
      void mySetEntry(long, long, const BigRat&) override  { CoCoA_ASSERT(0 && "ERR::ConstMatEntry"); }
    private: // data member
      const ring& myR;
    };

    ZeroBlock::ZeroBlock(const ring& R): MatrixView(new ZeroBlockImpl(R)) {}


    class ConstBlockMat2x2Impl: public ConstMatrixViewBase
    {
    public:
//    friend ConstMatrixView BlockMat2x2(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D); // pseudo-ctor
      ConstBlockMat2x2Impl(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D);
      ConstBlockMat2x2Impl(ConstMatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, ConstMatrixView D);
      ConstBlockMat2x2Impl(ZeroBlockIndicator, ConstMatrixView B, ConstMatrixView C, ZeroBlockIndicator);
      // default dtor works fine
    public: // disable default copy ctor and assignment
      ConstBlockMat2x2Impl(const ConstBlockMat2x2Impl&) = delete;
      ConstBlockMat2x2Impl& operator=(const ConstBlockMat2x2Impl&) = delete;

    public: // functions every matrix type must implement
      typedef std::vector<RingElem> vec;
      const ring& myRing() const override;
      long myNumRows() const override;
      long myNumCols() const override;
      RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
      bool myIsZeroRow(long i) const override;
      bool myIsZeroCol(long j) const override;
// Use default defns of the following two functions
//    void myDet(RingElem& d) const override;
//    long myRank() const override;

    private: // data members
      const ring& myR;
      long myNumRowsValue;
      long myNumColsValue;
      long myRowSplit;
      long myColSplit;
      ZeroBlock myZeroBlock; // useful only if ctor was called with "zeroes"
      ConstMatrixView myA, myB, myC, myD;
    };


    /***************************************************************************/
    // Functions for ConstBlockMat2x2Impl -- a matrix of the form  (A B)
    //                                                             (C D)
    // This implementation requires that the submatrices correspond
    // to two straight cuts -- several changes would be needed to handle
    // a more general assembly of submatrices.

    ConstBlockMat2x2Impl::ConstBlockMat2x2Impl(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D):
        myR(RingOf(A)),
        myNumRowsValue(NumRows(A)+NumRows(C)),
        myNumColsValue(NumCols(A)+NumCols(B)),
        myRowSplit(NumRows(A)),
        myColSplit(NumCols(A)),
        myZeroBlock(myR),
        myA(A),
        myB(B),
        myC(C),
        myD(D)
    {
      // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
      // only function which calls this ctor.
    }


    ConstBlockMat2x2Impl::ConstBlockMat2x2Impl(ConstMatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, ConstMatrixView D):
        myR(RingOf(A)),
        myNumRowsValue(NumRows(A)+NumRows(D)),
        myNumColsValue(NumCols(A)+NumCols(D)),
        myRowSplit(NumRows(A)),
        myColSplit(NumCols(A)),
        myZeroBlock(myR),
        myA(A),
        myB(myZeroBlock),
        myC(myZeroBlock),
        myD(D)
    {
      // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
      // only function which calls this ctor.
    }

    ConstBlockMat2x2Impl::ConstBlockMat2x2Impl(ZeroBlockIndicator, ConstMatrixView B, ConstMatrixView C, ZeroBlockIndicator):
        myR(RingOf(B)),
        myNumRowsValue(NumRows(B)+NumRows(C)),
        myNumColsValue(NumCols(B)+NumCols(C)),
        myRowSplit(NumRows(B)),
        myColSplit(NumCols(C)),
        myZeroBlock(myR),
        myA(myZeroBlock),
        myB(B),
        myC(C),
        myD(myZeroBlock)
    {
      // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
      // only function which calls this ctor.
    }


    const ring& ConstBlockMat2x2Impl::myRing() const
    { return myR; }


    long ConstBlockMat2x2Impl::myNumRows() const
    { return myNumRowsValue; }


    long ConstBlockMat2x2Impl::myNumCols() const
    { return myNumColsValue; }


    RingElemAlias ConstBlockMat2x2Impl::myEntry(long i, long j) const
    {
      CoCoA_ASSERT(i >= 0 && j >= 0);
      CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
      if (i < myRowSplit && j < myColSplit)
        return myA(i, j);
      if (i < myRowSplit)
        return myB(i, j-myColSplit);
      if (j < myColSplit)
        return myC(i-myRowSplit, j);
      return myD(i-myRowSplit, j-myColSplit);
    }


    bool ConstBlockMat2x2Impl::myIsZeroRow(long i) const
    {
      CoCoA_ASSERT(i < myNumRows());
      if (i < myRowSplit)
        return myA->myIsZeroRow(i) && myB->myIsZeroRow(i);
      return myC->myIsZeroRow(i-myRowSplit) && myD->myIsZeroRow(i-myRowSplit);
    }


    bool ConstBlockMat2x2Impl::myIsZeroCol(long j) const
    {
      CoCoA_ASSERT(j >= 0 && j < myNumCols());
      if (j < myColSplit)
        return myA->myIsZeroCol(j) && myC->myIsZeroCol(j);
      return myB->myIsZeroCol(j-myColSplit) && myD->myIsZeroCol(j-myColSplit);
    }

  } // end of anonymous namespace

  
  ConstMatrixView BlockMat2x2(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D)
  {
    {
      // Check that the rings are all the same
      const ring R = RingOf(A);
      if (RingOf(B) != R || RingOf(C) != R || RingOf(D) != R)
        CoCoA_THROW_ERROR1(ERR::MixedRings);
    }
    // Check that the matrix sizes are suitably matched.
    if (NumRows(A) != NumRows(B) ||
        NumRows(C) != NumRows(D) ||
        NumCols(A) != NumCols(C) ||
        NumCols(B) != NumCols(D))
      CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    if (NumCols(A) > numeric_limits<long>::max() -NumCols(B)) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "num cols");
    if (NumRows(A) > numeric_limits<long>::max() -NumRows(C)) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "num rows");
    return ConstMatrixView(new ConstBlockMat2x2Impl(A, B, C, D));
  }


  ConstMatrixView BlockMat2x2(ConstMatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, ConstMatrixView D)
  {
    // Check that the rings are the same
    if (RingOf(A) != RingOf(D))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return ConstMatrixView(new ConstBlockMat2x2Impl(A, zeroes, zeroes, D));
  }

  
  ConstMatrixView BlockMat2x2(ZeroBlockIndicator, ConstMatrixView B, ConstMatrixView C,  ZeroBlockIndicator)
  {
    // Check that the rings are the same
    if (RingOf(B) != RingOf(C))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return ConstMatrixView(new ConstBlockMat2x2Impl(zeroes, B, C, zeroes));
  }


/***************************************************************************/

  class BlockMat2x2Impl: public MatrixViewBase
  {
  public:
//    friend MatrixView BlockMat2x2(MatrixView A, MatrixView B, MatrixView C, MatrixView D); // pseudo-ctor
    BlockMat2x2Impl(MatrixView A, MatrixView B, MatrixView C, MatrixView D);
    BlockMat2x2Impl(MatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, MatrixView D);
    BlockMat2x2Impl(ZeroBlockIndicator, MatrixView B, MatrixView C, ZeroBlockIndicator);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    BlockMat2x2Impl(const BlockMat2x2Impl&) = delete;
    BlockMat2x2Impl& operator=(const BlockMat2x2Impl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
// Use default defns of the following two functions
//    void myDet(RingElem& d) const override;
//    long myRank() const override;

    // Non-const member functions
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    long myRowSplit;
    long myColSplit;
    ZeroBlock myZeroBlock; // useful only if ctor was called with "zeroes"
    MatrixView myA, myB, myC, myD;
  };


  /***************************************************************************/
  // Functions for BlockMat2x2Impl -- a matrix of the form  (A B)
  //                                                        (C D)
  // This implementation requires that the submatrices correspond
  // to two straight cuts -- several changes would be needed to handle
  // a more general assembly of submatrices.

  BlockMat2x2Impl::BlockMat2x2Impl(MatrixView A, MatrixView B, MatrixView C, MatrixView D):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(C)),
      myNumColsValue(NumCols(A)+NumCols(B)),
      myRowSplit(NumRows(A)),
      myColSplit(NumCols(A)),
      myZeroBlock(myR),
      myA(A),
      myB(B),
      myC(C),
      myD(D)
  {
    // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
    // only function which calls this ctor.
  }

  BlockMat2x2Impl::BlockMat2x2Impl(MatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, MatrixView D):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(D)),
      myNumColsValue(NumCols(A)+NumCols(D)),
      myRowSplit(NumRows(A)),
      myColSplit(NumCols(A)),
      myZeroBlock(myR),
      myA(A),
      myB(myZeroBlock),
      myC(myZeroBlock),
      myD(D)
  {
    // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
    // only function which calls this ctor.
  }

  BlockMat2x2Impl::BlockMat2x2Impl(ZeroBlockIndicator, MatrixView B, MatrixView C, ZeroBlockIndicator):
      myR(RingOf(B)),
      myNumRowsValue(NumRows(B)+NumRows(C)),
      myNumColsValue(NumCols(B)+NumCols(C)),
      myRowSplit(NumRows(B)),
      myColSplit(NumCols(C)),
      myZeroBlock(myR),
      myA(myZeroBlock),
      myB(B),
      myC(C),
      myD(myZeroBlock)
  {
    // Args are checked for compatibility by pseudo-ctor BlockMat2x2, the
    // only function which calls this ctor.
  }


  const ring& BlockMat2x2Impl::myRing() const
  { return myR; }


  long BlockMat2x2Impl::myNumRows() const
  { return myNumRowsValue; }


  long BlockMat2x2Impl::myNumCols() const
  { return myNumColsValue; }


  RingElemAlias BlockMat2x2Impl::myEntry(long i, long j) const
  {
// This one-liner doesn't work:    return ConstRefRingElem(myR, myRawEntry(i,j));
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < myRowSplit && j < myColSplit)
      return myA(i, j);
    if (i < myRowSplit)
      return myB(i, j-myColSplit);
    if (j < myColSplit)
      return myC(i-myRowSplit, j);
    return myD(i-myRowSplit, j-myColSplit);
  }


  bool BlockMat2x2Impl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < myRowSplit)
      return myA->myIsZeroRow(i) && myB->myIsZeroRow(i);
    return myC->myIsZeroRow(i-myRowSplit) && myD->myIsZeroRow(i-myRowSplit);
  }


  bool BlockMat2x2Impl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    if (j < myColSplit)
      return myA->myIsZeroCol(j) && myC->myIsZeroCol(j);
    return myB->myIsZeroCol(j-myColSplit) && myD->myIsZeroCol(j-myColSplit);
  }


  bool BlockMat2x2Impl::myIsWritable(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (j >= myColSplit)
    {
      if (i >= myRowSplit)
        return myD->myIsWritable(i-myRowSplit, j-myColSplit);
      else
        return myB->myIsWritable(i, j-myColSplit);
    }
    if (i >= myRowSplit)
      return myC->myIsWritable(i-myRowSplit, j);
    else
      return myA->myIsWritable(i, j);
  }


  RingElemRawPtr BlockMat2x2Impl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < myRowSplit && j < myColSplit)
      return myA->myRawEntry(i, j);
    if (i < myRowSplit)
      return myB->myRawEntry(i, j-myColSplit);
    if (j < myColSplit)
      return myC->myRawEntry(i-myRowSplit, j);
    return myD->myRawEntry(i-myRowSplit, j-myColSplit);
  }


  void BlockMat2x2Impl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myR == owner(r));
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), raw(r)); // myRefEntry(i,j) = r;
  }

  void BlockMat2x2Impl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), n); //    myRefEntry(i,j) = n;
  }

  void BlockMat2x2Impl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), N); //    myRefEntry(i,j) = N;
  }


  void BlockMat2x2Impl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), Q); //    myRefEntry(i,j) = Q;
  }


  MatrixView BlockMat2x2(MatrixView A, MatrixView B, MatrixView C, MatrixView D)
  {
    const ring R = RingOf(A);
    if (RingOf(B) != R || RingOf(C) != R || RingOf(D) != R)
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    // Check that the matrix sizes are suitably matched.
    if (NumRows(A) != NumRows(B) ||
        NumRows(C) != NumRows(D) ||
        NumCols(A) != NumCols(C) ||
        NumCols(B) != NumCols(D))
      CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    if (NumCols(A) > numeric_limits<long>::max() -NumCols(B)) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "num col");
    if (NumRows(A) > numeric_limits<long>::max() -NumRows(C)) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "num row");
    return MatrixView(new BlockMat2x2Impl(A, B, C, D));
  }

  
  MatrixView BlockMat2x2(MatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, MatrixView D)
  {
    // Check that the rings are the same
    if (RingOf(A) != RingOf(D))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return MatrixView(new BlockMat2x2Impl(A, zeroes,   zeroes, D));
  }

  
  MatrixView BlockMat2x2(ZeroBlockIndicator, MatrixView B, MatrixView C,  ZeroBlockIndicator)
  {
    if (RingOf(B) != RingOf(C))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return MatrixView(new BlockMat2x2Impl(zeroes, B,   C, zeroes));
  }


  /****************************************************************************/

  class ConstConcatVerImpl: public ConstMatrixViewBase
  {
  private:
    friend ConstMatrixView ConcatVer(ConstMatrixView A, ConstMatrixView B); // pseudo-ctor
    ConstConcatVerImpl(ConstMatrixView A, ConstMatrixView B);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ConstConcatVerImpl(const ConstConcatVerImpl&) = delete;
    ConstConcatVerImpl& operator=(const ConstConcatVerImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
// Use default defns of the following two functions
//    void myDet(RingElem& d) const override;
//    long myRank() const override;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    ConstMatrixView myA, myB;
  };



  ConstConcatVerImpl::ConstConcatVerImpl(ConstMatrixView A, ConstMatrixView B):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(B)),
      myNumColsValue(NumCols(A)),
      myA(A),
      myB(B)
  {
    // Args are checked for compatibility by ConcatVer, the
    // only function which calls this ctor.
  }


  const ring& ConstConcatVerImpl::myRing() const
  { return myR; }


  long ConstConcatVerImpl::myNumRows() const
  { return myNumRowsValue; }


  long ConstConcatVerImpl::myNumCols() const
  { return myNumColsValue; }


  RingElemAlias ConstConcatVerImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA))
      return myA(i, j);
    return myB(i-NumRows(myA), j);
  }


  bool ConstConcatVerImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < NumRows(myA))
      return myA->myIsZeroRow(i);
    return myB->myIsZeroRow(i-NumRows(myA));
  }


  bool ConstConcatVerImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myA->myIsZeroCol(j) && myB->myIsZeroCol(j);
  }


  ConstMatrixView ConcatVer(ConstMatrixView A, ConstMatrixView B)
  {
    if (RingOf(A) != RingOf(B))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumCols(A) != NumCols(B))  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    if (NumRows(A) > numeric_limits<long>::max() -NumRows(B)) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "num rows");
    return ConstMatrixView(new ConstConcatVerImpl(A, B));
  }


  /***************************************************************************/

  class ConcatVerImpl: public MatrixViewBase
  {
  private:
    friend MatrixView ConcatVer(MatrixView A, MatrixView B); // pseudo-ctor
    ConcatVerImpl(MatrixView A, MatrixView B);
    // default dtor works fine
  public: // disable default copy ctor and assignment
    ConcatVerImpl(const ConcatVerImpl&) = delete;
    ConcatVerImpl& operator=(const ConcatVerImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
// Use default defns of the following two functions
//    void myMulByRow(vec& lhs, const vec& v) const override;
//    void myMulByCol(vec& lhs, const vec& v) const override;
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
// Use default defns of the following two functions
//    void myDet(RingElem& d) const override;
//    long myRank() const override;

    // Non-const member functions
    bool myIsWritable(long i, long j) const override;
    RingElemRawPtr myRawEntry(long i, long j) override;
    void mySetEntry(long i, long j, ConstRefRingElem r) override;
    void mySetEntry(long i, long j, const MachineInt& n) override;
    void mySetEntry(long i, long j, const BigInt& N) override;
    void mySetEntry(long i, long j, const BigRat& Q) override;

  private: // data members
    const ring myR;
    long myNumRowsValue;
    long myNumColsValue;
    MatrixView myA, myB;
  };


  /***************************************************************************/
  // Functions for ConcatVerImpl -- a matrix of the form  (A)
  //                                                      (B)
  // This implementation requires that the submatrices correspond
  // to two straight cuts -- several changes would be needed to handle
  // a more general assembly of submatrices.

  ConcatVerImpl::ConcatVerImpl(MatrixView A, MatrixView B):
      myR(RingOf(A)),
      myNumRowsValue(NumRows(A)+NumRows(B)),
      myNumColsValue(NumCols(A)),
      myA(A),
      myB(B)
  {
    // Args are checked for compatibility by ConcatVer, the
    // only function which calls this ctor.
  }


  const ring& ConcatVerImpl::myRing() const
  { return myR; }


  long ConcatVerImpl::myNumRows() const
  { return myNumRowsValue; }


  long ConcatVerImpl::myNumCols() const
  { return myNumColsValue; }


  RingElemAlias ConcatVerImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i < NumRows(myA))
      return myA(i, j);
    return myB(i-NumRows(myA), j);
  }


  bool ConcatVerImpl::myIsZeroRow(long i) const
  {
    CoCoA_ASSERT(i < myNumRows());
    if (i < NumRows(myA))
      return myA->myIsZeroRow(i);
    return myB->myIsZeroRow(i-NumRows(myA));
  }


  bool ConcatVerImpl::myIsZeroCol(long j) const
  {
    CoCoA_ASSERT(j >= 0 && j < myNumCols());
    return myA->myIsZeroCol(j) && myB->myIsZeroCol(j);
  }


  bool ConcatVerImpl::myIsWritable(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i >= NumRows(myA))
      return myB->myIsWritable(i-NumRows(myA), j);
    else
      return myA->myIsWritable(i, j);
  }


  RingElemRawPtr ConcatVerImpl::myRawEntry(long i, long j)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    if (i >= NumRows(myA))
      return myB->myRawEntry(i-NumRows(myA), j);
    else
      return myA->myRawEntry(i, j);
  }


  void ConcatVerImpl::mySetEntry(long i, long j, ConstRefRingElem r)
  {
    CoCoA_ASSERT(myR == owner(r));
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), raw(r)); //    myRefEntry(i,j) = r;
  }

  void ConcatVerImpl::mySetEntry(long i, long j, const MachineInt& n)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), n); //    myRefEntry(i,j) = n;
  }

  void ConcatVerImpl::mySetEntry(long i, long j, const BigInt& N)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), N); //    myRefEntry(i,j) = N;
  }


  void ConcatVerImpl::mySetEntry(long i, long j, const BigRat& Q)
  {
    CoCoA_ASSERT(myIsWritable(i,j));
    myR->myAssign(myRawEntry(i,j), Q); //    myRefEntry(i,j) = Q;
  }


  MatrixView ConcatVer(MatrixView A, MatrixView B)
  {
    // Check that the base rings are the same.
    if (RingOf(A) != RingOf(B))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumCols(A) != NumCols(B))  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    return MatrixView(new ConcatVerImpl(A, B));
  }


  /***************************************************************************/
  // ConcatHor implemented in terms of ConcatVer.
  // Very simple implementations... but may not be so efficient.

  ConstMatrixView ConcatHor(ConstMatrixView M1, ConstMatrixView M2)
  {
    if (RingOf(M1) != RingOf(M2))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumRows(M1) != NumRows(M2))  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    return transpose(ConcatVer(transpose(M1), transpose(M2)));
  }

  
  MatrixView ConcatHor(MatrixView M1, MatrixView M2)
  {
    if (RingOf(M1) != RingOf(M2))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumRows(M1) != NumRows(M2))  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    return transpose(ConcatVer(transpose(M1), transpose(M2)));
  }


  /***************************************************************************/

  MatrixView ConcatDiag(MatrixView M1, MatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return BlockMat2x2(M1, zeroes,   zeroes, M2);
  }

  
  ConstMatrixView ConcatDiag(ConstMatrixView M1, ConstMatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return BlockMat2x2(M1, zeroes,   zeroes, M2);
  }


  MatrixView ConcatAntiDiag(MatrixView M1, MatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return BlockMat2x2(zeroes, M1,   M2, zeroes);
  }

  ConstMatrixView ConcatAntiDiag(ConstMatrixView M1, ConstMatrixView M2)
  {
    ring R = RingOf(M1);
    if (RingOf(M2) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    return BlockMat2x2(zeroes, M1,   M2, zeroes);
  }



} // End of namespace CoCoA
