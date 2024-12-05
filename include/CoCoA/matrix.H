#ifndef CoCoA_matrix_H
#define CoCoA_matrix_H

//   Copyright (c)  2005,2008,2015  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigInt.H"
#include "CoCoA/ring.H"
#include "CoCoA/error.H"

#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;


namespace CoCoA
{

  class MachineInt; // fwd decl -- defined in MachineInt.H
  class OpenMathOutput; // fwd decl -- defined in OpenMath.H

  // Operations on vectors -- should be ModuleElems?  (of a FreeModule?)
  bool IsZero(const std::vector<RingElem>& v);

  template <typename T>
  bool IsRectangular(const std::vector< std::vector<T> >& vv)
  {
    const std::size_t nrows = vv.size();
    if (nrows == 0) return true;
    const std::size_t ncols = vv[0].size();
    for (std::size_t n=1; n < nrows; ++n)
      if (vv[n].size() != ncols) return false;
    return  true;
  }

  class matrix;               // forward declaration
  class MatrixBase;           // forward declaration
  class MatrixViewBase;       // forward declaration
  class ConstMatrixViewBase;  // forward declaration
  class ConstMatrixBase;      // forward declaration

  class ConstMatrixView
  {
  public:
    explicit ConstMatrixView(const ConstMatrixViewBase* MatPtr);
    // Default copy ctor is OK.
    // Default assignment is OK.
    // Default dtor is OK.
    const ConstMatrixViewBase* operator->() const;
    RingElemAlias operator()(const MachineInt& i, const MachineInt& j) const;
    bool operator==(ConstMatrixView M) const;
    bool operator!=(ConstMatrixView M) const { return !(*this==M); }

  protected: // data member
    SmartPtrIRC<ConstMatrixViewBase> mySmartPtr;
    friend void swap(matrix& M1, matrix& M2);
    friend bool IsZeroMatImpl(const ConstMatrixView& M);
    friend bool IsIdentityMatImpl(const ConstMatrixView& M);
    friend bool IsDiagMatImpl(ConstMatrixView M);
    friend bool IsConstDiagMatImpl(ConstMatrixView M);
  };


  // Entries are constant!
  class ConstMatrix: public ConstMatrixView
  {
    // Data member inherited from ConstMatrixView
  public:
    explicit ConstMatrix(ConstMatrixBase* MatPtr);
    ConstMatrix(const ConstMatrix& M);
    // Default assignment is OK.
    // Default dtor is OK.
    const ConstMatrixBase* operator->() const;
  };


  class MatrixView: public ConstMatrixView
  {
    // Data member inherited from ConstMatrixView
  public:
    explicit MatrixView(MatrixViewBase* MatPtr);
    // Default copy ctor is OK.
    // Default assignment is OK.
    // Default dtor is OK.
    MatrixViewBase* operator->() const;
  };


  class matrix: public MatrixView
  {
    // Data member inherited from ConstMatrixView via MatrixView
  public:
    explicit matrix(MatrixBase* MatPtr);
    matrix(const matrix& M);
    matrix& operator=(const matrix&) = default; // to avoid compiler warning
    matrix& operator=(ConstMatrixView Mrhs) { ConstMatrixView::operator=(Mrhs); return *this; }
    // Default dtor is OK.
    MatrixBase* operator->() const;
  };


  //-------------------------------------------------------


  class ConstMatrixViewBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<ConstMatrixViewBase>;

  protected: // Hide default ctor & dtor from public.
    ConstMatrixViewBase() {}
    virtual ~ConstMatrixViewBase() {}

  public:
    typedef std::vector<RingElem> vec; ///???? is this the right thing to do???
    virtual const ring& myRing() const = 0;
    virtual long myNumRows() const = 0;
    virtual long myNumCols() const = 0;
    virtual RingElemAlias myEntry(long i, long j) const = 0;
    virtual bool IamEqual(ConstMatrixView M) const;       ///< default tests all entries
    virtual bool IamSymmetric() const;                    ///< default tests all entries
    virtual bool IamAntiSymmetric() const;                ///< default tests all entries
    virtual bool IamDiagonal() const;                     ///< default tests all entries
    virtual bool IamUpperTriangular() const;              ///< default tests all entries
    virtual bool IamLowerTriangular() const;              ///< default tests all entries
    virtual bool IhaveNegEntry() const;                   ///< default tests all entries
//???    void myCheckRowIndex(long i, const char* where="") const;  // throws if index is too large or neg
///    long myCheckRowIndex(const MachineInt& i, const char* where) const;// throws if index is too large or neg; result is row index
    long myCheckRowIndex(const MachineInt& i, const ErrorContext& ErrCtx) const;// throws if index is too large or neg; result is row index
///????    void myCheckColIndex(long j, const char* where="") const;  // throws if index is too large or neg
///    long myCheckColIndex(const MachineInt& j, const char* where) const;// throws if index is too large or neg; result is col index
    long myCheckColIndex(const MachineInt& j, const ErrorContext& ErrCtx) const;// throws if index is too large or neg; result is col index
    virtual bool myIsZeroRow(long i) const;                ///< tests whether row i is zero -- default "dense" definition
    virtual bool myIsZeroCol(long j) const;                ///< tests whether col j is zero -- default "dense" definition
    virtual void myMulByRow(vec& lhs, const vec& v) const; ///< lhs = v*M -- vector-by-matrix product, default "dense" definition
    virtual void myMulByCol(vec& lhs, const vec& v) const; ///< lhs = M*v -- matrix-by-vector product, default "dense" definition
    virtual void myDet(RingElem& d) const;                 ///< default method uses Gaussian reduction
    virtual long myRank() const;                           ///< default method uses Gaussian reduction
    virtual void myOutputSelf(std::ostream& out) const;    ///< has a default "dense" definition
    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const;///< has a default "dense" definition
  };


  class ConstMatrixBase: public ConstMatrixViewBase
  {
  protected: // Hide default ctor & dtor from public.
    ConstMatrixBase() {}
    virtual ~ConstMatrixBase() {}

  public:
    virtual ConstMatrixBase* myClone() const = 0;
  };

  
  class MatrixViewBase: public ConstMatrixViewBase
  {
  protected: // Hide default ctor & dtor from public.
    MatrixViewBase() {}
    virtual ~MatrixViewBase() {}

  public:
    virtual bool myIsWritable(long i, long j) const = 0; ///< true iff location (i,j) is writable
    virtual RingElemRawPtr myRawEntry(long i, long j) = 0; ///< may be used only if (i,j) is writable
    virtual void mySetEntry(long i, long j, ConstRefRingElem r) = 0;
    virtual void mySetEntry(long i, long j, const MachineInt& n) = 0;
    virtual void mySetEntry(long i, long j, const BigInt& N) = 0;
    virtual void mySetEntry(long i, long j, const BigRat& Q) = 0;
  };


  class MatrixBase: public MatrixViewBase
  {
  protected: // Hide default ctor & dtor from public.
    MatrixBase() {}
    virtual ~MatrixBase() {}

  public:
    virtual MatrixBase* myClone() const = 0;
    virtual MatrixBase* myZeroClone(const ring& R, long NumRows, long NumCols) const = 0;
    virtual void myResize(long NumRows, long NumCols) = 0;
    virtual void myRowMul(long i, ConstRefRingElem c) = 0; ///< row(i) = c*row(i)
    virtual void myColMul(long j, ConstRefRingElem c) = 0; ///< col(j) = c*col(j)
    virtual void myAddRowMul(long i1, long i2, ConstRefRingElem c) = 0; ///< row(i1) += c*row(i2)
    virtual void myAddColMul(long j1, long j2, ConstRefRingElem c) = 0; ///< col(j1) += c*col(j2)
    virtual void mySwapRows(long i1, long i2) = 0;
    virtual void mySwapCols(long j1, long j2) = 0;
  };

  //-------------------------------------------------------


  std::ostream& operator<<(std::ostream& out, const ConstMatrixView& M);
  void swap(matrix& M1, matrix& M2);
  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, ConstRefRingElem r);
  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, const MachineInt& n);
  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, const BigInt& N);
  void SetEntry(MatrixView& M, const MachineInt& i, const MachineInt& j, const BigRat& Q);
  void AddRowMul(matrix& M, long i1, long i2, ConstRefRingElem c);
  void AddColMul(matrix& M, long j1, long j2, ConstRefRingElem c);
  void SwapRows(matrix& M1, const MachineInt& i1, const MachineInt& i2);
  void SwapCols(matrix& M1, const MachineInt& i1, const MachineInt& i2);
  void DeleteRow(matrix& M, const MachineInt& i);
  void DeleteCol(matrix& M, const MachineInt& i);


  /////////////////////////////////////////////////////////////////////////////
  ConstMatrix ZeroMat(const ring& R, const MachineInt& NumRows, const MachineInt& NumCols); // NumCols, NumCols > 0
  ConstMatrix IdentityMat(const ring& R, const MachineInt& dim); // dim > 0

  bool IsZeroMatImpl(const ConstMatrixView& M);
  bool IsIdentityMatImpl(const ConstMatrixView& M);

  /////////////////////////////////////////////////////////////////////////////
  // Inline functions

  inline ConstMatrixView::ConstMatrixView(const ConstMatrixViewBase* MatPtr):
      mySmartPtr(const_cast<ConstMatrixViewBase*>(MatPtr))
  {}


  inline const ConstMatrixViewBase* ConstMatrixView::operator->() const
  { return mySmartPtr.operator->(); }


  inline bool ConstMatrixView::operator==(ConstMatrixView M) const
  {
    if (mySmartPtr == M.mySmartPtr) return true;
    return mySmartPtr->IamEqual(M);
  }


  inline ConstMatrix::ConstMatrix(ConstMatrixBase* MatPtr):
      ConstMatrixView(MatPtr)
  {}

  inline const ConstMatrixBase* ConstMatrix::operator->() const
  { return dynamic_cast<const ConstMatrixBase*>(ConstMatrixView::operator->()); }


  inline MatrixView::MatrixView(MatrixViewBase* MatPtr):
      ConstMatrixView(MatPtr)
  {}


  inline MatrixViewBase* MatrixView::operator->() const
  { return dynamic_cast<MatrixViewBase*>(const_cast<ConstMatrixViewBase*>(ConstMatrixView::operator->())); }



  inline matrix::matrix(MatrixBase* MatPtr):
      MatrixView(MatPtr)
  {}


  inline MatrixBase* matrix::operator->() const
  { return dynamic_cast<MatrixBase*>(MatrixView::operator->()); }


  inline const ring& RingOf(const ConstMatrixView& M)
  {
    return M->myRing();
  }


  inline long NumRows(const ConstMatrixView& M)
  {
    return M->myNumRows();
  }


  inline long NumCols(const ConstMatrixView& M)
  {
    return M->myNumCols();
  }


} // end of namespace CoCoA

#endif
