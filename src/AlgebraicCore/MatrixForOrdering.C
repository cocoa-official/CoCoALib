//   Copyright (c)  2008-2017,2021  John Abbott, Anna M. Bigatti
//   Authors:  2008-2012,2015  John Abbott, Anna M. Bigatti

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

#include "CoCoA/MatrixForOrdering.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    matrix MakeCopyOverZZ(const ConstMatrixView& M, const ErrorContext& ErrCtx)
    {
      const ring& R = RingOf(M);
      if (!(IsZZ(R) || IsQQ(R)))
        CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadRing, "Matrix must be over ZZ or QQ", ErrCtx);
      const long r = NumRows(M);
      const long c = NumCols(M);
      matrix ans = NewDenseMat(RingZZ(), r, c);
      for (long i=0; i < r; ++i)
        for (long j=0; j < c; ++j)
        {
          BigInt Mij;
          if (!IsInteger(Mij, M(i,j)))
            CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadArg, "Matrix must have integer entries", ErrCtx);
          SetEntry(ans,i,j, Mij);
        }
      return ans;
    }

    enum ColCheckFlag { WithoutZeroCol, AllowZeroCols };

    // Check that first non-zero in each col is positive,
    // and there are no null columns if WithoutZeroCol
    bool ColCheck(const ConstMatrixView& M, ColCheckFlag flag)
    {
      const long nrows = NumRows(M);
      if (nrows == 0) return false; // ???? (flag == AllowZeroCols)
      const long ncols = NumCols(M);
      for (long col=0; col < ncols; ++col)
      {
        long row=0;
        for (; row < nrows; ++row)
        {
          const int s = sign(M(row,col));
          if ( s<0 ) return false;
          if ( s>0 ) break;
        }
        if ((row == nrows) && (flag == WithoutZeroCol)) // found zero-column
          return false;
      }
      return true;
    }


    bool IsNonNegGrading(const ConstMatrixView& M)
    {
      if  (!IsZZ(RingOf(M)))
        return IsNonNegGrading(MakeCopyOverZZ(M, CoCoA_ERROR_CONTEXT));

      if (!ColCheck(M, AllowZeroCols)) return false;
      if (NumRows(M)==NumCols(M))    return (!IsZeroDet(M));
      return (rk(M) == NumRows(M));
    }

  } // end of anonymous namespace


  bool IsPositiveGrading(const ConstMatrixView& M)
  {
    if  (!IsZZ(RingOf(M)))
      return IsPositiveGrading(MakeCopyOverZZ(M, CoCoA_ERROR_CONTEXT));

    //    std::cout << "IsPositiveGrading" << std::endl;
    if (!ColCheck(M, WithoutZeroCol)) return false;
    if (NumRows(M)==NumCols(M))    return (!IsZeroDet(M));
    return (rk(M) == NumRows(M));
  }


  bool IsTermOrdering(const ConstMatrixView& M)
  {
    //    std::cout << "IsTermOrdering" << std::endl;
    if (NumCols(M) != NumRows(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    return IsPositiveGrading(M);
  }


  namespace // anonymous
  {

/**
   expects a matrix with entries in an ordered ring: very likely only ZZ (maybe QQ?)
   returns a matrix with positive entries which defines an equivalent ordering
*/

    // Stupid algm, but it "obviously works"
    // ***ASSUMES***  input M has no zero columns!!!
    matrix MakeNonNeg(const ConstMatrixView& M)
    {
      //      std::cout << "--MakeNonNeg" << std::endl;
      CoCoA_ASSERT(IsZZ(RingOf(M)));
      // if  (!IsZZ(RingOf(M)) && !IsQQ(RingOf(M)))
      //   CoCoA_THROW_ERROR2(ERR::BadRing, "matrix must be over RingZZ() or RingQQ()");
      // if ( !IsTermOrdering(M) )  CoCoA_THROW_ERROR1(ERR::ReqTermOrdering);

      matrix PosMat(NewDenseMat(M));
      const long nrows = NumRows(M);
      const long ncols = NumCols(M);
      for (long row=0; row < nrows; ++row)
        for (long col=0; col < ncols; ++col)
          if (PosMat(row,col) < 0)
          {
            long PosRow=0;
            // Loop to find positive entry in col
            while (IsZero(PosMat(PosRow,col)))  ++PosRow; // safe because assumed no zero cols
            CoCoA_ASSERT(PosRow < row);
            const BigInt q = ceil(BigRat(-ConvertTo<BigInt>(PosMat(row,col)), ConvertTo<BigInt>(PosMat(PosRow,col))));
            PosMat->myAddRowMul(row, PosRow, RingElem(RingZZ(),q));
          }
      return PosMat;
    }


    BigInt CommonDenomOfRow(const ConstMatrixView& M, long row)
    {
      CoCoA_ASSERT(IsZZ(RingOf(M)) || IsQQ(RingOf(M)));
      BigInt D(1);
      if (IsZZ(RingOf(M))) return D;
      for (long col=0; col < NumCols(M); ++col)
        D = lcm(D, den(ConvertTo<BigRat>(M(row,col))));
      return D;
    }
  

    // Input matrix may have rational entries, but first GrDim rows must have integer entries.
    // Rescale rows so that output matrix has integer entries.
    matrix ClearDenomByRow(const ConstMatrixView& M, long GrDim)
    {
      CoCoA_ASSERT(IsZZ(RingOf(M)) || IsQQ(RingOf(M)));
      // if  (IsZZ(RingOf(M)))  return NewDenseMat(M);

      const long nrows = NumRows(M);
      const long ncols = NumCols(M);
      matrix IntMat = NewDenseMat(RingZZ(), nrows, ncols);
      for (long row=0; row < nrows; ++row)
      {
        const BigInt D = CommonDenomOfRow(M, row);
        if (row < GrDim && D != 1)  CoCoA_THROW_ERROR2(ERR::BadArg, "row < GrDim && D != 1");
        for (long col=0; col < ncols; ++col)
          SetEntry(IntMat, row,col, ConvertTo<BigInt>(D*M(row,col)));
      }
      return IntMat;
    }


    // Assume full column rank (??? 2022-07-19 ???)
    matrix RemoveRedundantRows(const ConstMatrixView& M)
    {
      CoCoA_ASSERT(IsZZ(RingOf(M)));
      const long nrows = NumRows(M);
      const long ncols = NumCols(M);
      CoCoA_ASSERT(nrows >= ncols);

      if (!IsZeroDet(submat(M,LongRange(0,ncols-1),LongRange(0,ncols-1))))
        return NewDenseMat(submat(M,LongRange(0,ncols-1),LongRange(0,ncols-1)));
      
      matrix MM(NewDenseMat(RingZZ(), ncols, ncols)); // ncols-by-ncols !!
      long row=0;
      for (long r=0; r<nrows; ++r)
      {
        for (long col=0; col<ncols; ++col)
          SetEntry(MM, row, col, M(r, col));
        if (rk(MM) > row) ++row;
        if (row == ncols) return MM;
      }
      CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
      return MM; // just to keep compiler quiet
    }

  } // end of namespace anonymous

  
  /**************************************************************************/
  // ConstMatrix for lex ordering

  ConstMatrix LexMat(const MachineInt& n)
  {
    if (IsNegative(n) || IsZero(n))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    return IdentityMat(RingZZ(), AsSignedLong(n));
  }


  /***************************************************************************/
  // ConstMatrix for xel ordering incl. auxiliary class

  class XelMatImpl: public ConstMatrixBase
  {
  private:
    friend ConstMatrix XelMat(const MachineInt& dim); // pseudo-ctor, uses RingQQ
    XelMatImpl(const ring& R, long dim);
    // default dtor is fine
  public: // disable default copy ctor and assignment
    XelMatImpl(const XelMatImpl&) = delete;
    XelMatImpl& operator=(const XelMatImpl&) = delete;

  public: // functions every matrix type must implement
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
//NY default    virtual bool IamEqual(const ConstMatrixView& M) const override;
    bool IamSymmetric() const override  {return true;}
    bool IamAntiSymmetric() const override  {return (myDim==0);}
    bool IamDiagonal() const override  {return (myDim<=1);}
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    ConstMatrixBase* myClone() const override;

  private: // data members
    const ring myR;
    long myDim;
  };


  XelMatImpl::XelMatImpl(const ring& R, long dim):
      ConstMatrixBase(),
      myR(R),
      myDim(dim)
  {}


  const ring& XelMatImpl::myRing() const
  {
    return myR;
  }


  long XelMatImpl::myNumRows() const
  {
    return myDim;
  }


  long XelMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias XelMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i+j == myDim-1) return one(myR);
    return zero(myR);
  }


  //BUG: should use FreeModule elems!
  void XelMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    // ASSUMES NO ALIASING
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    const long n = len(v);
    for (long i=0; i < n; ++i)
      lhs[i] = v[n-i-1];
  }


  void XelMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    const long n = len(v);
    for (long i=0; i < n; ++i)
      lhs[i] = v[n-i-1];
  }

  // bool XelMatImpl::IamEqual(const ConstMatrixView& M) const
  // {
  //   if (myRing() != RingOf(M)) return false;
  //   if (myNumRows() != NumRows(M)) return false;
  //   if (myNumCols() != NumCols(M)) return false;
  //   if (IsZeroMatImpl(M)) return NumCols(M) == 0;
  //   if (IsXelMatImpl(M)) return true;
  //   if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
  //   {
  //     //      std::cout << "XelMatImpl::IamEqual - diag" << std::endl;
  //     for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
  //     return true;
  //   }
  //   return ConstMatrixViewBase::IamEqual(M);
  // }


  bool XelMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= i && i < myNumRows());
    return false;
  }


  bool XelMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return false;
  }


  void XelMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    if (IsEven(myDim/2))
      d = 1;
    else
      d = -1;
  }


  long XelMatImpl::myRank() const
  {
    return myNumRows();
  }


  ConstMatrixBase* XelMatImpl::myClone() const
  {
    return new XelMatImpl(myR, myDim);
  }


  ConstMatrix XelMat(const MachineInt& n)
  {
    if (IsNegative(n) || IsZero(n))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    return ConstMatrix(new XelMatImpl(RingZZ(), AsSignedLong(n)));
  }



  /***************************************************************************/
  // ConstMatrix for RevLex ordering incl. auxiliary class

  class RevLexMatImpl: public ConstMatrixBase
  {
  private:
    friend ConstMatrix RevLexMat(const MachineInt& dim); // pseudo-ctor, uses RingQQ
    RevLexMatImpl(const ring& R, long dim);
    // default dtor is fine
  public: // disable default copy ctor and assignment
    RevLexMatImpl(const RevLexMatImpl&) = delete;
    RevLexMatImpl& operator=(const RevLexMatImpl&) = delete;

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
//NY default    bool IamEqual(const ConstMatrixView& M) const override;
    bool IamSymmetric() const override  {return myNumRows()<=1;}
    bool IamAntiSymmetric() const override  {return myNumRows()==0;}
    bool IamDiagonal() const override  {return myNumRows()<=1;}
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    ConstMatrixBase* myClone() const override;

  private: // data members
    const ring myR;
    long myDim;
    const RingElem myMinusOne;
  };



  RevLexMatImpl::RevLexMatImpl(const ring& R, long dim):
      ConstMatrixBase(),
      myR(R),
      myDim(dim),
      myMinusOne(myR,-1)
  {}


  const ring& RevLexMatImpl::myRing() const
  {
    return myR;
  }


  long RevLexMatImpl::myNumRows() const
  {
    return myDim;
  }


  long RevLexMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias RevLexMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i+j == myDim-1) return myMinusOne;
    return zero(myR);
  }


  //BUG: should use FreeModule elems!
  void RevLexMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    // ASSUMES NO ALIASING
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    const long n = len(v);
    for (long i=0; i < n; ++i)
      lhs[i] = -v[n-i-1];
  }


  void RevLexMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    const long n = len(v);
    for (long i=0; i < n; ++i)
      lhs[i] = -v[n-i-1];
  }

  // bool RevLexMatImpl::IamEqual(const ConstMatrixView& M) const
  // {
  //   if (myRing() != RingOf(M)) return false;
  //   if (myNumRows() != NumRows(M)) return false;
  //   if (myNumCols() != NumCols(M)) return false;
  //   if (IsZeroMatImpl(M)) return NumCols(M) == 0;
  //   if (IsRevLexMatImpl(M)) return true;
  //   if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
  //   {
  //     //      std::cout << "RevLexMatImpl::IamEqual - diag" << std::endl;
  //     for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
  //     return true;
  //   }
  //   return ConstMatrixViewBase::IamEqual(M);
  // }


  bool RevLexMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= i && i < myNumRows());
    return false;
  }


  bool RevLexMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return false;
  }


  void RevLexMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    const int DimMod4 = myDim%4;
    if (DimMod4 < 3)
      d = 1;
    else
      d = -1;
  }


  long RevLexMatImpl::myRank() const
  {
    return myNumRows();
  }


  ConstMatrixBase* RevLexMatImpl::myClone() const
  {
    return new RevLexMatImpl(myR, myDim);
  }


  ConstMatrix RevLexMat(const MachineInt& n)
  {
    if (IsNegative(n) || IsZero(n))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    return ConstMatrix(new RevLexMatImpl(RingZZ(), AsSignedLong(n)));
  }



  /***************************************************************************/
  // ConstMatrix for StdDegLex ordering incl. auxiliary class

  class StdDegLexMatImpl: public ConstMatrixBase
  {
  private:
    friend ConstMatrix StdDegLexMat(const MachineInt& dim); // pseudo-ctor, uses RingQQ
    StdDegLexMatImpl(const ring& R, long dim);
    // default dtor is fine
  public: // disable default copy ctor and assignment
    StdDegLexMatImpl(const StdDegLexMatImpl&) = delete;
    StdDegLexMatImpl& operator=(const StdDegLexMatImpl&) = delete;

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
//NY default    bool IamEqual(const ConstMatrixView& M) const override;
    bool IamSymmetric() const override  {return myNumRows()<=1;}
    bool IamAntiSymmetric() const override  {return myNumRows()==0;}
    bool IamDiagonal() const override  {return myNumRows()<=1;}
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    ConstMatrixBase* myClone() const override;

  private: // data members
    const ring myR;
    long myDim;
  };



  StdDegLexMatImpl::StdDegLexMatImpl(const ring& R, long dim):
      ConstMatrixBase(),
      myR(R),
      myDim(dim)
  {}


  const ring& StdDegLexMatImpl::myRing() const
  {
    return myR;
  }


  long StdDegLexMatImpl::myNumRows() const
  {
    return myDim;
  }


  long StdDegLexMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias StdDegLexMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
    if (i == 0 || i == j+1) return one(myR);
    return zero(myR);
  }


  //BUG: should use FreeModule elems!
  void StdDegLexMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    // ASSUMES NO ALIASING
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    if (myDim == 0) return;
    const long n = len(v);
    lhs[n-1] = v[n-1];
    for (long i=n-2; i >= 0; --i)
    {
      lhs[i] = v[0]+v[i+1];
    }
  }


  void StdDegLexMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    const long n = len(v);
    RingElem sum = v[0];
    for (long i=n-1; i > 0; --i)
    {
      sum += v[i];
      lhs[i] = v[i-1];
    }
    lhs[0] = sum;
  }

  // bool StdDegLexMatImpl::IamEqual(const ConstMatrixView& M) const
  // {
  //   if (myRing() != RingOf(M)) return false;
  //   if (myNumRows() != NumRows(M)) return false;
  //   if (myNumCols() != NumCols(M)) return false;
  //   if (IsZeroMatImpl(M)) return NumCols(M) == 0;
  //   if (IsStdDegLexMatImpl(M)) return true;
  //   if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
  //   {
  //     //      std::cout << "StdDegLexMatImpl::IamEqual - diag" << std::endl;
  //     for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
  //     return true;
  //   }
  //   return ConstMatrixViewBase::IamEqual(M);
  // }


  bool StdDegLexMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= i && i < myNumRows());
    return false;
  }


  bool StdDegLexMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return false;
  }


  void StdDegLexMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    d = 1;
  }


  long StdDegLexMatImpl::myRank() const
  {
    return myNumRows();
  }


  ConstMatrixBase* StdDegLexMatImpl::myClone() const
  {
    return new StdDegLexMatImpl(myR, myDim);
  }


  ConstMatrix StdDegLexMat(const MachineInt& n)
  {
    if (IsNegative(n) || IsZero(n))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    return ConstMatrix(new StdDegLexMatImpl(RingZZ(), AsSignedLong(n)));
  }



  /***************************************************************************/
  // ConstMatrix for StdDegRevLex ordering incl. auxiliary class

  class StdDegRevLexMatImpl: public ConstMatrixBase
  {
  private:
    friend ConstMatrix StdDegRevLexMat(const MachineInt& dim); // pseudo-ctor
    StdDegRevLexMatImpl(const ring& R, long dim);
    // default dtor is fine
  public: // disable default copy ctor and assignment
    StdDegRevLexMatImpl(const StdDegRevLexMatImpl&) = delete;
    StdDegRevLexMatImpl& operator=(const StdDegRevLexMatImpl&) = delete;

  public:
    typedef std::vector<RingElem> vec;
    const ring& myRing() const override;
    long myNumRows() const override;
    long myNumCols() const override;
    RingElemAlias myEntry(long i, long j) const override;
    void myMulByRow(vec& lhs, const vec& v) const override;
    void myMulByCol(vec& lhs, const vec& v) const override;
//NY default    bool IamEqual(const ConstMatrixView& M) const override;
    bool IamSymmetric() const override  {return myNumRows()<=1;}
    bool IamAntiSymmetric() const override  {return myNumRows()==0;}
    bool IamDiagonal() const override  {return myNumRows()<=1;}
    bool myIsZeroRow(long i) const override;
    bool myIsZeroCol(long j) const override;
    void myDet(RingElem& d) const override;
    long myRank() const override;

    ConstMatrixBase* myClone() const override;

  private: // data members
    const ring myR;
    long myDim;
    const RingElem myMinusOne;
  };



  StdDegRevLexMatImpl::StdDegRevLexMatImpl(const ring& R, long dim):
      ConstMatrixBase(),
      myR(R),
      myDim(dim),
      myMinusOne(R,-1)
  { CoCoA_ASSERT(dim > 0); }


  const ring& StdDegRevLexMatImpl::myRing() const
  {
    return myR;
  }


  long StdDegRevLexMatImpl::myNumRows() const
  {
    return myDim;
  }


  long StdDegRevLexMatImpl::myNumCols() const
  {
    return myDim;
  }


  RingElemAlias StdDegRevLexMatImpl::myEntry(long i, long j) const
  {
    CoCoA_ASSERT(i >= 0 && j >= 0);
    CoCoA_ASSERT(i < myNumRows() && j < myNumCols());
//  Uncomment next 2 lines for the non-neg matrix for StdDegRevLex
//    if (i < myDim-j) return one(myR);
//    return zero(myR);
    if (i == 0) return one(myR);
    if (i+j == myDim) return myMinusOne;
    return zero(myR);
  }


  //BUG: should use FreeModule elems!
  void StdDegRevLexMatImpl::myMulByRow(vec& lhs, const vec& v) const
  {
    // ASSUMES NO ALIASING
    CoCoA_ASSERT(len(lhs) == myNumCols() && len(v) == myNumRows());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    if (myDim == 0) return;
    const long n = len(v);
    RingElem sum = v[0];
    for (long i=1; i < n; ++i)
    {
      lhs[myDim-i] = sum;
      sum += v[i];
    }
    lhs[0] = sum;
  }


  void StdDegRevLexMatImpl::myMulByCol(vec& lhs, const vec& v) const
  {
    CoCoA_ASSERT(len(lhs) == myNumRows() && len(v) == myNumCols());
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(v[0]));
    CoCoA_ASSERT(myNumRows() == 0 || myRing() == owner(lhs[0]));
//BUG???    CoCoA_ASSERT(myRing() == RingOf(v));
    if (myDim == 0) return;
    const long n = len(v);
    RingElem sum = v[0];
    for (long i=1; i < n; ++i)
    {
      lhs[myDim-i] = sum;
      sum += v[i];
    }
    lhs[0] = sum;
  }

  // bool StdDegRevLexMatImpl::IamEqual(const ConstMatrixView& M) const
  // {
  //   if (myRing() != RingOf(M)) return false;
  //   if (myNumRows() != NumRows(M)) return false;
  //   if (myNumCols() != NumCols(M)) return false;
  //   if (IsZeroMatImpl(M)) return NumCols(M) == 0;
  //   if (IsStdDegRevLexMatImpl(M)) return true;
  //   if (IsDiagMatImpl(M) || IsConstDiagMatImpl(M))
  //   {
  //     //      std::cout << "StdDegRevLexMatImpl::IamEqual - diag" << std::endl;
  //     for (long i=0; i < myNumRows(); ++i) if (!IsOne(M(i,i))) return false;
  //     return true;
  //   }
  //   return ConstMatrixViewBase::IamEqual(M);
  // }


  bool StdDegRevLexMatImpl::myIsZeroRow(long i) const
  {
    (void)(i); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= i && i < myNumRows());
    return false;
  }


  bool StdDegRevLexMatImpl::myIsZeroCol(long j) const
  {
    (void)(j); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(0 <= j && j < myNumCols());
    return false;
  }


  void StdDegRevLexMatImpl::myDet(RingElem& d) const
  {
    CoCoA_ASSERT(owner(d) == myRing());
    const int DimMod4 = myDim%4;
    if (DimMod4 < 3)
      d = 1;
    else
      d = -1;
  }


  long StdDegRevLexMatImpl::myRank() const
  {
    return myNumRows();
  }


  ConstMatrixBase* StdDegRevLexMatImpl::myClone() const
  {
    return new StdDegRevLexMatImpl(myR, myDim);
  }


  ConstMatrix StdDegRevLexMat(const MachineInt& n)
  {
    if (IsNegative(n) || IsZero(n))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    return ConstMatrix(new StdDegRevLexMatImpl(RingZZ(), AsSignedLong(n)));
  }



  /**********************************************************************/
  // Sundry functions

  matrix MakeTermOrdMat(ConstMatrixView M)
  {
    return MakeTermOrdMat(M, 0);
  }
  
  matrix MakeTermOrdMat(ConstMatrixView M, const MachineInt& GradingDim)
  {
    //    std::cout << "----MakeTermOrd" << std::endl;
    if (IsNegative(GradingDim))  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const long GrDim = AsSignedLong(GradingDim);
    if (GrDim > NumRows(M))  CoCoA_THROW_ERROR2(ERR::BadIndex, "row index");
    if (!IsZZ(RingOf(M)))
      return MakeTermOrdMat(MakeCopyOverZZ(ClearDenomByRow(M,GrDim),
                                           CoCoA_ERROR_CONTEXT));
    
    matrix TORow = NewDenseMat(RingZZ(), 1,NumCols(M));
    for (long c=0; c<NumCols(M); ++c)
      if (IsZeroCol(M,c))  SetEntry(TORow, 0,c, 1);
    if (!IsZeroRow(TORow,0))
      M = ConcatVer(M, TORow);
    matrix NonNegMat = MakeNonNeg(M);
    if (!ColCheck(M, WithoutZeroCol))
      CoCoA_THROW_ERROR2(ERR::BadArg, "Topmost non-zero entry in each col must be positive");
    if (NumRows(NonNegMat)==NumCols(NonNegMat) && !IsZeroDet(NonNegMat))
      return NonNegMat;
    return RemoveRedundantRows(ConcatVer(NonNegMat, RevLexMat(NumCols(M))));
  }


  //--- matrices for elimination -----------------------------

  namespace // anonymous namespace for file local auxiliary funcs and defs
  {
    
    matrix ElimRow(const std::vector<long>& IndetsToElim, long NumIndets)
    {
      matrix M = NewDenseMat(RingZZ(), 1, NumIndets);
      long s = len(IndetsToElim);
      for (long i=0; i < s; ++i)
      {
        if (IndetsToElim[i]<0 || IndetsToElim[i]>=NumIndets)
          CoCoA_THROW_ERROR1(ERR::BadIndex);
        SetEntry(M,  0, IndetsToElim[i],  1);
      }
      return M;
    }

  } // end of namespace anonymous
  
  
  matrix ElimMat(const std::vector<long>& IndetsToElim, const MachineInt& NumIndets)
  {
    if (IsNegative(NumIndets) || IsZero(NumIndets) || !IsSignedLong(NumIndets))
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    const long n = AsSignedLong(NumIndets);
    return MakeTermOrdMat(ConcatVer(ElimRow(IndetsToElim, n),
                                 RowMat(vector<RingElem>(n, one(RingZZ())))));
  }


  matrix ElimMat(const std::vector<long>& IndetsToElim,
                 const ConstMatrixView& GradM)
  {
    if (NumRows(GradM)==0)
      return ElimMat(IndetsToElim, NumCols(GradM)); // with row of 1's
    if (!IsNonNegGrading(GradM))  CoCoA_THROW_ERROR1(ERR::ReqNonNegativeGrading);
    if (!IsZZ(RingOf(GradM)))
      return ElimMat(IndetsToElim, MakeCopyOverZZ(GradM, CoCoA_ERROR_CONTEXT));
    return MakeTermOrdMat(ConcatVer(ElimRow(IndetsToElim, NumCols(GradM)),
                                 GradM));
  }


  matrix ElimHomogMat(const std::vector<long>& IndetsToElim,
                      const ConstMatrixView& GradM)
  {
    if (NumRows(GradM)==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (!IsNonNegGrading(GradM))  CoCoA_THROW_ERROR1(ERR::ReqNonNegativeGrading);
    if (!IsZZ(RingOf(GradM)))
      return ElimHomogMat(IndetsToElim,
                          MakeCopyOverZZ(GradM, CoCoA_ERROR_CONTEXT));
    return MakeTermOrdMat(ConcatVer(GradM,
                                    ElimRow(IndetsToElim, NumCols(GradM))));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixForOrdering.C,v 1.58 2024/08/28 15:23:06 bigatti Exp $
// $Log: MatrixForOrdering.C,v $
// Revision 1.58  2024/08/28 15:23:06  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2 and CoCoA_THROW_ERROR_WITH_CONTEXT2/3
//
// Revision 1.57  2024/07/15 17:03:56  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.56  2024/07/03 15:40:41  bigatti
// Summary: Error codes: changed BadIndetIndex into BadIndex
//
// Revision 1.55  2024/07/02 15:31:12  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.54  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.53  2024/03/15 18:56:03  abbott
// Summary: Improved a comment
//
// Revision 1.52  2024/03/12 14:25:03  bigatti
// Summary: now using IsZeroDet(M) instead of rk/det
//
// Revision 1.51  2023/10/05 17:57:33  abbott
// Summary: Improved error mesg (redmine 1761)
//
// Revision 1.50  2023/08/06 10:14:41  abbott
// Summary: Added some spaces
//
// Revision 1.49  2022/08/05 20:43:53  abbott
// Summary: Took NumRows/NumCols out of for-loop exit conds
//
// Revision 1.48  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.47  2022/02/04 21:41:06  abbott
// Summary: Changed name MakeTermOrd to MakeTermOrdMat (redine 854)
//
// Revision 1.46  2021/10/30 16:48:27  abbott
// Summary: Used keywords override and delete (redmine 1625 & 1627)
//
// Revision 1.45  2021/10/04 08:53:51  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.44  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.43  2020/02/13 10:14:00  abbott
// Summary: Added const somewhere
//
// Revision 1.42  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.41  2019/09/25 14:25:05  bigatti
// -- just some resorting of the code for IsPositiveGrading
// -- IsNonNegGrading is now hidden in anonymous namespace
//
// Revision 1.40  2018/08/28 12:37:39  abbott
// Summary: Changed ERR::BadArg into ERR::NotSquareMatrix
//
// Revision 1.39  2018/07/04 13:09:27  bigatti
// -- minor optimization for MakeTermOrd
//
// Revision 1.38  2018/06/15 08:46:26  abbott
// Summary: Added IsNonNegGrading; ElimMat, ElimHomogMat req grading to be non-neg (and mat over ZZ or QQ)
//
// Revision 1.37  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.36  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.35  2018/05/17 15:35:34  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.34  2017/09/06 11:56:28  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.33  2017/05/17 15:56:36  bigatti
// -- added check for det of square matrix (MakeTermOrd)
//
// Revision 1.32  2017/01/25 13:00:48  abbott
// Summary: Removed cruft
//
// Revision 1.31  2016/09/22 15:33:37  bigatti
// -- renamed HomogElimMat into ElimHomogMat
// -- improved readability for ElimHomogMat/ElimMat (removed auxiliary functions)
//
// Revision 1.30  2015/12/11 13:03:20  bigatti
// -- added IsUpperTriangular, IhaveNegEntry
// -- removed useless checks
//
// Revision 1.29  2015/12/10 09:18:45  abbott
// Summary: Added new "hidden" fn ColCheck
//
// Revision 1.28  2015/12/08 14:05:11  abbott
// Summary: Renamed NewMatCompleteOrd to MakeTermOrd
//
// Revision 1.27  2015/12/01 16:53:15  abbott
// Summary: Added new MatCompleteOrd with 2nd arg
//
// Revision 1.26  2015/12/01 15:57:36  abbott
// Summary: Rewriting the fns to normalize/simplify an partial order matrix
//
// Revision 1.25  2015/12/01 13:34:44  abbott
// Summary: Changed arg order in ElimMat and HomogElimMat; doc is out-of-date!!
//
// Revision 1.24  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.23  2015/07/28 13:27:53  bigatti
// -- added check for GrDim in IsPositiveGrading
//
// Revision 1.22  2015/04/13 15:35:24  abbott
// Summary: Changed "rank" --> "rk"
// Author: JAA
//
// Revision 1.21  2014/07/30 14:06:06  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.20  2014/04/17 13:38:41  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.19  2014/04/11 15:44:27  abbott
// Summary: Renamed MatrixArith to MatrixOperations (in includes)
// Author: JAA
//
// Revision 1.18  2014/04/08 12:56:12  bigatti
// -- cambiata ElimMat per non usare FilledMat
//
// Revision 1.17  2013/02/14 17:34:35  bigatti
// -- cleaned up code for elimination matrices
//
// Revision 1.16  2012/11/21 09:46:33  bigatti
// -- NewMatCompleteOrd(M) now returns a term-ordering (if M suitable)
//
// Revision 1.15  2012/03/30 17:28:09  bigatti
// -- added NewIntegerOrdMat
// -- accepting and returning matrices over QQ
//
// Revision 1.14  2012/02/10 10:28:08  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.13  2012/02/08 17:20:43  bigatti
// -- changed: Z,Q -> ZZ,QQ
// -- code reorganization
//
// Revision 1.12  2011/05/26 11:58:05  bigatti
// -- added IsPositiveGrading with one arg
//
// Revision 1.11  2011/04/26 10:11:06  bigatti
// -- added NewMatCompleteOrd
//
// Revision 1.10  2011/03/23 17:29:54  bigatti
// -- added NewDenseStdDegLexMat
//
// Revision 1.9  2011/03/21 07:50:51  bigatti
// -- added NewDenseMatXel, NewDenseMatStdDegRevLex
//
// Revision 1.8  2011/03/10 11:26:26  bigatti
// -- using len(v) instead of v.size()
//
// Revision 1.7  2011/03/09 09:08:25  bigatti
// -- removing signed/unsigned warnings
// -- added samity checks for indets to elim
//
// Revision 1.6  2011/03/08 16:10:16  abbott
// Changed size_t into long.
//
// Revision 1.5  2011/02/16 14:58:40  bigatti
// -- fixed typo
//
// Revision 1.4  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.3  2009/09/22 13:35:55  bigatti
// -- following coding conventions in function names Matrix --> Mat
// -- forced all matrices to be over RingZZ
//
// Revision 1.2  2008/05/30 12:44:14  abbott
// Moved "ordering matrices" into their ownn special file.
//
// Revision 1.1  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
//
