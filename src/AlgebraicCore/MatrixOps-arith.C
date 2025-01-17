//   Copyright (c)  2005,2008,2016  John Abbott and Anna M. Bigatti

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

#include "CoCoA/MatrixOps.H"

#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/ring-AutomaticConversion.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"

#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{

  // Naive dense matrix multiplication
  // Currently just creates a DenseMat to contain the answer.
  // BUG: must make this behave "intelligently" when multiplying two sparse matrices
  // (what does "sparse" mean?  e.g. consider  transpose(SparseMat)*SparseMat
  //  or even DiagMat(...)*AnyMat, or even ZeroMat*AnyMat,... lots of cases!!!??? BUG
  matrix operator*(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const char* const FnName = "Mat*Mat";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,FnName);
    if (NumCols(Mleft) != NumRows(Mright))
      CoCoA_THROW_ERROR(ERR::BadMatrixSize, FnName);
    const ring& Rleft = RingOf(Mleft);
    const ring& Rright = RingOf(Mright);
    if (Rleft != Rright)
    {
      const RingHom promote = AutomaticConversionHom(Rleft,Rright,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == Rleft)
        return Mleft * promote(Mright);
      return promote(Mleft) * Mright;
    }

    const long N = NumCols(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mright);
    matrix ans = NewDenseMat(Rleft, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
      {
        CheckForInterrupt("Matrix multiplication");
        RingElem tmp(Rleft);
        for (long k=0; k < N; ++k)
          tmp += Mleft(i,k)*Mright(k,j);
        SetEntry(ans, i, j, tmp);
      }
    return ans;
  }


  matrix operator+(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const char* const FnName = "Mat+Mat";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,FnName);
    const ring& Rleft = RingOf(Mleft);
    const ring& Rright = RingOf(Mright);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mleft);
    if (NumRows(Mright) != Nrows) CoCoA_THROW_ERROR(ERR::BadMatrixSize, FnName);
    if (NumCols(Mright) != Ncols) CoCoA_THROW_ERROR(ERR::BadMatrixSize, FnName);
    if (Rleft != Rright)
    {
      const RingHom promote = AutomaticConversionHom(Rleft,Rright,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == Rleft)
        return Mleft + promote(Mright);
      return promote(Mleft) + Mright;
    }

    matrix ans = NewDenseMat(Rleft, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, Mleft(i,j)+Mright(i,j));
    return ans;
  }


  matrix operator-(ConstMatrixView Mleft, ConstMatrixView Mright)
  {
    const char* const FnName = "Mat-Mat";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,FnName);
    const ring& Rleft = RingOf(Mleft);
    const ring& Rright = RingOf(Mleft);
    const long Nrows = NumRows(Mleft);
    const long Ncols = NumCols(Mleft);
    if (NumRows(Mright) != Nrows) CoCoA_THROW_ERROR(ERR::BadMatrixSize, FnName);
    if (NumCols(Mright) != Ncols) CoCoA_THROW_ERROR(ERR::BadMatrixSize, FnName);
    if (Rleft != Rright)
    {
      const RingHom promote = AutomaticConversionHom(Rleft,Rright,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == Rleft)
        return Mleft + promote(Mright);
      return promote(Mleft) + Mright;
    }

    matrix ans = NewDenseMat(Rleft, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, Mleft(i,j)-Mright(i,j));
    return ans;
  }


  // // Innermost loop is ugly but rather faster than a "clean" implementation.
  // void mul(matrix& lhs, ConstMatrixView M1, ConstMatrixView M2)
  // {
  //   const ring& R = RingOf(lhs);
  //   if (NumCols(M1) != NumRows(M2))
  //     CoCoA_THROW_ERROR(ERR::BadMatrixSize, "mul(Mat1, Mat2)");
  //   if ( NumRows(lhs) != NumRows(M1) || NumCols(lhs) != NumCols(M2) )
  //     CoCoA_THROW_ERROR(ERR::BadMatrixSize, "mul(Mat1, Mat2)");
  //   if (RingOf(M1) != R || RingOf(M2) != R)
  //     CoCoA_THROW_ERROR(ERR::MixedRings, "mul(Mat1, Mat2)");
  //   // Use of the temporary ans avoids aliasing problems and makes the code exception safe.
  //   matrix ans(lhs->myZeroClone(RingOf(lhs), NumRows(M1), NumCols(M2)));
  //   RingElem sum(R), prod(R);
  //   for (long i=0; i < NumRows(M1); ++i)
  //     for (long j=0; j < NumCols(M2); ++j)
  //     {
  //       sum = 0;
  //       for (long k=0; k < NumCols(M1); ++k)
  //       {
  //         // Next 2 lines just do:  sum += M1(i,k) * M2(k,j);
  //         R->myMul(raw(prod), raw(M1(i,k)), raw(M2(k,j)));
  //         R->myAdd(raw(sum), raw(sum), raw(prod));
  //       }
  //       SetEntry(ans, i, j, sum);
  //     }

  //   // The answer is in ans; now swap it with the entries of lhs.
  //   swap(lhs, ans);
  // }

  matrix operator*(ConstRefRingElem x, ConstMatrixView M)
  {
    static const char* const FnName = "RingElem*Mat";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings, FnName);
    const ring& Rx = owner(x);
    const ring& R = RingOf(M);
    if (Rx != R)
    {
      const RingHom promote = AutomaticConversionHom(Rx,R,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == Rx)
        CoCoA_THROW_ERROR("MixedRings: no automatic promotion mapping for matrices", FnName);
/////      return x*promote(M);
      return promote(x) * M;
    }
    // Case: x and M are in same ring
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);

    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, x*M(i,j));
    return ans;
  }
  
  matrix operator*(const BigRat& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }

  matrix operator*(const BigInt& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }

  matrix operator*(const MachineInt& x, ConstMatrixView M)
  { return RingElem(RingOf(M),x) * M; }


  // Separate impl in case ring is not commutative
  matrix operator*(ConstMatrixView M, ConstRefRingElem x)
  {
    static const char* const FnName = "Mat*RingElem";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings, FnName);
    const ring& Rx = owner(x);
    const ring& R = RingOf(M);
    if (Rx != R)
    {
      const RingHom promote = AutomaticConversionHom(Rx,R,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == Rx)
        CoCoA_THROW_ERROR("MixedRings: no automatic promotion mapping for matrices", FnName);
///        return promote(M) * x;
      return M * promote(x);
    }
    // Case: x and M are in same ring
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, M(i,j)*x);
    return ans;
  }

  matrix operator*(ConstMatrixView M, const BigRat& x)
  { return M * RingElem(RingOf(M),x); }

  matrix operator*(ConstMatrixView M, const BigInt& x)
  { return M * RingElem(RingOf(M),x); }

  matrix operator*(ConstMatrixView M, const MachineInt& x)
  { return M * RingElem(RingOf(M),x); }


  matrix operator-(const ConstMatrixView& M)
  {
    return RingElem(RingOf(M),-1)*M;
  }


  matrix operator/(ConstMatrixView M, ConstRefRingElem x)
  {
    static const char* const FnName = "Mat/RingElem";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,FnName);
    if (IsZeroDivisor(x)) CoCoA_THROW_ERROR(ERR::DivByZero, FnName);
    const ring& Rx = owner(x);
    const ring& R = RingOf(M);
    if (Rx != R)
    {
      const RingHom promote = AutomaticConversionHom(Rx,R,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == Rx)
        CoCoA_THROW_ERROR("MixedRings: no automatic promotion mapping for matrices", FnName);
///        return promote(M) / x;
      return M / promote(x);
    }
    // Case: x and M are in same ring
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix ans = NewDenseMat(R, Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i, j, M(i,j)/x);
    return ans;
  }

  matrix operator/(ConstMatrixView M, const BigRat& x)
  { return M / RingElem(RingOf(M),x); }

  matrix operator/(ConstMatrixView M, const BigInt& x)
  { return M / RingElem(RingOf(M),x); }

  matrix operator/(ConstMatrixView M, const MachineInt& x)
  { return M / RingElem(RingOf(M),x); }


  namespace // anonymous
  {
    matrix PowerOfDiagMat(ConstMatrixView M, long n)
    {
      CoCoA_ASSERT(IsDiagonal(M));
//???      CoCoA_ASSERT(n > 0);
      const long sz = NumRows(M); // same as NumCols(M);
      matrix ans = NewDenseMat(IdentityMat(RingOf(M), sz));
      for (int i=0; i < sz; ++i)
        SetEntry(ans,i,i, power(M(i,i),n));
      return ans;
    }
    
  } // end of anon namespace

  matrix power(ConstMatrixView M, long n)
  {
    const ring R = RingOf(M);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (Nrows != Ncols) CoCoA_THROW_ERROR(ERR::BadMatrixSize, "power(M,n)");
    if (n == numeric_limits<long>::min()) CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(M,n)");
    if (n < 0) return power(inverse(M), -n); // cannot overflow because we have excluded n == MinLong
    if (n == 0) return NewDenseMat(IdentityMat(R,Nrows));
    if (IsDiagonal(M)) return PowerOfDiagMat(M,n);

    // An iterative implementation of binary powering.
    long bit = 1; while (bit <= n/2) bit <<= 1;
    matrix ans = NewDenseMat(M);
    while (bit > 1)
    {
      ans = ans*ans;/////mul(ans,ans,ans);//ans *= ans;
      bit >>= 1;
      if (n&bit) ans = ans*M;/////mul(ans,ans,M);//ans *= M;
    }
    return ans;
  }


  matrix power(ConstMatrixView M, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(M,N)");
    return power(M, n);
  }


  

} // end of namespace CoCoA
