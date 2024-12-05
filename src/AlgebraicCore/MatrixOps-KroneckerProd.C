//   Copyright (c)  2020  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/error.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/ring-AutomaticConversion.H"
#include "CoCoA/RingHom.H"

#include <limits>
using std::numeric_limits;

namespace CoCoA
{

  //  Simple "dense" impl
  matrix KroneckerProd(ConstMatrixView M1, ConstMatrixView M2)
  {
    const char* const FnName = "KroneckerProd";
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,FnName);
    const ring& R1 = RingOf(M1);
    const ring& R2 = RingOf(M2);
    if (R1 != R2)
    {
      const RingHom promote = AutomaticConversionHom(R1,R2,ErrMixed); // throws ErrMixed if auto-conv not possible
      if (codomain(promote) == R1)
        return KroneckerProd(M1, promote(M2));
      return KroneckerProd(promote(M1), M2);
    }
    // Here M1 and M2 are over the same ring.
    const long r1 = NumRows(M1);  // int or long???
    const long c1 = NumCols(M1);  //
    const long r2 = NumRows(M2);  //
    const long c2 = NumCols(M2);  //
    if (c1 > 1 && c2 > 1 && c1 > numeric_limits<long>::max() / c2) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "col index"); // ERR::TooBig???
    if (r1 > 1 && r2 > 1 && r1 > numeric_limits<long>::max() / r2) //avoid overflow  
      CoCoA_THROW_ERROR2(ERR::BadIndex, "row index"); // ERR::TooBig???

    const long r = r1*r2;
    const long c = c1*c2;
    
    matrix M = NewDenseMat(R1, r,c);

    for (long i1=0; i1 < r1; ++i1)
      for (long i2=0; i2 < r2; ++i2)
        for (long j1=0; j1 < c1; ++j1)
          for (long j2=0; j2 < c2; ++j2)
            SetEntry(M,i1*r2+i2, j1*c2+j2, M1(i1,j1)*M2(i2,j2));

    return M;
  }
  
  // ORIG DEFN; better local var names?
  // //---------  tensor matrix

  // matrix TensorMat(ConstMatrixView A, ConstMatrixView B)
  // {
  //   if (RingOf(A) != RingOf(B))
  //     CoCoA_THROW_ERROR(ERR::MixedRings, "TensorMat(A,B)");
  //   if (NumCols(A) > numeric_limits<long>::max() /NumCols(B)) //avoid overflow  
  //     CoCoA_THROW_ERROR(ERR::BadColIndex, "TensorMat(A,B)");
  //   if (NumRows(A) > numeric_limits<long>::max() /NumRows(B)) //avoid overflow  
  //     CoCoA_THROW_ERROR(ERR::BadRowIndex, "TensorMat(A,B)");
  //   matrix ans = NewDenseMat(RingOf(A), NumRows(A)*NumRows(B), NumCols(A)*NumCols(B));
  //   for (long iA=0; iA < NumRows(A); ++iA)
  //     for (long jA=0; jA < NumCols(A); ++jA)
  //       for (long iB=0; iB < NumRows(B); ++iB)
  //         for (long jB=0; jB < NumCols(B); ++jB)
  //           SetEntry(ans, iA*NumRows(B)+iB, jA*NumCols(B)+jB, A(iA, jA)*B(iB, jB));
  //   return ans;
  // }



} // end of namespace CoCoA
