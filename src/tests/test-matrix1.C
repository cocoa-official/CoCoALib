//   Copyright (c)  2007,2008  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/time.H"


#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for matrix
// functions: IdentityMat, submat, BlockMat2x2, ConcatVer, ConcatHor
//            operator==, IsSymmetric, IsAntiSymmetric, IsDiagonal
//----------------------------------------------------------------------

namespace CoCoA
{

  void TestZeroMat(ring R)
  {
    //  std::cout << " --TestZeroMat" << std::endl;
    // A few quick checks on a (large) zero matrix.
    ConstMatrixView ZM = ZeroMat(R, 12345,67890);
    CoCoA_ASSERT_ALWAYS(IsZero(ZM));
    for (long i=0; i < NumRows(ZM); ++i)
      CoCoA_ASSERT_ALWAYS(IsZeroRow(ZM, i));
    for (long j=0; j < NumCols(ZM); ++j)
      CoCoA_ASSERT_ALWAYS(IsZeroCol(ZM, j));
    for (long i=0; i < min(NumRows(ZM), NumCols(ZM)); ++i)
      CoCoA_ASSERT_ALWAYS(IsZero(ZM(i, i)));
    CoCoA_ASSERT_ALWAYS(ZM == ZM);
    CoCoA_ASSERT_ALWAYS(IsSymmetric(ZeroMat(R, 1000000000,1000000000)));
    CoCoA_ASSERT_ALWAYS(IsAntiSymmetric(ZeroMat(R, 1000000000,1000000000)));
    CoCoA_ASSERT_ALWAYS(IsDiagonal(ZeroMat(R, 1000000000,1000000000)));
  }


  void TestIdentityMat(ring R)
  {
    //  std::cout << " --TestIdentityMat" << std::endl;
    // An IdentityMatrix is constant; you cannot change it, but you can copy it.
    ConstMatrixView Id = IdentityMat(R, 10);
    CoCoA_ASSERT_ALWAYS(NumCols(Id) == 10);
    CoCoA_ASSERT_ALWAYS(NumCols(Id) == NumRows(Id));
    // SetEntry(Id,1,1,1); <-- as expected, this does not compile because I10 is constant
    for (long i=0; i < NumRows(Id); ++i)
      for (long j=0; j < NumCols(Id); ++j)
      {
        if (i == j) CoCoA_ASSERT_ALWAYS(IsOne(Id(i, j)));
        else CoCoA_ASSERT_ALWAYS(IsZero(Id(i, j)));
      }
    CoCoA_ASSERT_ALWAYS(Id != ZeroMat(R, 10, 10));
    CoCoA_ASSERT_ALWAYS(ZeroMat(R, 10, 10) != Id);
    CoCoA_ASSERT_ALWAYS(Id == IdentityMat(R, 10));
    CoCoA_ASSERT_ALWAYS(IsSymmetric(IdentityMat(R, 1000000000)));
    CoCoA_ASSERT_ALWAYS(!IsAntiSymmetric(IdentityMat(R, 1000000000)));
    CoCoA_ASSERT_ALWAYS(IsDiagonal(IdentityMat(R, 1000000000)));

    matrix DMId = NewDenseMat(Id);
    CoCoA_ASSERT_ALWAYS(Id == DMId);
    CoCoA_ASSERT_ALWAYS(IsSymmetric(DMId));
    CoCoA_ASSERT_ALWAYS(!IsAntiSymmetric(DMId));
    CoCoA_ASSERT_ALWAYS(IsDiagonal(DMId));

    // Create a submatrix view of I10; since the row and column
    // indices are the same lists, the result will be another
    // identity matrix.
    vector<long> RowIndices;
    RowIndices.push_back(2);
    RowIndices.push_back(7);
    RowIndices.push_back(1);
    RowIndices.push_back(8);
    vector<long> ColIndices = RowIndices;
    ConstMatrixView SubM = submat(Id, RowIndices, ColIndices);
    CoCoA_ASSERT_ALWAYS(SubM == IdentityMat(R, 4));
  }


  void TestDiagMat(ring R)
  {
    //  std::cout << " --TestDiagMat" << std::endl;
    std::vector<RingElem> v(10, one(R));
    ConstMatrixView D = DiagMat(v);
    CoCoA_ASSERT_ALWAYS(!IsZero(D));
    CoCoA_ASSERT_ALWAYS(D == D);
    CoCoA_ASSERT_ALWAYS(D == NewDenseMat(D));
    CoCoA_ASSERT_ALWAYS(D == IdentityMat(R, 10));
    CoCoA_ASSERT_ALWAYS(IdentityMat(R, 10) == D);
    CoCoA_ASSERT_ALWAYS(IsSymmetric(D));
    CoCoA_ASSERT_ALWAYS(!IsAntiSymmetric(D));
    CoCoA_ASSERT_ALWAYS(IsDiagonal(D));
    // SetEntry(D,1,1,1); <-- as expected, this does not compile because D is constant
    std::vector<RingElem> v0(10, zero(R));
    MatrixView D0 = DiagMat(v0);
    CoCoA_ASSERT_ALWAYS(D0 == ZeroMat(R, 10, 10));
    CoCoA_ASSERT_ALWAYS(IsAntiSymmetric(D0));
    CoCoA_ASSERT_ALWAYS(ZeroMat(R, 10, 10) == D0);
    // SetEntry(D0,1,2,5); //<-- as expected, this throws at run-time
  }


  void TestBlockMat(ring R)
  {
    //  std::cout << " --TestBlockMat" << std::endl;
    // Create a 6x6 matrix M containing the numbers 0 to 35
    matrix M(NewDenseMat(R, 6, 6));
    for (long i=0; i < 6; ++i)
      for (long j=0; j < 6; ++j)
        SetEntry(M, i, j, 6*i+j);

    CoCoA_ASSERT_ALWAYS(!IsSymmetric(M));
    CoCoA_ASSERT_ALWAYS(!IsAntiSymmetric(M));
    CoCoA_ASSERT_ALWAYS(!IsDiagonal(M));
    CoCoA_ASSERT_ALWAYS(IsSymmetric(M*transpose(M)));
  
    vector<long> R1, R2, R3;
    R1 = LongRange(0,0);
    R2 = LongRange(1,2);
    R3 = LongRange(3,5);
    vector<long> C1 = R1;
    vector<long> C2 = R2;
    vector<long> C3 = R3;

    // Cut M up into 9 submatrices
    ConstMatrixView M11(submat(M, R1, C1));
    ConstMatrixView M12(submat(M, R1, C2));
    ConstMatrixView M13(submat(M, R1, C3));

    ConstMatrixView M21(submat(M, R2, C1));
    ConstMatrixView M22(submat(M, R2, C2));
    ConstMatrixView M23(submat(M, R2, C3));

    ConstMatrixView M31(submat(M, R3, C1));
    ConstMatrixView M32(submat(M, R3, C2));
    ConstMatrixView M33(submat(M, R3, C3));

    // Reassemble the submatrices into a single matrix (method 1)
    ConstMatrixView Mcopy1 = BlockMat2x2(BlockMat2x2(M11, M12, M21, M22),
                                         ConcatVer(M13, M23),
                                         ConcatHor(M31, M32),
                                         M33);

    // Reassemble the submatrices into a single matrix (method 2)
    ConstMatrixView Mcopy2 = BlockMat2x2(M11,
                                         ConcatHor(M12, M13),
                                         ConcatVer(M21, M31),
                                         BlockMat2x2(M22, M23, M32, M33));

    // Check that the reassembled matrices are equal to M.
    CoCoA_ASSERT_ALWAYS(M == Mcopy1);
    CoCoA_ASSERT_ALWAYS(M == Mcopy2);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    const ring ZZ = RingZZ();

    TestZeroMat(ZZ);
    TestZeroMat(RingQQ());
    TestZeroMat(NewZZmod(2));

    TestIdentityMat(ZZ);
    TestIdentityMat(RingQQ());
    TestIdentityMat(NewZZmod(7));

    TestDiagMat(ZZ);

    TestBlockMat(ZZ);

    matrix M(NewDenseMat(ZZ, 2, 3));
    SetEntry(M,1,0,BigInt(3));
    SetEntry(M,1,1,1);
    SetEntry(M,1,2,BigRat(4,1));

    CoCoA_ASSERT_ALWAYS(IsZero(M(0,0)));

    CoCoA_ASSERT_ALWAYS(M->myIsZeroRow(0));
    CoCoA_ASSERT_ALWAYS(!M->myIsZeroRow(1));

    M->myAddRowMul(0, 1, RingElem(ZZ, 2));
    CoCoA_ASSERT_ALWAYS(!M->myIsZeroRow(0));
    CoCoA_ASSERT_ALWAYS(M(0,1) == 2);
    M->myAddRowMul(0, 1, M(0,0)); // alias test
    CoCoA_ASSERT_ALWAYS(M(0,0) == 8*M(1,0));
    CoCoA_ASSERT_ALWAYS(M(0,1) == 8*M(1,1));
    CoCoA_ASSERT_ALWAYS(M(0,2) == 8*M(1,2));

    // Some simple tests on the transposed view of a matrix.
    MatrixView TrM = transpose(M);
    CoCoA_ASSERT_ALWAYS(NumRows(TrM) == NumCols(M));
    CoCoA_ASSERT_ALWAYS(NumRows(M) == NumCols(TrM));
    MatrixView TrTrM = transpose(TrM);
    CoCoA_ASSERT_ALWAYS(TrTrM == M);
    for (long i = 0; i < NumRows(M); ++i)
      for (long j = 0 ; j < NumCols(M); ++j)
        SetEntry(TrTrM,i,j, 2*M(i,j));
    CoCoA_ASSERT_ALWAYS(TrTrM == M);
  }

} // end of namespace CoCoA


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
