//   Copyright (c)  2011  John Abbott, Anna Bigatti

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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixSpecial.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-resultant.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"


#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for special matrices:
// functions: JacobianMat
//----------------------------------------------------------------------
namespace CoCoA
{

void TestSylvesterMat(const PolyRing& P)
{
  RingElem f(P, "x^2+y^2-1");
  RingElem g(P, "x^4+y^4-1");
  RingElem y(P, "y");

  const int degf = deg(f,1); // y = indet(P,1)
  const int degg = deg(g,1); // y = indet(P,1)

  matrix M = SylvesterMat(f,g,y);
  CoCoA_ASSERT_ALWAYS(NumRows(M)==degf+degg && NumCols(M)==NumRows(M));
  CoCoA_ASSERT_ALWAYS(det(M) == resultant(f,g,1)); // y = indet(P,1)
}


void TestJacobianMat(const PolyRing& P)
{
  std::vector<RingElem> f;
  long n = NumIndets(P);
  
  matrix M = JacobianMat(indets(P));
  CoCoA_ASSERT_ALWAYS(NumRows(M)==n && NumCols(M)==n);
  CoCoA_ASSERT_ALWAYS(IsOne(M(0,0)) && IsOne(M(1,1)));
  CoCoA_ASSERT_ALWAYS(IsZero(M(0,1)) && IsZero(M(1,0)));

  M = JacobianMat(indets(P), indets(P));
  CoCoA_ASSERT_ALWAYS(NumRows(M)==n && NumCols(M)==n);
  CoCoA_ASSERT_ALWAYS(IsOne(M(0,0)) && IsOne(M(1,1)));
  CoCoA_ASSERT_ALWAYS(IsZero(M(0,1)) && IsZero(M(1,0)));

  M = JacobianMat(indets(P), f);
  CoCoA_ASSERT_ALWAYS(NumRows(M)==n && NumCols(M)==0);

  M = JacobianMat(f, indets(P));
  CoCoA_ASSERT_ALWAYS(NumRows(M)==0 && NumCols(M)==n);
  //  std::cout << JacobianMat(f, f) << endl; ---> CoCoA_ERROR
}


void TestKroneckerProd(const ring& R)
{
  {
    matrix M = KroneckerProd(IdentityMat(RingZZ(), 2), IdentityMat(RingZZ(), 2));
    CoCoA_ASSERT_ALWAYS(M == IdentityMat(RingZZ(), 4));
  }
  {
    ConstMatrixView I1x1 = IdentityMat(R, 1);
    ConstMatrixView A = BlockMat2x2(2*I1x1, 3*I1x1, 5*I1x1, 7*I1x1);
    matrix M = KroneckerProd(A, A);
    std::vector<long> ZeroOne;
    std::vector<long> TwoThree;
    ZeroOne.push_back(0);    ZeroOne.push_back(1);
    TwoThree.push_back(2);   TwoThree.push_back(3);
    CoCoA_ASSERT_ALWAYS(submat(M, ZeroOne, ZeroOne) == 2 * A);
    CoCoA_ASSERT_ALWAYS(submat(M, ZeroOne, TwoThree) == 3 * A);
    CoCoA_ASSERT_ALWAYS(submat(M, TwoThree, ZeroOne) == 5 * A);
    CoCoA_ASSERT_ALWAYS(submat(M, TwoThree, TwoThree) == 7 * A);
  }
  
}


void program()
{
  GlobalManager CoCoAFoundations;

  PolyRing P = NewPolyRing(RingQQ(), symbols("x,y"));
  TestSylvesterMat(P);
  TestJacobianMat(P);
  TestKroneckerProd(P);
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
