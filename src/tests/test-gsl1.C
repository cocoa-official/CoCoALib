//   Copyright (c)  2011 Anna Bigatti

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#ifdef CoCoA_WITH_GSL
#include "CoCoA/BigRat.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/ExternalLibs-GSL.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/matrix.H"
#endif


#include <iostream>
using std::cerr;
using std::endl;
#ifdef CoCoA_WITH_GSL
#include <vector>
using std::vector;
#endif
//----------------------------------------------------------------------
// First test for GSL library.
// Trivial computation to test proper integration.
//----------------------------------------------------------------------
namespace CoCoA
{

void program()
{
  GlobalManager CoCoAFoundations;

#ifdef CoCoA_WITH_GSL

  const int m = 4;
  const int n = 5;
  
  // C matrix
  int C_matrix[m][n] = {{1,0,0,0,2},
                        {0,0,3,0,0},
                        {0,0,0,0,0},
                        {0,4,0,0,0}};

  matrix A(NewDenseMat(RingQQ(),n,m));
  for (int i=0; i < m; i++)
    for (int j=0; j < n; j++)
      SetEntry(A, j, i, C_matrix[i][j]);

  // Calculate SVD
  vector<matrix> r = GslSVD(A);
  //  std::cout << "A = " << A << std::endl;
  //  std::cout << "U = " << r[0] << std::endl;
  //  std::cout << "Elements of V in untransposed form " << r[2] << std::endl;
  //  std::cout << "Singular Values " << r[1] << std::endl;
  CoCoA_ASSERT_ALWAYS(r[1](0,0)==4);
  CoCoA_ASSERT_ALWAYS(r[1](0,1)==3);
  CoCoA_ASSERT_ALWAYS(r[1](0,2)==BigRatFromString("629397181890197/281474976710656"));
  CoCoA_ASSERT_ALWAYS(r[1](0,3)==0);
#endif
}

} // end of namespace CoCoA

//----------------------------------------------------------------------
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
