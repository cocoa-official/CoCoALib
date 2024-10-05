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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-gsl1.C,v 1.7 2018/08/28 12:32:52 abbott Exp $
// $Log: test-gsl1.C,v $
// Revision 1.7  2018/08/28 12:32:52  abbott
// Summary: Updated to use BigRatFromString
//
// Revision 1.6  2016/09/16 16:36:41  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.5  2015/05/07 15:00:40  abbott
// Summary: Moved test code into namespace CoCoA
// Author: JAA
//
// Revision 1.4  2014/09/05 16:12:30  bigatti
// -- removed printed output
//
// Revision 1.3  2014/09/05 14:24:44  bigatti
// -- this compiles ;-)
//
// Revision 1.2  2014/09/02 16:38:33  bigatti
// -- added prefix "Gsl"
//
// Revision 1.1  2011/03/21 17:40:16  bigatti
// -- first import
//
// Revision 1.7  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.6  2010/11/26 15:32:16  bigatti
// -- renamed Dimension into dimension (coding conventions)
//
// Revision 1.5  2010/11/22 17:39:15  bigatti
// -- updated "TmpFrobby" --> "ExternalLibs-Frobby"
//
// Revision 1.4  2010/02/04 09:38:28  bigatti
// -- new syntax for frobby (more CoCoALib-like)
//
// Revision 1.3  2010/02/03 18:00:22  bigatti
// -- more functions from frobby
//
// Revision 1.2  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
//
// Revision 1.1  2009/02/11 15:08:26  bigatti
// -- first import
//

