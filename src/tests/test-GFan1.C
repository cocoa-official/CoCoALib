// Copyright (c) 2015 John Abbott, Anna Bigatti, Anders Nedergaard Jensen

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
#include "CoCoA/ExternalLibs-GFan.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/matrix.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// First test for GFan library.
// Simple computation to test proper integration.
// This tests is for the functions working directly on cone.
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

#ifdef CoCoA_WITH_GFAN
    using namespace CoCoA::GFan;

    matrix CM(NewDenseMat(RingQQ(),2,3));
    SetEntry(CM,0,0, 300);  //     M[0][0] = 300; 
    SetEntry(CM,0,1, 2);    //     M[0][1] = 2;   
    SetEntry(CM,0,2, -2);   //     M[0][2] = -2;  
    SetEntry(CM,1,0, 30);   //     M[1][0] = 30;  
    SetEntry(CM,1,1, 5);    //     M[1][1] = 5;   
    SetEntry(CM,1,2, 1);    //     M[1][2] = 1;

    matrix CM2(NewDenseMat(RingQQ(),1,3));
    SetEntry(CM,0,0, 4);
    SetEntry(CM,0,1, 1);
    SetEntry(CM,0,2, 5);

    GFan::cone CC(CM, CM2);
    //    cout << "The GFan cone is " << CC << endl;

    using std::vector;
    CoCoA_ASSERT_ALWAYS(equations(CC) == ZeroMat(RingZZ(), 0, 3));
    matrix M = NewDenseMat(RingZZ(), 
                           vector<vector<long>>{{4L,1L,5L},{30L,5L,1L}});
    CoCoA_ASSERT_ALWAYS(inequalities(CC) == M);
    
    CoCoA_ASSERT_ALWAYS(transpose(RelativeInteriorPoint(CC)) ==
                        NewDenseMat(RingZZ(), 
                                    vector<vector<long>>{{-2L,13L,0L}}));
#endif // CoCoA_WITH_GFAN
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-GFan1.C,v 1.6 2020/01/24 21:39:59 abbott Exp $
// $Log: test-GFan1.C,v $
// Revision 1.6  2020/01/24 21:39:59  abbott
// Summary: Replaced CoCoAVector by braced ctor
//
// Revision 1.5  2018/05/17 16:08:42  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.4  2017/03/08 15:58:53  bigatti
// -- cleaned up
//
// Revision 1.3  2016/09/21 14:57:22  abbott
// Summary: Added missing includes
//
// Revision 1.2  2016/09/19 06:37:57  bigatti
// -- simplified (CoCoA_ASSERT_ALWAYS)
//
// Revision 1.1  2016/09/19 06:01:55  bigatti
// -- first import
//
