// Copyright (c) 2012 John Abbott, Anna Bigatti

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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PPMonoidEvOv.H"
#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/PPMonoidOv.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
//#include "CoCoA/DenseMatrix.H"
//#include "CoCoA/matrix.H"

#include<iostream>
using std::cout;
using std::cerr;
using std::endl;
#include<vector>
using std::vector;

//----------------------------------------------------------------------
// First basic test for GeneralHom and RestrictionHom
//----------------------------------------------------------------------
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    PPMonoid PPM1 = NewPPMonoid(symbols("x,y,z"), StdDegRevLex);
    PPMonoid PPM2 = NewPPMonoid(symbols("alpha,beta,gamma,delta"), lex);

    const vector<PPMonoidElem>& x = indets(PPM1);
    const vector<PPMonoidElem>& alpha = indets(PPM2);

    vector<PPMonoidElem> images;
    images.push_back(alpha[0]*alpha[1]);
    images.push_back(alpha[1]*alpha[2]);
    images.push_back(alpha[2]*power(alpha[3],4));

    PPMonoidHom phi = GeneralHom(PPM1, images);
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
        {
          PPMonoidElem t = power(x[0],i) * power(x[1],j) * power(x[2],k);
          CoCoA_ASSERT_ALWAYS(phi(t) == power(images[0],i) * power(images[1],j) * power(images[2],k));
        }

    PPMonoidHom psi = RestrictionHom(PPM1, vector<long>(1,1));  
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
        {
          PPMonoidElem t = power(x[0],i) * power(x[1],j) * power(x[2],k);
          CoCoA_ASSERT_ALWAYS(psi(t) == power(x[1],j));
        }
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-PPMonoidHom1.C,v 1.7 2018/05/18 12:31:54 bigatti Exp $
// $Log: test-PPMonoidHom1.C,v $
// Revision 1.7  2018/05/18 12:31:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.6  2016/09/16 16:36:39  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.5  2016/08/02 16:44:53  abbott
// Summary: Removed include of obsolete header file (PPMonoidEvZZ.H)
//
// Revision 1.4  2015/07/27 12:45:47  bigatti
// -- now using "symbols(string)" comma separated symbols
//
// Revision 1.3  2015/05/07 15:00:39  abbott
// Summary: Moved test code into namespace CoCoA
// Author: JAA
//
// Revision 1.2  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.1  2012/02/14 16:01:00  bigatti
// -- fist import
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2010/07/09 17:03:03  abbott
// Example for PPMonoid homs
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
