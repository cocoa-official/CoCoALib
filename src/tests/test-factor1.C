//   Copyright (c)  2006,2017  John Abbott,  Anna M. Bigatti

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



//----------------------------------------------------------------------
// WARNING: the user interface to the factorizer is not final!
// Of course, this is only a minimal test.
//----------------------------------------------------------------------

#include "CoCoA/BuildInfo.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/factor.H"


#include <vector>
using std::vector;
#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  void TestQQ()
  {
    const ring QQ = RingQQ();
    SparsePolyRing P = NewPolyRing(QQ, symbols("x,y"));

    const vector<RingElem>& x = indets(P);
    const RingElem f = 5*(power(2*x[0],96) - power(3*x[1],96));
    const factorization<RingElem> FacInfo = factor(f);
    const vector<RingElem>& facs = FacInfo.myFactors();   // handy alias
    const vector<long>& mult = FacInfo.myMultiplicities();// handy alias
    const int NumFacs = len(facs);
    CoCoA_ASSERT_ALWAYS(NumFacs == 12);
    RingElem prod = FacInfo.myRemainingFactor();
    for (int i = 0; i < NumFacs; ++i)
    {
      CoCoA_ASSERT_ALWAYS(mult[i] == 1);
      prod *= facs[i];
    }
    CoCoA_ASSERT_ALWAYS(prod == f);
  }

  void TestZZ()
  {
    const ring ZZ = RingZZ();
    SparsePolyRing P = NewPolyRing(ZZ, symbols("x,y"));

    const vector<RingElem>& x = indets(P);
    const RingElem f = 5*(power(2*x[0],96) - power(3*x[1],96));
    const factorization<RingElem> FacInfo = factor(f);
    const vector<RingElem>& facs = FacInfo.myFactors();   // handy alias
    const vector<long>& mult = FacInfo.myMultiplicities();// handy alias
    const int NumFacs = len(facs);
    CoCoA_ASSERT_ALWAYS(NumFacs == 12);
    RingElem prod = FacInfo.myRemainingFactor();
    for (int i = 0; i < NumFacs; ++i)
    {
      CoCoA_ASSERT_ALWAYS(mult[i] == 1);
      prod *= facs[i];
    }
    CoCoA_ASSERT_ALWAYS(prod == f);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    TestQQ();
    TestZZ();
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-factor1.C,v 1.14 2022/02/11 09:54:00 abbott Exp $
// $Log: test-factor1.C,v $
// Revision 1.14  2022/02/11 09:54:00  abbott
// Summary: Updated copyright (redmine 855)
//
// Revision 1.13  2018/09/28 15:54:05  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.12  2017/06/22 11:10:33  abbott
// Summary: Added test for factorization over ZZ; made slightly harder (non-triv. coeffs)
//
// Revision 1.11  2016/09/16 16:36:41  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.10  2015/05/07 15:00:40  abbott
// Summary: Moved test code into namespace CoCoA
// Author: JAA
//
// Revision 1.9  2014/04/11 15:08:24  abbott
// Summary: Renamed TmpFactor to factor in include
// Author: JAA
//
// Revision 1.8  2014/03/24 12:11:34  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.7  2012/10/05 09:33:20  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.6  2012/02/10 11:57:12  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.5  2012/02/08 17:36:48  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.3  2009/07/24 14:23:03  abbott
// Modified as result of changed interface to factorizer.
//
// Revision 1.2  2007/10/30 17:14:05  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.6  2007/03/05 13:42:58  bigatti
// -- changed: no longer depends on order of factors (removed test-factor1.out)
//
// Revision 1.5  2007/03/03 14:13:21  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.3  2007/02/26 17:11:58  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1  2006/11/11 14:51:43  cocoa
// Added simple test for the factorizer.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
