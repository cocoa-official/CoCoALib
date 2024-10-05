//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingWeyl.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;



//----------------------------------------------------------------------
// Basic operations on RingWeyl, the interface is not quite settled yet
// functions: *
//----------------------------------------------------------------------
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    vector<symbol> names;
    names.push_back(symbol("x"));
    names.push_back(symbol("y"));
    names.push_back(symbol("z"));
    names.push_back(symbol("t"));  

    long NumInds = names.size();
    vector<long> ElimIndets;
    SparsePolyRing WA = NewWeylAlgebra(RingQQ(), names, ElimIndets);
    CoCoA_ASSERT_ALWAYS(!IsCommutative(WA));
  
    RingElem x = indet(WA, 0);
    RingElem y = indet(WA, 1);
    RingElem z = indet(WA, 2);
    RingElem t = indet(WA, 3);

    RingElem dx = indet(WA, 0+NumInds);
    RingElem dy = indet(WA, 1+NumInds);
    RingElem dz = indet(WA, 2+NumInds);
    RingElem dt = indet(WA, 3+NumInds);

    CoCoA_ASSERT_ALWAYS(dx*x == x*dx + 1);
    CoCoA_ASSERT_ALWAYS(dy*y == y*dy + 1);
    CoCoA_ASSERT_ALWAYS(dz*z == z*dz + 1);
    CoCoA_ASSERT_ALWAYS(dt*t == t*dt + 1);

    const vector<RingElem> GB = TidyGens(ideal(x, dx));
    
    CoCoA_ASSERT_ALWAYS("GBasis" && len(GB)==1);
    CoCoA_ASSERT_ALWAYS("GBasis" && GB[0]==1);
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

//output
//gens(I)=[
//-x[0]*x[1] +(1),
//x[0]*x[5] +-x[2]^3 +-x[3]^3 +-x[4]^3,
//(3)*x[1]*x[2]^2*d[5] +d[2],
//(3)*x[1]*x[3]^2*d[5] +d[3],
//(3)*x[1]*x[4]^2*d[5] +d[4]]
 
//TidyGens(I)=[
//(3)*x[1]*x[4]^2*d[5] +d[4],
//(3)*x[1]*x[3]^2*d[5] +d[3],
//-x[4]^2*d[3] +x[3]^2*d[4],
//(3)*x[1]*x[2]^2*d[5] +d[2],
//-x[4]^2*d[2] +x[2]^2*d[4],
//x[3]^2*d[2]*d[4] +-x[2]^2*d[3]*d[4],
//(-2)*x[3]^2*x[4]*d[2] +(2)*x[2]^2*x[4]*d[3],
//-x[3]^2*d[2] +x[2]^2*d[3],
//x[0]*x[5] +-x[2]^3 +-x[3]^3 +-x[4]^3,
//(-3)*x[0]*x[1]*x[4]^2 +(3)*x[4]^2,
//x[0]*d[4] +(-3)/(-1) *x[4]^2*d[5],
//-x[2]^3*d[4] +-x[3]^3*d[4] +-x[4]^3*d[4] +(3)/(-1) *x[4]^2*x[5]*d[5] 
//+(-3)*x[4]^2,
//(-2)*x[2]^3*x[4]*d[3] +(-2)*x[3]^3*x[4]*d[3] +(-2)*x[3]^2*x[4]^2*d[4] 
//+(-6)*x[3]^2*x[4]*x[5]*d[5] +(6)/(-1) *x[3]^2*x[4],
//x[2]^3*d[3] +x[3]^3*d[3] +x[3]^2*x[4]*d[4] +(-3)/(-1) *x[3]^2*x[5]*d[5] +(3)*x[3]^2,
//(-2)*x[2]^3*x[4]*d[2] +(2)/(-1) *x[2]^2*x[3]*x[4]*d[3] +(-2)*x[2]^2*x[4]^2*d[4]
//+(-6)*x[2]^2*x[4]*x[5]*d[5] +(6)/(-1) *x[2]^2*x[4],
//x[2]^3*d[2] +x[2]^2*x[3]*d[3] +x[2]^2*x[4]*d[4] +(-3)/(-1) *x[2]^2*x[5]*d[5] 
//+(3)*x[2]^2,
//x[2]*d[2]*d[4] +x[3]*d[3]*d[4] +x[4]*d[4]^2 +(3)*x[5]*d[4]*d[5] +(-4)/(-1) *d[4],
//(-2)*x[2]*x[4]*d[2] +(-2)*x[3]*x[4]*d[3] +(-2)*x[4]^2*d[4] +(-6)*x[4]*x[5]*d[5]
//+(6)/(-1) *x[4],
//-x[2]*d[2] +-x[3]*d[3] +-x[4]*d[4] +(3)/(-1) *x[5]*d[5] +(-3),
//(3)/(-1) *x[1]*x[2]*x[3]*d[3]*d[5] +(3)/(-1) *x[1]*x[2]*x[4]*d[4]*d[5] 
//+(9)/(-1) *x[1]*x[2]*x[5]*d[5]^2 +(12)/(-1) *x[1]*x[2]*d[5] +d[2]^2,
//(-2)*x[1]*x[3]*x[4]*d[3]*d[4]*d[5] +(-6)*x[1]*x[3]*x[5]*d[3]*d[5]^2 
//+(6)/(-1) *x[1]*x[4]*x[5]*d[4]*d[5]^2 +(9)/(-1) *x[1]*x[5]^2*d[5]^3 
//+(-6)*x[1]*x[3]*d[3]*d[5] +(6)/(-1) *x[1]*x[4]*d[4]*d[5] +(-36)*x[1]*x[5]*d[5]^2 
//+(16)/(-1) *x[1]*d[5] +(1)/(-3) *d[2]^3 +(-1)/(-3) *d[3]^3 +(-1)/(-3) *d[4]^3,
//(3)/(-1) *x[1]*x[2]*x[3]*x[4]*d[4]*d[5] +(9)/(-1) *x[1]*x[2]*x[3]*x[5]*d[5]^2 
//+(6)/(-1) *x[1]*x[2]*x[3]*d[5] +x[3]*d[2]^2 +x[2]*d[3]^2,
//(9)/(-1) *x[1]*x[3]*x[4]*x[5]*d[3]*d[5]^2 +(-27)/(2) *x[1]*x[4]*x[5]^2*d[5]^3 
//+(3)/(-1) *x[1]*x[3]*x[4]*d[3]*d[5] +(36)/(-1) *x[1]*x[4]*x[5]*d[5]^2 
//+(-6)*x[1]*x[4]*d[5] +(-1)/(2) *x[4]*d[2]^3 +(1)/(2) *x[4]*d[3]^3 +x[3]*d[3]*d[4]^2 
//+(1)/(2) *x[4]*d[4]^3 +(3)*x[5]*d[4]^2*d[5] +(3)*d[4]^2,
//(-9)*x[1]*x[3]*x[4]*x[5]*d[4]*d[5]^2 +(-27)/(2) *x[1]*x[3]*x[5]^2*d[5]^3 
//+(-3)*x[1]*x[3]*x[4]*d[4]*d[5] +(-36)*x[1]*x[3]*x[5]*d[5]^2 +(6)/(-1) *x[1]*x[3]*d[5] 
//+(-1)/(2) *x[3]*d[2]^3 +(1)/(2) *x[3]*d[3]^3 +x[4]*d[3]^2*d[4] +(1)/(2) *x[3]*d[4]^3 
//+(-3)/(-1) *x[5]*d[3]^2*d[5] +(-3)/(-1) *d[3]^2,
//(9)/(-1) *x[1]*x[2]*x[4]*x[5]*d[4]*d[5]^2 +(-27)/(2) *x[1]*x[2]*x[5]^2*d[5]^3 
//+(3)/(-1) *x[1]*x[2]*x[4]*d[4]*d[5] +(36)/(-1) *x[1]*x[2]*x[5]*d[5]^2 
//+(6)/(-1) *x[1]*x[2]*d[5] +(-1)/(2) *x[3]*d[2]^2*d[3] +(1)/(-2) *x[2]*d[3]^3 
//+(1)/(2) *x[4]*d[2]^2*d[4] +(1)/(2) *x[2]*d[4]^3 +(3)/(2) *x[5]*d[2]^2*d[5] 
//+(1)/(2) *d[2]^2,
//(9)/(-1) *x[1]*x[2]*x[3]*x[4]*x[5]*d[5]^2 +x[3]*x[4]*d[2]^2 +x[2]*x[4]*d[3]^2 
//+x[2]*x[3]*d[4]^2,
//(-6)*x[1]*x[3]*x[5]^2*d[3]*d[5]^3 +(6)/(-1) *x[1]*x[5]^3*d[5]^4 
//+(12)/(-1) *x[1]*x[3]*x[5]*d[3]*d[5]^2 +(36)/(-1) *x[1]*x[5]^2*d[5]^3 
//+(4)/(-3) *x[1]*x[3]*d[3]*d[5] +(-112)/(3) *x[1]*x[5]*d[5]^2 +(8)/(-3) *x[1]*d[5] 
//+(2)/(-27) *x[3]*d[2]^3*d[3] +(-2)/(-27) *x[3]*d[3]^4 +(4)/(27) *x[4]*d[2]^3*d[4] 
//+(-2)/(9) *x[3]*d[3]*d[4]^3 +(-4)/(27) *x[4]*d[4]^4 +(2)/(-9) *x[5]*d[2]^3*d[5] 
//+(-2)/(-3) *x[5]*d[3]^3*d[5] +(2)/(-3) *x[5]*d[4]^3*d[5] +(-4)/(-9) *d[3]^3 +(-8)/(9) *d[4]^3,
//(-9)/(2) *x[1]*x[4]*x[5]^2*d[4]*d[5]^3 +(-9)/(2) *x[1]*x[5]^3*d[5]^4 
//+(-9)*x[1]*x[4]*x[5]*d[4]*d[5]^2 +(27)/(-1) *x[1]*x[5]^2*d[5]^3 +-x[1]*x[4]*d[4]*d[5] 
//+(-28)*x[1]*x[5]*d[5]^2 +(-2)*x[1]*d[5] +(-1)/(-9) *x[3]*d[2]^3*d[3] +(1)/(-9) *x[3]*d[3]^4 
//+(-1)/(18) *x[4]*d[2]^3*d[4] +(1)/(-6) *x[4]*d[3]^3*d[4] +(1)/(18) *x[4]*d[4]^4 
//+(-1)/(6) *x[5]*d[2]^3*d[5] +(-1)/(2) *x[5]*d[3]^3*d[5] +(1)/(2) *x[5]*d[4]^3*d[5] 
//+(2)/(-3) *d[3]^3 +(-1)/(-3) *d[4]^3,
//(9)/(-2) *x[1]*x[3]*x[4]*x[5]^2*d[5]^3 +(6)/(-1) *x[1]*x[3]*x[4]*x[5]*d[5]^2 
//+(1)/(-6) *x[3]*x[4]*d[2]^3 +(-1)/(-6) *x[3]*x[4]*d[3]^3 +(1)/(3) *x[3]^2*d[3]*d[4]^2 
//+(-1)/(-6) *x[3]*x[4]*d[4]^3 +x[4]*x[5]*d[3]^2*d[5] +x[3]*x[5]*d[4]^2*d[5] 
//+(1)/(3) *x[4]*d[3]^2 +x[3]*d[4]^2,
//(-9)/(2) *x[1]*x[2]*x[4]*x[5]^2*d[5]^3 +(6)/(-1) *x[1]*x[2]*x[4]*x[5]*d[5]^2 
//+(-1)/(6) *x[3]*x[4]*d[2]^2*d[3] +(1)/(-6) *x[2]*x[4]*d[3]^3 
//+(-1)/(6) *x[2]*x[3]*d[3]*d[4]^2 +(1)/(2) *x[4]*x[5]*d[2]^2*d[5] 
//+(-1)/(-2) *x[2]*x[5]*d[4]^2*d[5] +(1)/(-6) *x[4]*d[2]^2 +(1)/(-6) *x[2]*d[4]^2,
//(9)/(-2) *x[1]*x[2]*x[3]*x[5]^2*d[5]^3 +(6)/(-1) *x[1]*x[2]*x[3]*x[5]*d[5]^2 
//+(1)/(-6) *x[3]*x[4]*d[2]^2*d[4] +(1)/(-6) *x[2]*x[4]*d[3]^2*d[4] 
//+(-1)/(6) *x[2]*x[3]*d[4]^3 +(-1)/(-2) *x[3]*x[5]*d[2]^2*d[5] 
//+(-1)/(-2) *x[2]*x[5]*d[3]^2*d[5] +(1)/(-6) *x[3]*d[2]^2 +(-1)/(6) *x[2]*d[3]^2,
//(3)/(-1) *x[1]*x[4]*x[5]^3*d[5]^4 +(-12)*x[1]*x[4]*x[5]^2*d[5]^3 
//+(-20)/(3) *x[1]*x[4]*x[5]*d[5]^2 +(2)/(27) *x[3]*x[4]*d[2]^3*d[3] 
//+(-2)/(27) *x[3]*x[4]*d[3]^4 +(4)/(-27) *x[3]^2*d[3]^2*d[4]^2 
//+(-2)/(27) *x[3]*x[4]*d[3]*d[4]^3 +(1)/(-9) *x[4]*x[5]*d[2]^3*d[5] 
//+(1)/(-3) *x[4]*x[5]*d[3]^3*d[5] +(-2)/(9) *x[3]*x[5]*d[3]*d[4]^2*d[5] 
//+(-1)/(-9) *x[4]*x[5]*d[4]^3*d[5] +(2)/(3) *x[5]^2*d[4]^2*d[5]^2 +(2)/(27) *x[4]*d[2]^3 
//+(2)/(-9) *x[4]*d[3]^3 +(-20)/(27) *x[3]*d[3]*d[4]^2 +(2)/(-27) *x[4]*d[4]^3 
//+(-8)/(-9) *x[5]*d[4]^2*d[5] +(4)/(-9) *d[4]^2,
//(9)/(-2) *x[1]*x[3]*x[5]^3*d[5]^4 +(-18)*x[1]*x[3]*x[5]^2*d[5]^3 
//+(-10)*x[1]*x[3]*x[5]*d[5]^2 +(-1)/(-9) *x[3]*x[4]*d[2]^3*d[4] 
//+(-1)/(9) *x[3]*x[4]*d[3]^3*d[4] +(-2)/(9) *x[3]^2*d[3]*d[4]^3 
//+(1)/(-9) *x[3]*x[4]*d[4]^4 +(1)/(-6) *x[3]*x[5]*d[2]^3*d[5] 
//+(-1)/(-6) *x[3]*x[5]*d[3]^3*d[5] +(1)/(-3) *x[4]*x[5]*d[3]^2*d[4]*d[5] 
//+(1)/(-2) *x[3]*x[5]*d[4]^3*d[5] +x[5]^2*d[3]^2*d[5]^2 +(-1)/(-9) *x[3]*d[2]^3 
//+(-1)/(9) *x[3]*d[3]^3 +(2)/(-9) *x[4]*d[3]^2*d[4] +(-7)/(9) *x[3]*d[4]^3 
//+(4)/(3) *x[5]*d[3]^2*d[5] +(2)/(-9) *d[3]^2,
//(-9)/(2) *x[1]*x[2]*x[5]^3*d[5]^4 +(-18)*x[1]*x[2]*x[5]^2*d[5]^3 
//+(10)/(-1) *x[1]*x[2]*x[5]*d[5]^2 +(-1)/(-9) *x[3]*x[4]*d[2]^2*d[3]*d[4] 
//+(-1)/(-9) *x[2]*x[4]*d[3]^3*d[4] +(-1)/(-9) *x[2]*x[3]*d[3]*d[4]^3 
//+(1)/(-6) *x[3]*x[5]*d[2]^2*d[3]*d[5] +(1)/(-6) *x[2]*x[5]*d[3]^3*d[5] 
//+(1)/(-6) *x[4]*x[5]*d[2]^2*d[4]*d[5] +(-1)/(6) *x[2]*x[5]*d[4]^3*d[5] 
//+(-1)/(-2) *x[5]^2*d[2]^2*d[5]^2 +(1)/(9) *x[3]*d[2]^2*d[3] +(-1)/(-9) *x[2]*d[3]^3 
//+(-1)/(-9) *x[4]*d[2]^2*d[4] +(-1)/(-9) *x[2]*d[4]^3 +(1)/(3) *x[5]*d[2]^2*d[5] 
//+(-1)/(-9) *d[2]^2,
//(-9)/(2) *x[1]*x[5]^4*d[5]^5 +(-36)*x[1]*x[5]^3*d[5]^4 +(64)/(-1) *x[1]*x[5]^2*d[5]^3 
//+(-20)*x[1]*x[5]*d[5]^2 +(-1)/(9) *x[3]*x[4]*d[2]^3*d[3]*d[4] 
//+(1)/(9) *x[3]*x[4]*d[3]^4*d[4] +(-2)/(-9) *x[3]^2*d[3]^2*d[4]^3 
//+(1)/(9) *x[3]*x[4]*d[3]*d[4]^4 +(-1)/(-9) *x[3]*x[5]*d[2]^3*d[3]*d[5] 
//+(1)/(-9) *x[3]*x[5]*d[3]^4*d[5] +(-1)/(-9) *x[4]*x[5]*d[2]^3*d[4]*d[5] 
//+(1)/(3) *x[4]*x[5]*d[3]^3*d[4]*d[5] +(1)/(3) *x[3]*x[5]*d[3]*d[4]^3*d[5] 
//+(1)/(-9) *x[4]*x[5]*d[4]^4*d[5] +(-1)/(6) *x[5]^2*d[2]^3*d[5]^2 
//+(-1)/(2) *x[5]^2*d[3]^3*d[5]^2 +(-1)/(2) *x[5]^2*d[4]^3*d[5]^2 
//+(-1)/(9) *x[3]*d[2]^3*d[3] +(1)/(9) *x[3]*d[3]^4 +(-1)/(9) *x[4]*d[2]^3*d[4] 
//+(-1)/(-3) *x[4]*d[3]^3*d[4] +(11)/(9) *x[3]*d[3]*d[4]^3 +(-1)/(-9) *x[4]*d[4]^4 
//+(-2)/(3) *x[5]*d[3]^3*d[5] +(2)/(-3) *x[5]*d[4]^3*d[5] +(-1)/(9) *d[2]^3 
//+(-1)/(-3) *d[3]^3 +(7)/(9) *d[4]^3,
//(-2)*x[0]*x[4]*d[3] +(6)/(-1) *x[3]^2*x[4]*d[5],
//-x[0]*d[3] +(-3)*x[3]^2*d[5],
//(-2)*x[0]*x[4]*d[2] +(6)/(-1) *x[2]^2*x[4]*d[5],
//-x[0]*d[2] +(-3)*x[2]^2*d[5],
//(-6)*x[0]*x[1]*x[4] +(6)*x[4],
//-x[0]*x[1] +(-1)/(-1) ,
//-x[1]*x[2]^3 +-x[1]*x[3]^3 +-x[1]*x[4]^3 +x[5]]

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-RingWeyl1.C,v 1.15 2022/02/11 09:53:59 abbott Exp $
// $Log: test-RingWeyl1.C,v $
// Revision 1.15  2022/02/11 09:53:59  abbott
// Summary: Updated copyright (redmine 855)
//
// Revision 1.14  2018/05/25 09:28:28  abbott
// Summary: Added const (just to be cleaner)
//
// Revision 1.13  2017/04/21 20:27:58  bigatti
// -- fixed GB interreduction --> changed output for tidyGens
//
// Revision 1.12  2016/09/16 16:36:40  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.11  2015/05/11 13:00:52  abbott
// Summary: Put all code into namespace CoCoA (redmine #642)
// Author: JAA
//
// Revision 1.10  2012/02/10 10:40:58  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.9  2012/02/08 17:34:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.8  2011/03/08 18:00:36  bigatti
// -- changed size_t into long
//
// Revision 1.7  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.6  2010/05/14 11:12:22  bigatti
// -- fixed output
//
// Revision 1.5  2008/11/18 17:25:53  bigatti
// -- added small test for TidyGens
//
// Revision 1.4  2007/10/30 17:14:05  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/31 16:01:45  bigatti
// -- default implementation for IamField, myCharacteristic  in PolyRing
// -- added !IsCommutative in test
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.6  2007/03/08 12:18:26  bigatti
// -- added #include "CoCoA/symbol.H"
//
// Revision 1.5  2007/03/03 14:13:21  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.3  2007/02/26 17:08:12  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1  2006/11/14 17:47:59  cocoa
// -- added first test for RingWeyl
//
// Revision 1.4  2006/08/30 15:22:12  cocoa
// -- added example from WeylAlgebra9.C
//
// Revision 1.3  2006/08/17 10:03:24  cocoa
// -- creation of RingWeyl with given list of symbols
// -- updated header
//
