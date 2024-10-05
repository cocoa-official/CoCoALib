//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingElemInput.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// This test checks whether ReadExpr works correctly for integers,
// rationals and "decimals" (i.e. rationals written with a decimal point)
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void ReadExprShouldWork(const string& s, const BigRat& q)
  {
    CoCoA_ASSERT_ALWAYS(RingElem(RingQQ(), s) == q);
    CoCoA_ASSERT_ALWAYS(RingElem(RingQQ(), s+" ") == q);
    CoCoA_ASSERT_ALWAYS(RingElem(RingQQ(), " "+s) == q);
    CoCoA_ASSERT_ALWAYS(RingElem(RingQQ(), " "+s+" ") == q);
  }

  void ReadExprShouldFail(const string& s)
  {
    try
    {
      RingElem(RingQQ(), s);
      throw "ReadExprShouldFail: should never get here";
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      // ignore CoCoA error
    }
  }

  void RingElemVecShouldWork(const string& s, const vector<RingElem>& v)
  {
    CoCoA_ASSERT_ALWAYS(!v.empty());
    const ring& R = owner(v[0]);
    CoCoA_ASSERT_ALWAYS(RingElemVec(R, s) == v);
    CoCoA_ASSERT_ALWAYS(RingElemVec(R, s+" ") == v);
    CoCoA_ASSERT_ALWAYS(RingElemVec(R, " "+s) == v);
    CoCoA_ASSERT_ALWAYS(RingElemVec(R, " "+s+" ") == v);
  }

  void RingElemVecShouldFail(const ring& R, const string& s)
  {
    try
    {
      RingElemVec(R, s);
      throw "ReadExprShouldFail: should never get here";
    }
    catch (const CoCoA::ErrorInfo& err) { /* ignore CoCoA error*/ }
  }


  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    // integer values
    ReadExprShouldWork("0", BigRat(0));
    ReadExprShouldWork("007", BigRat(7));
    ReadExprShouldWork("1", BigRat(1));
    ReadExprShouldWork("100000000000000000000000000000000000000000000000000", power(BigRat(10),50));
    ReadExprShouldWork("10^50", power(BigRat(10),50));
    ReadExprShouldWork("11739085287969531650666649599035831993898213898723001", power(BigRat(11),50));
    ReadExprShouldWork("+0", BigRat(0));
    ReadExprShouldWork("+007", BigRat(7));
    ReadExprShouldWork("+1", BigRat(1));
    ReadExprShouldWork("+100000000000000000000000000000000000000000000000000", power(BigRat(10),50));
    ReadExprShouldWork("+10^50", power(BigRat(10),50));
    ReadExprShouldWork("+11739085287969531650666649599035831993898213898723001", power(BigRat(11),50));
    ReadExprShouldWork("-0", BigRat(0));
    ReadExprShouldWork("-007", BigRat(-7));
    ReadExprShouldWork("-1", BigRat(-1));
    ReadExprShouldWork("-100000000000000000000000000000000000000000000000000", -power(BigRat(10),50));
    ReadExprShouldWork("-10^50", -power(BigRat(10),50));
    ReadExprShouldWork("-11739085287969531650666649599035831993898213898723001", -power(BigRat(11),50));
    ReadExprShouldWork("- 0", BigRat(0));
    ReadExprShouldWork("- 007", BigRat(-7));
    ReadExprShouldWork("- 1", BigRat(-1));
    ReadExprShouldWork("- 100000000000000000000000000000000000000000000000000", -power(BigRat(10),50));
    ReadExprShouldWork("- 10^50", -power(BigRat(10),50));
    ReadExprShouldWork("- 11739085287969531650666649599035831993898213898723001", -power(BigRat(11),50));


    // Rationals as num/den
    ReadExprShouldWork("0/1", BigRat(0));    
    ReadExprShouldWork("1/1", BigRat(1));    
    ReadExprShouldWork("-1/1", BigRat(-1));    
    ReadExprShouldWork("2/2", BigRat(1));    
    ReadExprShouldWork("200/300", BigRat(2,3));    
    ReadExprShouldWork("1/2/3", BigRat(1,6));    

    ReadExprShouldWork("0/ 1", BigRat(0));    
    ReadExprShouldWork("1/ 1", BigRat(1));    
    ReadExprShouldWork("-1/ 1", BigRat(-1));    
    ReadExprShouldWork("2/ 2", BigRat(1));    
    ReadExprShouldWork("200/ 300", BigRat(2,3));    
    ReadExprShouldWork("1/ 2/ 3", BigRat(1,6));    

    ReadExprShouldWork("0 / 1", BigRat(0));    
    ReadExprShouldWork("1 / 1", BigRat(1));    
    ReadExprShouldWork("-1 / 1", BigRat(-1));    
    ReadExprShouldWork("2 / 2", BigRat(1));    
    ReadExprShouldWork("200 / 300", BigRat(2,3));    
    ReadExprShouldWork("1 / 2 / 3", BigRat(1,6));    

    ReadExprShouldWork("5/4^3", BigRat(5, 64));
    ReadExprShouldWork("5/4 ^ 3", BigRat(5, 64));


    // Rationals expressed as decimals

    ReadExprShouldWork("0.", BigRat(0));
    ReadExprShouldWork("0.0", BigRat(0));
    ReadExprShouldWork("0.5", BigRat(1,2));
    ReadExprShouldWork("1.5", BigRat(3,2));
    ReadExprShouldWork("1.500", BigRat(3,2));
    ReadExprShouldWork("1.001", BigRat(1001,1000));

    ReadExprShouldWork("-0.", BigRat(0));
    ReadExprShouldWork("-0.0", BigRat(0));
    ReadExprShouldWork("-0.5", BigRat(-1,2));
    ReadExprShouldWork("-1.5", BigRat(-3,2));
    ReadExprShouldWork("-1.500", BigRat(-3,2));
    ReadExprShouldWork("-1.001", BigRat(-1001,1000));

    ReadExprShouldWork("- 0.", BigRat(0));
    ReadExprShouldWork("- 0.0", BigRat(0));
    ReadExprShouldWork("- 0.5", BigRat(-1,2));
    ReadExprShouldWork("- 1.5", BigRat(-3,2));
    ReadExprShouldWork("- 1.500", BigRat(-3,2));
    ReadExprShouldWork("- 1.001", BigRat(-1001,1000));

    ReadExprShouldWork("1.25^3", BigRat(125,64)); // NB different from "5/4^3"!
    ReadExprShouldWork("1.25 ^ 3", BigRat(125,64)); // NB different from "5/4^3"!

    ReadExprShouldWork("2^3", BigRat(8));
    ReadExprShouldWork("2 ^ 3", BigRat(8));
    ReadExprShouldWork("2^(3)", BigRat(8));
    ReadExprShouldWork("2 ^ ( 3 )", BigRat(8));
    ReadExprShouldWork("2^(+3)", BigRat(8));
    ReadExprShouldWork("2 ^ ( + 3 )", BigRat(8));
    ReadExprShouldWork("2^(-3)", BigRat(1,8));
    ReadExprShouldWork("2 ^ ( -3 )", BigRat(1,8));
    ReadExprShouldWork("-2^2", BigRat(-4));
    ReadExprShouldWork("-2^(-2)", BigRat(-1,4));

    ReadExprShouldFail("");
    ReadExprShouldFail("1 2");
    ReadExprShouldFail("0/-1");
    ReadExprShouldFail("0 .");
    ReadExprShouldFail(".5");
    ReadExprShouldFail("0 .5");
    ReadExprShouldFail("0. 5");
    ReadExprShouldFail("1..2");
    ReadExprShouldFail("1.2.3");
    ReadExprShouldFail("--1");
    ReadExprShouldFail("++1");
    ReadExprShouldFail("+-1");
    ReadExprShouldFail("-+1");

    ReadExprShouldFail("1^-1");
    ReadExprShouldFail("1^((-1))");

    //--------------------------------------------
    // Tests for RingElemVec
    ring P = NewPolyRing(RingQQ(), symbols("x,y[1]"));
    const RingElem x = RingElem(P, "x");
    const RingElem y1 = RingElem(P, "y[1]");
    vector<RingElem> v; v.push_back(x); v.push_back(x+y1);

    CoCoA_ASSERT_ALWAYS(RingElemVec(P, "[]").empty());
    CoCoA_ASSERT_ALWAYS(RingElemVec(P, " [ ] ").empty());
    
    RingElemVecShouldWork("[x,x+y[1]]", v );
    RingElemVecShouldWork("[ x , x + y[ 1 ] ]", v );

    RingElemVecShouldFail(P, "");
    RingElemVecShouldFail(P, "[");
    RingElemVecShouldFail(P, "]");
    RingElemVecShouldFail(P, "[],"); // junk after last char
    RingElemVecShouldFail(P, "x");
    RingElemVecShouldFail(P, "[x");
    RingElemVecShouldFail(P, "[x+");
    RingElemVecShouldFail(P, "[y]");
    RingElemVecShouldFail(P, "[y[1]");
    RingElemVecShouldFail(P, "[x,,x]");
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
