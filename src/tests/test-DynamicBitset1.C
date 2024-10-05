//   Copyright (c)  2010 Anna Bigatti

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
#include "CoCoA/DynamicBitset.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/VectorOps.H"  // for product
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
#include "CoCoA/time.H"


#include<iostream>
using std::cout;
using std::cerr;
using std::endl;
#include<vector>
using std::vector;

//----------------------------------------------------------------------
// Test for DynamicBitset
// functions: 
//----------------------------------------------------------------------

namespace CoCoA
{

  vector<long> VecCtor(long a1, long a2, long a3, long a4)
  {
    vector<long> v;
    v.push_back(a1); v.push_back(a2); v.push_back(a3); v.push_back(a4);
    return v;
  }


//-- TestDynamicBitset ------------------------------------------------------
// behaviour on different lengths

  void TestDynamicBitset(long n)
  {
//  cout << "TESTING: " << n << endl << endl;

    PPMonoid PPM = NewPPMonoid(SymbolRange("x", 0, n-1), lex);
    PPMonoidElem t1 = product(indets(PPM));
    PPMonoidElem t2 = indet(PPM,0) * indet(PPM,n-1);
    PPMonoidElem t3 = indet(PPM,0) * indet(PPM,1);
    DynamicBitset b1(t1);
    DynamicBitset b2(t2);
    DynamicBitset b3(t3);    //  exponents(v, t3);

//   cout << "Given b1 = " << b1 << endl;
//   cout << "  and b2 = " << b2 << endl;
//   cout << "  and b3 = " << b3 << endl << endl;

    CoCoA_ASSERT_ALWAYS(NewPP(PPM, b2) == t2);
    CoCoA_ASSERT_ALWAYS(NewPP(PPM, b3) == t3);

    CoCoA_ASSERT_ALWAYS(b1.IamAll1s());
    CoCoA_ASSERT_ALWAYS(!b2.IamAll1s());
    CoCoA_ASSERT_ALWAYS(!b3.IamAll1s());

    CoCoA_ASSERT_ALWAYS(b2 < b1);
    CoCoA_ASSERT_ALWAYS(b3 < b2);
  
//   cout << "-- not checking length compatibility:" << endl;
//   cout << "b2 | b3   gives  " << (b2 | b3) << endl;
//   cout << "b2 & b3   gives  " << (b2 & b3) << endl;
//   cout << "b1 - b2   gives  " << (b1 - b2) << endl;
//   cout << "-- checking length compatibility:" << endl;
//   cout << "union(b2,b3)   gives  " << union(b2,b3) << endl;
//   cout << "intersection(b2,b3)   gives  " << intersection(b2,b3) << endl;
//   cout << "difference(b1,b2)   gives  " << difference(b1,b2) << endl;
//   cout << "IsSubset(b3, b1)   gives  " << IsSubset(b3, b1) << endl;
//   cout << "contains(b1, b3)   gives  " << contains(b1, b3) << endl;

//   PPMonoidElem t1t2 = t1*t2;
//   PPMonoidElem t2t1 = t2*t1;
//   CoCoA_ASSERT_ALWAYS(t1t2 == t2t1);

    CoCoA_ASSERT_ALWAYS( !b2.IamAll0s() );
    CoCoA_ASSERT_ALWAYS( Is1At(b1, 2) );
  }


  //-- program --------------------------------------------------------------
  // we run TestPolyRing on predefined and user-defined orderings

  void program()
  {
    GlobalManager CoCoAFoundations;

    TestDynamicBitset(4);
    TestDynamicBitset(32);
    TestDynamicBitset(64);
    TestDynamicBitset(64*5 - 1);
    TestDynamicBitset(4001);
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
