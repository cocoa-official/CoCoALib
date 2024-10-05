//   Copyright (c)  2010  John Abbott,  Anna Bigatti

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


#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows the use of DynamicBitsets.  \n";

const string LongDescription =
  "We show how to convert to and from PPMonoidElem.             \n"
  "We list the available operations which are a lot faster than \n"
  "the corresponding ones on PPMonoidElems.                     \n"
  "By default DynamicBitsets are printed as STL bitsets, i.e.   \n"
  "0th bit in the rightmost position.  We show other options.   \n"
  "\n"
  "Note that the ordering is xel, i.e. like lex but starting    \n"
  "the comparison from the position with highest index (as they \n"
  "are printed by default)";

//----------------------------------------------------------------------

namespace CoCoA
{

  void ExDynamicBitset(long n)
  {
    cout << "DynamicBitset of length " << n << endl << endl;

    PPMonoid PPM = NewPPMonoid(SymbolRange("x", 0, n-1), lex);
    const vector<PPMonoidElem>& x = indets(PPM);
    DynamicBitset b1(product(indets(PPM)));
    DynamicBitset b2(x[0] * x[n-1]);
    DynamicBitset b3(x[0] * x[1]);
    DynamicBitset b(n);  // of length n, full of 0s

    cout << "-- Styles for printing: default is \"clean\"" << endl;
    cout << "-- Note that 0th component is last (as for bitset)" << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::clean;
    cout << "ourOutputStyle = clean          b3 = " << b3 << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::WithSeparators;
    cout << "ourOutputStyle = WithSeparators b3 = " << b3 << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::AsRevVecOfLong;
    cout << "ourOutputStyle = AsRevVecOfLong b3 = " << b3 << endl;
    cout << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::WithSeparators;

    cout << "ourOutputStyle = WithSeparators --> " << endl;
    cout << " b1 = " << b1 << endl;
    cout << " b2 = " << b2 << endl;
    cout << " b3 = " << b3 << endl;
    cout << endl;

    cout << "NewPP(PPM, b2) = " << NewPP(PPM, b2) << endl;
    cout << endl;

    cout << "b1.IamAll1s() = " << b1.IamAll1s() << endl;
    cout << "b2.IamAll0s() = " << b2.IamAll0s() << endl;
    cout << "count(b1) = " << count(b1) << endl;
    cout << "count(b2) = " << count(b2) << endl;
    cout << "count(b3) = " << count(b3) << endl;
    cout << endl;

    cout << "-- ordered by xel: Oth position is smallest" << endl;
    cout << "-- (like lex, as printed)" << endl;
    cout << "b1 < b2 = " << (b1 < b2) << endl;
    cout << "b2 < b3 = " << (b2 < b3) << endl;
    cout << endl;

    cout << "-- functions checking length compatibility:" << endl;
    cout << "b2 | b3   gives  " << (b2 | b3) << endl;
    cout << "b2 & b3   gives  " << (b2 & b3) << endl;
    cout << "b1 - b3   gives  " << (b1 - b3) << endl;
    cout << "IsSubset(b3, b1)   gives  " << IsSubset(b3, b1) << endl;
    cout << "IsDisjoint(b3, b1) gives  " << IsDisjoint(b3, b1) << endl;
    cout << "Is1At(b3, 1)       gives  " << Is1At(b3, 1)  << endl;
    cout << endl;

    cout << "-- functions not checking length compatibility:" << endl;
    b = b2; b |= b3;
    cout << "b = b2; b |= b3;   gives  " << b << endl;
    b = b2; b &= b3;
    cout << "b = b2; b &= b3;   gives  " << b << endl;
    b = b1; b -= b3;
    cout << "b = b1; b -= b3;   gives  " << b << endl;
    cout << "b.mySet(0)         gives  " << b.mySet(0)   << endl;
    cout << "b.mySet(2,0)       gives  " << b.mySet(2,0) << endl;
    cout << "b.mySet(2,1)       gives  " << b.mySet(2,1) << endl;
    cout << "b3.IamSubset(b1)   gives  " << b3.IamSubset(b1) << endl;
    cout << "b3.IamDisjoint(b1) gives  " << b3.IamDisjoint(b1) << endl;
    cout << "b3.Iam1At(1)       gives  " << b3.Iam1At(1)  << endl;
    cout << endl;

    cout << "------------------------------------------------" << endl << endl;
  }


//-- program --------------------------------------------------------------
// we run TestPolyRing on predefined and user-defined orderings

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << "DynamicBitset::ourNumBitsInBlock = "
         <<  DynamicBitset::ourNumBitsInBlock << endl << endl;
  
    ExDynamicBitset(10);
    ExDynamicBitset(66);
    ExDynamicBitset(201);
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
