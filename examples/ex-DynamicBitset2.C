//   Copyright (c)  2014  John Abbott, Anna Bigatti

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
  "We show how to convert to and from RingElem.             \n"
  "We create an ideal I and convert its generator into DynamicBitsets.\n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void ExDynamicBitset(long n)
  {
    cout << "DynamicBitset of length " << n << endl << endl;

    PolyRing P = NewPolyRing(RingQQ(), SymbolRange("x", 0, n-1), lex);
    const vector<RingElem>& x = indets(P);
    ideal I = x[0] * x[2] * ideal(x);  // NB not square-free

    const vector<RingElem>& g = gens(I);
    vector<DynamicBitset> DBGens;
    for (long i=0; i<len(g); ++i)
      DBGens.push_back(DynamicBitset(LPP(g[i]))); // LPP(g[i]) is a PPMonoidElem

    cout << "Conversion of " << g << endl;
    cout << "-- Styles for printing: default is \"clean\"" << endl;
    cout << "-- Note that 0th component is last (as for bitset)" << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::clean;
    cout << "ourOutputStyle = clean" << endl << DBGens << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::WithSeparators;
    cout << "ourOutputStyle = WithSeparators" << endl << DBGens << endl;
    DynamicBitset::ourOutputStyle = DynamicBitset::OutputStyle::AsRevVecOfLong;
    cout << "ourOutputStyle = AsRevVecOfLong" << endl << DBGens << endl;
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
  
    ExDynamicBitset(6);
    ExDynamicBitset(10);
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
