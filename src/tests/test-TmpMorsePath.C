//   Copyright (c)  2013-2015  John Abbott,  Anna M. Bigatti
//   Original author: 2013-2015  Mario Albert

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
#include "CoCoA/TmpMorsePaths.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/error.H"

using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

DynamicBitset VectorBoolToDynamicBitset(const std::vector<bool>& bools)
{
  DynamicBitset result(bools.size());
  long counter(0);
  for (std::vector<bool>::const_iterator i = bools.begin(); i != bools.end(); ++i)
  {
    result.mySet(counter, *i);
    ++counter;
  }
  return result;
}

void test_myAddPath()
{
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[0]);
  PolInput.push_back(x[1]);
  PolInput.push_back(x[2]);
  PolInput.push_back(x[3]);
  PolInput.push_back(x[4]);
  JBMill mill(JBMill::Builder().setInput(PolInput));

  std::vector<MorseElement::JBElem> TestBasis;
  TestBasis.push_back(MorseElement::JBElem(x[4], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[4]))), 1));
  TestBasis.push_back(MorseElement::JBElem(x[3], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[3]))), 2));
  TestBasis.push_back(MorseElement::JBElem(x[2], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[2]))), 3));
  TestBasis.push_back(MorseElement::JBElem(x[1], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[1]))), 4));
  TestBasis.push_back(MorseElement::JBElem(x[0], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[0]))), 5));

  std::map<MorseElement, MorsePaths> TestResolution;
  for (std::vector<MorseElement::JBElem>::iterator i = TestBasis.begin(); i != TestBasis.end(); ++i)
  {
    TestResolution.insert(std::pair<MorseElement, MorsePaths>(MorseElement(DynamicBitset(5), i), MorsePaths()));
  }
  MorsePaths test;
  ConstResIter iter(TestResolution.begin());
  test.myAddPath(iter, x[0]);
  CoCoA_ASSERT_ALWAYS(x[0] == test.myGetPath(iter));
  ++iter;
  test.myAddPath(iter, x[1]);

  CoCoA_ASSERT_ALWAYS(x[1] == test.myGetPath(iter));

  --iter;
  test.myAddPath(iter, -x[0]);
  PathMap paths(test.myGetPaths());
  CoCoA_ASSERT_ALWAYS(paths.size() == 1);
  ++iter;
  CoCoA_ASSERT_ALWAYS(paths[iter] == x[1]);
}

void program()
{
  GlobalManager CoCoAFoundations;
  test_myAddPath();

}


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
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

  BuildInfo::PrintAll(cerr);
  return 1;
}
