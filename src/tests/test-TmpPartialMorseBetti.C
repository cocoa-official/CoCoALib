//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Original author: 2015  Mario Albert

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
#include "CoCoA/TmpJBMill.H"
#include "CoCoA/TmpPartialMorseBetti.H"
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/error.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/library.H"

using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

const bool WRITE_VECS = false;
const bool WRITE_NAME = false;

struct LPPCompare
{
  bool operator()(const RingElem& r1, const RingElem& r2)
  {
    return LPP(r1) < LPP(r2);
  }
};

void test_one(const std::vector<RingElem>& res)
{
  CoCoA_ASSERT_ALWAYS(res.size() == 1);
  CoCoA_ASSERT_ALWAYS(IsOne(res[0]));
}

void test_principal_ideal(const std::vector<RingElem>& to_test, RingElem result )
{
  CoCoA_ASSERT_ALWAYS(to_test.size() == 1);
  CoCoA_ASSERT_ALWAYS(to_test[0] == result);
}

void test_equal_vectors(std::vector<RingElem>& vec1, std::vector<RingElem>& vec2)
{
  std::sort(vec1.begin(), vec1.end(), LPPCompare());
  std::sort(vec2.begin(), vec2.end(), LPPCompare());
  SparsePolyRing ring(owner(vec1.front()));
  for (vector<RingElem>::iterator i = vec1.begin(); i != vec1.end(); ++i) {
    *i = (*i) / monomial(ring, LC(*i), LPP(*i)/LPP(*i));
  }
  for (vector<RingElem>::iterator i = vec2.begin(); i != vec2.end(); ++i) {
    *i = (*i) / monomial(ring, LC(*i), LPP(*i)/LPP(*i));
  }
  if (WRITE_VECS) {
    std::cout << "vec1 = ";
    for(std::vector<RingElem>::const_iterator i = vec1.begin(); i != vec1.end(); ++i)
    {
      std::cout << *i << "; ";
    }
    std::cout << std::endl;
    std::cout << "vec2 = ";
    for(std::vector<RingElem>::const_iterator i = vec2.begin(); i != vec2.end(); ++i)
    {
      std::cout << *i << "; ";
    }
    std::cout << std::endl;
  }
  CoCoA_ASSERT_ALWAYS(std::equal(vec1.begin(), vec1.end(), vec2.begin()));
}

void test_myPossibleWedgesOfLength()
{
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[0]);
  PolInput.push_back(x[1]);
  PolInput.push_back(x[2]);
  PolInput.push_back(x[3]);
  PolInput.push_back(x[4]);
  JBMill mill(JBMill::Builder().setInput(PolInput));


  PartialMorseBetti mg(mill);
  std::vector<DynamicBitset> TestResult;
  std::vector<long> input;

  TestResult = mg.myPossibleWedgesOfLength(input, 0);
  CoCoA_ASSERT_ALWAYS(TestResult.size() == 1);
  CoCoA_ASSERT_ALWAYS(TestResult[0].IamAll0s());


  input.push_back(0);
  input.push_back(1);
  input.push_back(2);
  TestResult = mg.myPossibleWedgesOfLength(input, 2);
  CoCoA_ASSERT_ALWAYS(TestResult.size() == 3);

  CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(0) == true);
  CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(1) == true);
  CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(2) == false);

  CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(0) == true);
  CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(1) == false);
  CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(2) == true);

  CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(0) == false);
  CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(1) == true);
  CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(2) == true);



  TestResult = mg.myPossibleWedgesOfLength(input, 4);
  CoCoA_ASSERT_ALWAYS(TestResult.size() == 0);
}

void test_myComputeGeneralPartialBasis()
{
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,1), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[0]);
  PolInput.push_back(x[1]);
  JBMill mill(JBMill::Builder().setInput(PolInput));

  PartialMorseBetti mg(mill);

  std::vector<MorseElement> elems(mg.myComputeGeneralColumnBasis(0,0));

  CoCoA_ASSERT_ALWAYS(elems.size() == 3);
  std::vector<MorseElement::JBElem > TestBasis;
  TestBasis.push_back(MorseElement::JBElem(x[0], DynamicBitset(LPP(x[0] * x[1])), 0));
  TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1])), 1));
  TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1])), 1));
  std::vector<MorseElement::JBElem >::iterator TBIter(TestBasis.begin());
  MorseElement m1(DynamicBitset(2), TBIter);
  // std::cout <<  "m1 = " << m1 << std::endl;
  ++TBIter;
  MorseElement m2(DynamicBitset(2), TBIter);
  // std::cout <<  "m2 = " << m2 << std::endl;
  ++TBIter;
  MorseElement m3(DynamicBitset(LPP(x[0])), TBIter);
  // std::cout <<  "m3 = " << m3 << std::endl;
  CoCoA_ASSERT_ALWAYS(elems[0] == m2);
  CoCoA_ASSERT_ALWAYS(elems[1] == m3);
  CoCoA_ASSERT_ALWAYS(elems[2] == m1);

  elems = mg.myComputeGeneralColumnBasis(2, 0);

  CoCoA_ASSERT_ALWAYS(elems.size() == 1);
  CoCoA_ASSERT_ALWAYS(elems[0]     == m3);

  elems = mg.myComputeGeneralColumnBasis(0, 1);

  CoCoA_ASSERT_ALWAYS(elems.size() == 1);
  CoCoA_ASSERT_ALWAYS(elems[0]     == m3);

  elems = mg.myComputeGeneralColumnBasis(2, 1);

  CoCoA_ASSERT_ALWAYS(elems.size() == 1);
  CoCoA_ASSERT_ALWAYS(elems[0]     == m3);

  elems = mg.myComputeGeneralColumnBasis(3, 1);

  CoCoA_ASSERT_ALWAYS(elems.size() == 0);
}

void test_BinomialIdealCycle()
{
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  vector<symbol> Xs(SymbolRange("x", 0, 3));
  vector<symbol> Ys(SymbolRange("y", 0, 3));
  vector<symbol> AllVars;
  AllVars.insert(AllVars.end(), Xs.begin(), Xs.end());
  AllVars.insert(AllVars.end(), Ys.begin(), Ys.end());
  PolyRing P = NewPolyRing(Q, AllVars);
  vector<RingElem> gens;
  // std::cout <<  "(RingElem(x[5]) > RingElem(x[4])) = " << (LPP(RingElem(P,"x[5]")) > LPP(RingElem(P,"x[4]"))) << std::endl;
  gens.push_back(RingElem(P, "x[0]*y[1]-x[1]*y[0]"));
  gens.push_back(RingElem(P, "x[1]*y[2]-x[2]*y[1]"));
  gens.push_back(RingElem(P, "x[2]*y[3]-x[3]*y[2]"));
  gens.push_back(RingElem(P, "x[0]*y[3]-x[3]*y[0]"));
  JBMill mill(JBMill::Builder().setInput(gens));
  // std::cout <<  "mill.myReturnJB() = " << mill.myReturnJB() << std::endl;
  PartialMorseBetti mg(mill);
  matrix Bettis(NewDenseMat(RingZZ(), 1, 1));
  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,0));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,1));

  SetEntry(Bettis, 0, 0, 9);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,2));

  SetEntry(Bettis, 0, 0, 8);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,3));

  SetEntry(Bettis, 0, 0, 2);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,4));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,5));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,6));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,7));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(2,8));


  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,0));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,1));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,2));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,3));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,4));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,5));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,6));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,7));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(3,8));


  matrix TwoBettis(NewDenseMat(RingZZ(), 2, 1));
  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(1,0));

  SetEntry(Bettis, 0, 0, 4);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(1,1));

  SetEntry(TwoBettis, 0, 0, 0);
  SetEntry(TwoBettis, 1, 0, 9);
  CoCoA_ASSERT_ALWAYS(TwoBettis == mg.myComputeDownToBettiNumber(1,2));

  SetEntry(TwoBettis, 0, 0, 0);
  SetEntry(TwoBettis, 1, 0, 8);
  CoCoA_ASSERT_ALWAYS(TwoBettis == mg.myComputeDownToBettiNumber(1,3));

  SetEntry(TwoBettis, 0, 0, 0);
  SetEntry(TwoBettis, 1, 0, 2);
  CoCoA_ASSERT_ALWAYS(TwoBettis == mg.myComputeDownToBettiNumber(1,4));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(1,5));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(1,6));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(1,7));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(1,8));

  matrix ThreeBettis(NewDenseMat(RingZZ(), 3, 1));
  SetEntry(Bettis, 0, 0, 1);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(0,0));

  SetEntry(TwoBettis, 0, 0, 0);
  SetEntry(TwoBettis, 1, 0, 4);
  CoCoA_ASSERT_ALWAYS(TwoBettis == mg.myComputeDownToBettiNumber(0,1));

  SetEntry(ThreeBettis, 0, 0, 0);
  SetEntry(ThreeBettis, 1, 0, 0);
  SetEntry(ThreeBettis, 2, 0, 9);
  CoCoA_ASSERT_ALWAYS(ThreeBettis == mg.myComputeDownToBettiNumber(0,2));

  SetEntry(ThreeBettis, 0, 0, 0);
  SetEntry(ThreeBettis, 1, 0, 0);
  SetEntry(ThreeBettis, 2, 0, 8);
  CoCoA_ASSERT_ALWAYS(ThreeBettis == mg.myComputeDownToBettiNumber(0,3));

  SetEntry(ThreeBettis, 0, 0, 0);
  SetEntry(ThreeBettis, 1, 0, 0);
  SetEntry(ThreeBettis, 2, 0, 2);
  CoCoA_ASSERT_ALWAYS(ThreeBettis == mg.myComputeDownToBettiNumber(0,4));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(0,5));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(0,6));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(0,7));

  SetEntry(Bettis, 0, 0, 0);
  CoCoA_ASSERT_ALWAYS(Bettis == mg.myComputeDownToBettiNumber(0,8));
}

void test_BinomialIdealCycleSingleBetti()
{
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  vector<symbol> Xs(SymbolRange("x", 0, 3));
  vector<symbol> Ys(SymbolRange("y", 0, 3));
  vector<symbol> AllVars;
  AllVars.insert(AllVars.end(), Xs.begin(), Xs.end());
  AllVars.insert(AllVars.end(), Ys.begin(), Ys.end());
  PolyRing P = NewPolyRing(Q, AllVars);
  vector<RingElem> gens;
  // std::cout <<  "(RingElem(x[5]) > RingElem(x[4])) = " << (LPP(RingElem(P,"x[5]")) > LPP(RingElem(P,"x[4]"))) << std::endl;
  gens.push_back(RingElem(P, "x[0]*y[1]-x[1]*y[0]"));
  gens.push_back(RingElem(P, "x[1]*y[2]-x[2]*y[1]"));
  gens.push_back(RingElem(P, "x[2]*y[3]-x[3]*y[2]"));
  gens.push_back(RingElem(P, "x[0]*y[3]-x[3]*y[0]"));
  JBMill mill(JBMill::Builder().setInput(gens));
  // std::cout <<  "mill.myReturnJB() = " << mill.myReturnJB() << std::endl;
  PartialMorseBetti mg(mill);

  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 0));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 1));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 2));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 3));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 4));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 5));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 6));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 7));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 3, 8));

  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 2, 0));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 2, 1));
  CoCoA_ASSERT_ALWAYS(9 == BettiNumber(mill, 2, 2));
  CoCoA_ASSERT_ALWAYS(8 == BettiNumber(mill, 2, 3));
  CoCoA_ASSERT_ALWAYS(2 == BettiNumber(mill, 2, 4));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 2, 5));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 2, 6));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 2, 7));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 2, 8));

  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 0));
  CoCoA_ASSERT_ALWAYS(4 == BettiNumber(mill, 1, 1));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 2));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 3));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 4));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 5));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 6));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 7));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 1, 8));

  CoCoA_ASSERT_ALWAYS(1 == BettiNumber(mill, 0, 0));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 1));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 2));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 3));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 4));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 5));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 6));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 7));
  CoCoA_ASSERT_ALWAYS(0 == BettiNumber(mill, 0, 8));
}

void program()
{
  GlobalManager CoCoAFoundations;
  test_myPossibleWedgesOfLength();
  test_myComputeGeneralPartialBasis();
  test_BinomialIdealCycle();
  test_BinomialIdealCycleSingleBetti();
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
