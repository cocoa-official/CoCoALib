//   Copyright (c)  2013-2015  John Abbott,  Anna M. Bigatti
//   Original author:  2013-2015 Mario Albert

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
#include "CoCoA/TmpPBMill.H"
#include "CoCoA/TmpMorseElement.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/error.H"

using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;


std::vector<bool> DynamicBitsetToBool(const DynamicBitset& DynBitset)
{
  std::vector<bool> result;
  long length(len(DynBitset));
  for (long i = 0; i < length; ++i)
  {
    result.push_back(DynBitset.Iam1At(i));
  }
  return result;
}

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



void test_myDynamicBitsetToLong()
{
  DynamicBitset bitset(5);
  std::vector<long> output;

  output = myDynamicBitsetToLong(bitset);
  CoCoA_ASSERT_ALWAYS(len(output) == 0);

  bitset.mySet(0);
  output = myDynamicBitsetToLong(bitset);
  CoCoA_ASSERT_ALWAYS(len(output) == 1);
  CoCoA_ASSERT_ALWAYS(output[0] == 0);

  bitset.mySet(3);
  output = myDynamicBitsetToLong(bitset);
  CoCoA_ASSERT_ALWAYS(len(output) == 2);
  CoCoA_ASSERT_ALWAYS(output[0] == 0);
  CoCoA_ASSERT_ALWAYS(output[1] == 3);

  bitset.mySet(1);
  bitset.mySet(2);
  bitset.mySet(4);
  output = myDynamicBitsetToLong(bitset);
  CoCoA_ASSERT_ALWAYS(len(output) == 5);
  CoCoA_ASSERT_ALWAYS(output[0] == 0);
  CoCoA_ASSERT_ALWAYS(output[1] == 1);
  CoCoA_ASSERT_ALWAYS(output[2] == 2);
  CoCoA_ASSERT_ALWAYS(output[3] == 3);
  CoCoA_ASSERT_ALWAYS(output[4] == 4);
}


void test_myVectorLongToDynamicBitset()
{
  long ArrayLongs1[] = {2};
  std::vector<long> longs1 (ArrayLongs1, ArrayLongs1 + sizeof(ArrayLongs1) / sizeof(long) );
  DynamicBitset bitset = myVectorLongToDynamicBitset(longs1, 5);
  CoCoA_ASSERT_ALWAYS(len(bitset) == 5);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,0) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,1) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,2) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,3) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,4) == false);

  long ArrayLongs2[] = {0, 2};
  std::vector<long> longs2 (ArrayLongs2, ArrayLongs2 + sizeof(ArrayLongs2) / sizeof(long) );
  bitset = myVectorLongToDynamicBitset(longs2, 5);
  CoCoA_ASSERT_ALWAYS(len(bitset) == 5);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,0) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,1) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,2) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,3) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,4) == false);

  long ArrayLongs3[] = {0,1,2,3,4};
  std::vector<long> longs3 (ArrayLongs3, ArrayLongs3 + sizeof(ArrayLongs3) / sizeof(long) );
  bitset = myVectorLongToDynamicBitset(longs3, 5);
  CoCoA_ASSERT_ALWAYS(len(bitset) == 5);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,0) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,1) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,2) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,3) == true);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,4) == true);

  std::vector<long> longs4;
  bitset = myVectorLongToDynamicBitset(longs4, 5);
  CoCoA_ASSERT_ALWAYS(len(bitset) == 5);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,0) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,1) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,2) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,3) == false);
  CoCoA_ASSERT_ALWAYS(Is1At(bitset,4) == false);
}


void test_operators()
{
  ring QQ = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(QQ, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem>input;
  DynamicBitset wedge(5);
  PPMonoidElem RightFactor(one(PPM(PolyRing)));

  RingElem BasisElement1(x[1] * x[2] * x[0]);
  DynamicBitset ncrit1(5);
  long lexpos1(1);
  RingElem BasisElement2(x[1] * x[2] * x[0]);
  DynamicBitset ncrit2(5);
  long lexpos2(1);
  std::vector<MorseElement::JBElem> TestBasis;
  TestBasis.push_back(MorseElement::JBElem(BasisElement1, ncrit1, lexpos1));
  TestBasis.push_back(MorseElement::JBElem(BasisElement2, ncrit2, lexpos2));
  std::vector<MorseElement::JBElem>::iterator TBIter1(TestBasis.begin());
  std::vector<MorseElement::JBElem>::iterator TBIter2(TestBasis.begin());
  ++TBIter2;

  MorseElement m1(wedge, RightFactor, TBIter1);
  MorseElement m2(wedge, RightFactor, TBIter2);
  //first test

  CoCoA_ASSERT_ALWAYS(m1 == m2);
  CoCoA_ASSERT_ALWAYS(m1 <= m2);
  CoCoA_ASSERT_ALWAYS(m1 >= m2);

  //second test
  wedge.mySet(0);
  m2.mySetWedgeProduct(wedge);
  //number of true in wedge: 0 < 1
  CoCoA_ASSERT_ALWAYS(m1 < m2);

  //third test
  wedge.mySet(0, false);
  wedge.mySet(1);
  m1.mySetWedgeProduct(wedge);
  // number of trues in wedge: 1 == 1, but x0 * x1^2 * x2 < x0^2 * x1 * x2
  CoCoA_ASSERT_ALWAYS(m1 < m2);

  //fourth test

  TBIter1->elem   = x[0] * x[2];
  TBIter1->lexPos = 1;
  TBIter2->elem   = x[1] * x[2];
  TBIter2->lexPos = 2;
  // number of trues in wedge: 1 == 1
  // x1 * 1 * (x0 * x2) == x0 * 1 * (x1 * x2)
  // m1.ncrit == m2.ncrit
  // m1.wedgeProduct = (0,1,0,0,0), m2.wedgeProduct = (1,0,0,0,0) => m1 < m2
  CoCoA_ASSERT_ALWAYS(m1 < m2);

  //fifth test
  ncrit1.mySet(4);
  ncrit1.mySet(3);
  ncrit1.mySet(2);
  TBIter1->multVars = ncrit1;
  // number of trues in wedge: 1 == 1
  // x1 * 1 * (x0 * x2) == x0 * 1 * (x1 * x2)
  // m1.ncrit == (0,0,1,1,1) > (0,0,0,0,0) == m2.ncrit => m1 < m2
  CoCoA_ASSERT_ALWAYS(m1 < m2);
  CoCoA_ASSERT_ALWAYS(m2 > m1);
}

void test_myEpsilon()
{
  ring QQ = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(QQ, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem>input;
  RingElem LeftFactor(x[0]+x[1]+x[2]+x[3]+x[4]);
  DynamicBitset wedge(5);
  PPMonoidElem RightFactor(one(PPM(PolyRing)));
  RingElem BasisElement(x[1] * x[2] * x[0]);
  DynamicBitset ncrit(5);
  std::vector<MorseElement::JBElem> TestBasis;
  TestBasis.push_back(MorseElement::JBElem(BasisElement, ncrit, 1));

  //first test
  MorseElement m1(wedge, RightFactor, TestBasis.begin());

  DynamicBitset bitset(5);
  m1.mySetWedgeProduct(bitset);
  CoCoA_ASSERT_ALWAYS( 1 == m1.myEpsilon(3, 4));

  bitset.mySet(3, true);
  bitset.mySet(2, true);
  m1.mySetWedgeProduct(bitset);
  CoCoA_ASSERT_ALWAYS( -1 == m1.myEpsilon(3, 4));
  CoCoA_ASSERT_ALWAYS( 1 == m1.myEpsilon(4, 4));
  CoCoA_ASSERT_ALWAYS( 1 == m1.myEpsilon(2, 4));
  CoCoA_ASSERT_ALWAYS( 1 == m1.myEpsilon(0, 4));
  CoCoA_ASSERT_ALWAYS( 1 == m1.myEpsilon(0, 0));
}

void test_myMaxTypeOne()
{
  ring QQ = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(QQ, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem>input;
  RingElem LeftFactor(x[0]);
  DynamicBitset wedge(5);
  PPMonoidElem RightFactor(one(PPM(PolyRing)));
  RingElem BasisElement(x[1] * x[2] * x[0]);
  DynamicBitset ncrit(5);
  std::vector<MorseElement::JBElem> TestBasis;
  TestBasis.push_back(MorseElement::JBElem(BasisElement, ncrit, 1));
  std::vector<MorseElement::JBElem>::iterator TBIter1(TestBasis.begin());
  //first test
  MorseElement m1(wedge, RightFactor, TestBasis.begin());
  CoCoA_ASSERT_ALWAYS(-1 == m1.myMaxTypeOne());

  wedge.mySet(0, true);
  wedge.mySet(1, true);
  wedge.mySet(2, false);
  wedge.mySet(3, false);
  wedge.mySet(4, false);
  m1.mySetWedgeProduct(wedge);
  m1.mySetRightFactor(LPP(x[0] * x[1]));
  ncrit.mySet(0, true);
  ncrit.mySet(1, true);
  ncrit.mySet(2, true);
  ncrit.mySet(3, true);
  ncrit.mySet(4, true);
  TBIter1->multVars = ncrit;
  CoCoA_ASSERT_ALWAYS(-1 == m1.myMaxTypeOne());

  wedge.mySet(1, false);
  m1.mySetWedgeProduct(wedge);
  CoCoA_ASSERT_ALWAYS(1 == m1.myMaxTypeOne());

  m1.mySetRightFactor(LPP(x[0] * x[1] * x[3]));
  CoCoA_ASSERT_ALWAYS(3 == m1.myMaxTypeOne());
}


void test_myMaxTypeTwo()
{
  ring QQ = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(QQ, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  RingElem LeftFactor(x[0]);
  DynamicBitset wedge(5);
  PPMonoidElem RightFactor(one(PPM(PolyRing)));
  RingElem BasisElement(x[1] * x[2] * x[0]);
  DynamicBitset ncrit(5);
  std::vector<MorseElement::JBElem> TestBasis;
  TestBasis.push_back(MorseElement::JBElem(BasisElement, ncrit, 1));
  std::vector<MorseElement::JBElem>::iterator TBIter1(TestBasis.begin());

  //first test
  MorseElement m1(wedge, RightFactor, TestBasis.begin());
  CoCoA_ASSERT_ALWAYS(-1 == m1.myMaxTypeTwo());

  wedge.mySet(0, true);
  wedge.mySet(1, true);
  wedge.mySet(2, false);
  wedge.mySet(3, false);
  wedge.mySet(4, false);
  m1.mySetWedgeProduct(wedge);
  m1.mySetRightFactor(LPP(x[0] * x[1]));
  ncrit.mySet(0, true);
  ncrit.mySet(1, true);
  ncrit.mySet(2, true);
  ncrit.mySet(3, true);
  ncrit.mySet(4, true);
  TBIter1->multVars = ncrit;
  CoCoA_ASSERT_ALWAYS(1 == m1.myMaxTypeTwo());

  wedge.mySet(1, false);
  m1.mySetWedgeProduct(wedge);
  CoCoA_ASSERT_ALWAYS(1 == m1.myMaxTypeTwo());

  m1.mySetRightFactor(LPP(x[0] * x[1] * x[3]));
  CoCoA_ASSERT_ALWAYS(3 == m1.myMaxTypeTwo());
}


void test_myComputeBasicMaps()
{
  ring QQ = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(QQ, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem>input;
  input.push_back(x[0]);
  input.push_back(x[1]);
  input.push_back(x[2]);
  input.push_back(x[3]);
  input.push_back(x[4]);
  JBMill mill(JBMill::Builder().setInput(input));

  //if there is no wedge product there should not be a map
  DynamicBitset wedge(5);
  PPMonoidElem RightFactor(one(PPM(PolyRing)));
  RingElem BasisElement(x[1] * x[2] * x[0]);
  DynamicBitset ncrit(5);
  StandardRepresentationContainer container(mill);
  std::map<PPMonoidElem, std::vector<bool> > NonMultVars(mill.myNonMultVars());
  std::vector<MorseElement::JBElem> TestBasis;
  TestBasis.push_back(MorseElement::JBElem(BasisElement, ncrit, 1));
  MorseElement m1(wedge, RightFactor, TestBasis.begin());
  std::pair<vector<MorseElement::JBElem>::const_iterator, vector<MorseElement::JBElem>::const_iterator> pair(TestBasis.begin(), TestBasis.end());
  CoCoA_ASSERT_ALWAYS(m1.myComputeBasicMaps(pair, container).size() == 0);


  //maps of 0/\1 * (1) * (x4) with respect to pommaret basis (x0, x1, x2, x3, x4)
  wedge.mySet(0, true);
  wedge.mySet(1, true);

  TestBasis.clear();
  TestBasis.push_back(MorseElement::JBElem(x[4], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[4]))), 1));
  TestBasis.push_back(MorseElement::JBElem(x[3], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[3]))), 2));
  TestBasis.push_back(MorseElement::JBElem(x[2], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[2]))), 3));
  TestBasis.push_back(MorseElement::JBElem(x[1], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[1]))), 4));
  TestBasis.push_back(MorseElement::JBElem(x[0], flip(VectorBoolToDynamicBitset(mill.myNonMultVarsOf(x[0]))), 5));
  std::vector<MorseElement::JBElem>::iterator TBIter(TestBasis.begin());
  MorseElement m2(wedge, RightFactor, TBIter);
  pair.first = TestBasis.begin();
  pair.second = TestBasis.end();
  std::vector<std::pair<MorseElement, RingElem> > maps(m2.myComputeBasicMaps(pair, container));
  CoCoA_ASSERT_ALWAYS(maps.size() == 3);


  MorseElement MorseMap1 = maps[0].first;
  RingElem map1 = maps[0].second;

  DynamicBitset ExpectedWedge1(5);
  ExpectedWedge1.mySet(1, true);
  DynamicBitset ExpectedNCrit1(5);
  ExpectedNCrit1.mySet(4, true);

  CoCoA_ASSERT_ALWAYS(ExpectedWedge1 == MorseMap1.myGetWedgeProduct());
  CoCoA_ASSERT_ALWAYS(RightFactor == MorseMap1.myGetRightFactor());
  CoCoA_ASSERT_ALWAYS(x[4] == MorseMap1.myGetBasisElement());
  CoCoA_ASSERT_ALWAYS(ExpectedNCrit1 == MorseMap1.myGetNCrit());
  CoCoA_ASSERT_ALWAYS(map1 == x[0]);


  // MorseElement MorseMap2 = maps[1].first;
  // RingElem map2 = maps[1].second;

  // DynamicBitset ExpectedWedge2(ExpectedWedge1);
  // DynamicBitset ExpectedNCrit2(5);
  // ExpectedNCrit2 = flip(ExpectedNCrit2);

  // CoCoA_ASSERT_ALWAYS(ExpectedWedge2 == MorseMap2.myGetWedgeProduct());
  // CoCoA_ASSERT_ALWAYS(LPP(x[4]) == MorseMap2.myGetRightFactor());
  // CoCoA_ASSERT_ALWAYS(x[0] == MorseMap2.myGetBasisElement());
  // CoCoA_ASSERT_ALWAYS(ExpectedNCrit2 == MorseMap2.myGetNCrit());
  // CoCoA_ASSERT_ALWAYS(IsOne(map2 * (-1)));


  MorseElement MorseMap3 = maps[1].first;
  RingElem map3 = maps[1].second;

  DynamicBitset ExpectedWedge3(5);
  ExpectedWedge3.mySet(0, true);
  DynamicBitset ExpectedNCrit3(ExpectedNCrit1);

  CoCoA_ASSERT_ALWAYS(ExpectedWedge3 == MorseMap3.myGetWedgeProduct());
  CoCoA_ASSERT_ALWAYS(RightFactor == MorseMap3.myGetRightFactor());
  CoCoA_ASSERT_ALWAYS(x[4] == MorseMap3.myGetBasisElement());
  CoCoA_ASSERT_ALWAYS(ExpectedNCrit3 == MorseMap3.myGetNCrit());
  CoCoA_ASSERT_ALWAYS(map3 == -x[1]);


  MorseElement MorseMap4 = maps[2].first;
  RingElem map4 = maps[2].second;

  DynamicBitset ExpectedWedge4(ExpectedWedge3);
  DynamicBitset ExpectedNCrit4(5);
  ExpectedNCrit4 = flip(ExpectedNCrit4);
  ExpectedNCrit4.mySet(0, false);

  CoCoA_ASSERT_ALWAYS(ExpectedWedge4 == MorseMap4.myGetWedgeProduct());
  CoCoA_ASSERT_ALWAYS(LPP(x[4]) == MorseMap4.myGetRightFactor());
  CoCoA_ASSERT_ALWAYS(x[1] == MorseMap4.myGetBasisElement());
  CoCoA_ASSERT_ALWAYS(ExpectedNCrit4 == MorseMap4.myGetNCrit());
  CoCoA_ASSERT_ALWAYS(IsOne(map4));
}

void test_StandardRepresentationContainer()
{
  ring QQ = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(QQ, SymbolRange("x",0,2));
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(power(x[0],2)+ 4*power(x[1],2) - 4);
  input.push_back(- x[0]+2*power(x[1],2));
  JBMill mill(JBMill::Builder().setInput(input));
  StandardRepresentationContainer container(mill);

  RingElem t1 = x[0] * x[0] * x[0] + x[1] * x[0] + 10;
  RingElem t2 = x[0] * x[0] * x[0] + x[2] * x[0] + 10;
  RingElem t3 = x[0] * x[0] * x[0] + x[2] * x[2] + 10;


  std::vector<RingElem> ResContainer(container.myComputeStandardRepresentation(t1));
  std::vector<std::pair<RingElem, RingElem> > ResMill(mill.myStandardRepresentationWithoutRest(t1));
  std::vector<std::pair<RingElem, RingElem> >::iterator i2(ResMill.begin());
  for (std::vector<RingElem>::iterator i = ResContainer.begin(); i != ResContainer.end(); ++i)
  {
    CoCoA_ASSERT_ALWAYS((*i) == i2->second);
    ++i2;
  }
  ResContainer = container.myComputeStandardRepresentation(t2);
  ResMill = mill.myStandardRepresentationWithoutRest(t2);
  i2 = ResMill.begin();
  for (std::vector<RingElem>::iterator i = ResContainer.begin(); i != ResContainer.end(); ++i)
  {
    CoCoA_ASSERT_ALWAYS((*i) == i2->second);
    ++i2;
  }
  ResContainer = container.myComputeStandardRepresentation(t3);
  ResMill = mill.myStandardRepresentationWithoutRest(t3);
  i2 = ResMill.begin();
  for (std::vector<RingElem>::iterator i = ResContainer.begin(); i != ResContainer.end(); ++i)
  {
    CoCoA_ASSERT_ALWAYS((*i) == i2->second);
    ++i2;
  }
}

void program()
{
  GlobalManager CoCoAFoundations;
  test_myDynamicBitsetToLong();
  test_myVectorLongToDynamicBitset();
  test_operators();
  test_myEpsilon();
  test_myMaxTypeOne();
  test_myMaxTypeTwo();
  test_myComputeBasicMaps();
  test_StandardRepresentationContainer();
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
