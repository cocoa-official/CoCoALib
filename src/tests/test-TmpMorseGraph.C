//   Copyright (c)  2013-2015  John Abbott,  Anna M. Bigatti
//   Original author:  2013-2015  Mario Albert

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
#include "CoCoA/TmpMorseBetti.H"
#include "CoCoA/TmpMorseResolution.H"
#include "CoCoA/ideal.H"
#include "CoCoA/library.H"
#include "CoCoA/error.H"

using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

// private
// void test_myVariationWithoutRepetition(bool WriteName)
// {
//   if(WriteName)
//   {
//     std::cout << "test_myVariationWithoutRepetition" << std::endl;
//   }
//   ring Q = RingQQ();
//   SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegRevLex);
//   const vector<RingElem> x = indets(PolyRing);
//   vector<RingElem> PolInput;
//   PolInput.push_back(x[0]);
//   PolInput.push_back(x[1]);
//   PolInput.push_back(x[2]);
//   PolInput.push_back(x[3]);
//   PolInput.push_back(x[4]);
//   jbmill mill(JBMill::Builder().setInput(PolInput));


//   MorseResolution mg(mill);
//   std::vector<std::vector<long> > TestResult;
//   std::vector<long> input;
//   input.push_back(0);
//   input.push_back(1);
//   input.push_back(2);
//   input.push_back(3);
//   input.push_back(4);

//   // test with length 0
//   mg.myVariationWithoutRepetition(TestResult, std::vector<long>(), input, 0);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].size() == 0);

//   // test with length 6
//   TestResult.clear();
//   mg.myVariationWithoutRepetition(TestResult, std::vector<long>(), input, 6);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 0);

//   // test with length 1
//   TestResult.clear();
//   mg.myVariationWithoutRepetition(TestResult, std::vector<long>(), input, 1);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 5);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[0][0] == 0);
//   CoCoA_ASSERT_ALWAYS(TestResult[1][0] == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[2][0] == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[3][0] == 3);
//   CoCoA_ASSERT_ALWAYS(TestResult[4][0] == 4);

//   // test with length 2
//   TestResult.clear();
//   mg.myVariationWithoutRepetition(TestResult, std::vector<long>(), input, 2);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 10);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[5].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[6].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[7].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[8].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[9].size() == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[0][0] == 0);
//   CoCoA_ASSERT_ALWAYS(TestResult[0][1] == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[1][0] == 0);
//   CoCoA_ASSERT_ALWAYS(TestResult[1][1] == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[2][0] == 0);
//   CoCoA_ASSERT_ALWAYS(TestResult[2][1] == 3);
//   CoCoA_ASSERT_ALWAYS(TestResult[3][0] == 0);
//   CoCoA_ASSERT_ALWAYS(TestResult[3][1] == 4);
//   CoCoA_ASSERT_ALWAYS(TestResult[4][0] == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[4][1] == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[5][0] == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[5][1] == 3);
//   CoCoA_ASSERT_ALWAYS(TestResult[6][0] == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[6][1] == 4);
//   CoCoA_ASSERT_ALWAYS(TestResult[7][0] == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[7][1] == 3);
//   CoCoA_ASSERT_ALWAYS(TestResult[8][0] == 2);
//   CoCoA_ASSERT_ALWAYS(TestResult[8][1] == 4);
//   CoCoA_ASSERT_ALWAYS(TestResult[9][0] == 3);
//   CoCoA_ASSERT_ALWAYS(TestResult[9][1] == 4);
// }

// private
// void test_myPossibleWedges(bool WriteName)
// {
//   if(WriteName)
//   {
//     std::cout << "test_myPossibleWedges" << std::endl;
//   }
//   ring Q = RingQQ();
//   SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegRevLex);
//   const vector<RingElem> x = indets(PolyRing);
//   vector<RingElem> PolInput;
//   PolInput.push_back(x[0]);
//   PolInput.push_back(x[1]);
//   PolInput.push_back(x[2]);
//   PolInput.push_back(x[3]);
//   PolInput.push_back(x[4]);
//   JBMill mill(JBMill::Builder().setInput(PolInput));


//   MorseResolution mg(mill);
//   std::vector<DynamicBitset> TestResult;
//   std::vector<long> input;

//   TestResult = mg.myPossibleWedges(input, 0);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 1);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].IamAll0s());


//   input.push_back(0);
//   input.push_back(1);
//   input.push_back(2);
//   TestResult = mg.myPossibleWedges(input, 3);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 8);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(2) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(2) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(2) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(2) == true);

//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(2) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(2) == true);

//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(2) == true);

//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(2) == true);


//   TestResult = mg.myPossibleWedges(input, 4);
//   CoCoA_ASSERT_ALWAYS(TestResult.size() == 8);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(2) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[0].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(2) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[1].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(2) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[2].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(2) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[3].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(2) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[4].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(1) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(2) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[5].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(0) == false);
//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(2) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[6].Iam1At(3) == false);

//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(0) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(1) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(2) == true);
//   CoCoA_ASSERT_ALWAYS(TestResult[7].Iam1At(3) == false);
// }

// void test_myComputeGeneralBasis(bool WriteName)
// {
//   if(WriteName)
//   {
//     std::cout << "test_myComputeGeneralBasis" << std::endl;
//   }
//   ring Q = RingQQ();
//   SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,1), StdDegRevLex);
//   const vector<RingElem> x = indets(PolyRing);
//   vector<RingElem> PolInput;
//   PolInput.push_back(x[0]);
//   PolInput.push_back(x[1]);
//   JBMill mill(JBMill::Builder().setInput(PolInput));

//   MorseResolution mg(mill);

//   std::vector<MorseElement> elems(mg.myComputeGeneralBasis());

//   CoCoA_ASSERT_ALWAYS(elems.size() == 3);
//   std::vector<MorseElement::JBElem > TestBasis;
//   TestBasis.push_back(MorseElement::JBElem(x[0], DynamicBitset(LPP(x[0] * x[1]))));
//   TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1]))));
//   TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1]))));
//   std::vector<MorseElement::JBElem >::iterator TBIter(TestBasis.begin());
//   MorseElement m1(DynamicBitset(2), TBIter);
//   ++TBIter;
//   MorseElement m2(DynamicBitset(2), TBIter);
//   ++TBIter;
//   MorseElement m3(DynamicBitset(LPP(x[0])), TBIter);
//   CoCoA_ASSERT_ALWAYS(elems[0] == m2);
//   CoCoA_ASSERT_ALWAYS(elems[1] == m3);
//   CoCoA_ASSERT_ALWAYS(elems[2] == m1);
// }

//  void test_myComputeBasicGraph(bool WriteName)
// {
//   if(WriteName)
//   {
//     std::cout << " test_myComputeBasicGraph" << std::endl;
//   }
//   ring Q = RingQQ();
//   SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,1), StdDegRevLex);
//   const vector<RingElem> x = indets(PolyRing);
//   vector<RingElem> PolInput;
//   PolInput.push_back(x[0]);
//   PolInput.push_back(x[1]);
//   JBMill mill(JBMill::Builder().setInput(PolInput));

//   MorseResolution mg(mill);
//   StandardRepresentationContainer container(mill);
//   mg.myComputeBasicGraph(mg.myComputeGeneralBasis(), container);

//   std::map<MorseElement, MorsePaths> res(mg.myGetBasicResolution());

//   CoCoA_ASSERT_ALWAYS(res.size() == 4);
//   std::vector<MorseElement::JBElem > TestBasis;
//   TestBasis.push_back(MorseElement::JBElem(x[0], DynamicBitset(LPP(x[0] * x[1]))));
//   TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1]))));
//   TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1]))));
//   std::vector<MorseElement::JBElem >::iterator TBIter(TestBasis.begin());

//   MorseElement m1(DynamicBitset(2), TBIter);
//   ++TBIter;
//   MorseElement m2(DynamicBitset(2), TBIter);
//   ++TBIter;
//   MorseElement m3(DynamicBitset(LPP(x[0])), TBIter);
//   TBIter = TestBasis.begin();
//   MorseElement m4(DynamicBitset(2), LPP(x[1]), TBIter);

//   CoCoA_ASSERT_ALWAYS(res.find(m1) != res.end());
//   CoCoA_ASSERT_ALWAYS(res.find(m2) != res.end());
//   CoCoA_ASSERT_ALWAYS(res.find(m3) != res.end());
//   CoCoA_ASSERT_ALWAYS(res.find(m4) != res.end());

//   CoCoA_ASSERT_ALWAYS(res[m1].IamEmpty() == true);
//   CoCoA_ASSERT_ALWAYS(res[m2].IamEmpty() == false);
//   CoCoA_ASSERT_ALWAYS(res[m3].IamEmpty() == true);
//   CoCoA_ASSERT_ALWAYS(res[m4].IamEmpty() == false);

//   CoCoA_ASSERT_ALWAYS(res[m2].myGetPaths().size() == 1);
//   CoCoA_ASSERT_ALWAYS(res[m4].myGetPaths().size() == 1);
// }

void outputResolution(std::map<MorseElement, MorsePaths> res)
{
  for (std::map<MorseElement, MorsePaths>::iterator i = res.begin(); i != res.end(); ++i)
  {
    std::cout << "=========================================" << std::endl;
    std::cout <<  "MorseElement = " << i->first << std::endl;
    PathMap maps = (i->second).myGetPaths();
    for (PathMap::iterator j = maps.begin(); j != maps.end(); ++j)
    {
      std::cout << ((j->first)->first) << " ----- " << (j->second) << " ----> " << i->first << std::endl;
    }
  }
}
void test_myDirectMorseReduction1(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_myDirectMorseReduction1" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,1), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[0]);
  PolInput.push_back(x[1]);
  JBMill mill(JBMill::Builder().setInput(PolInput));

  MorseResolution mg(mill);
  mg.myComputeResolution();

  std::map<MorseElement, MorsePaths> res(mg.myGetBasicResolution());
  CoCoA_ASSERT_ALWAYS(res.size() == 3);
  std::vector<MorseElement::JBElem > TestBasis;
  TestBasis.push_back(MorseElement::JBElem(x[0], DynamicBitset(LPP(x[0] * x[1])), 1));
  TestBasis.push_back(MorseElement::JBElem(x[1], DynamicBitset(LPP(x[1]))       , 0));
  std::vector<MorseElement::JBElem >::iterator TBIter1(TestBasis.begin());
  std::vector<MorseElement::JBElem >::iterator TBIter2(TestBasis.begin());
  std::advance(TBIter2, 1);

  MorseElement m1(DynamicBitset(2), TBIter1);
  MorseElement m2(DynamicBitset(2), TBIter2);
  MorseElement m3(DynamicBitset(LPP(x[0])), TBIter2);

  CoCoA_ASSERT_ALWAYS(res.find(m1) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m2) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m3) != res.end());

  CoCoA_ASSERT_ALWAYS(res[m1].IamEmpty() == false);
  CoCoA_ASSERT_ALWAYS(res[m2].IamEmpty() == false);
  CoCoA_ASSERT_ALWAYS(res[m3].IamEmpty() == true);

  CoCoA_ASSERT_ALWAYS(res[m1].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m2].myGetPaths().size() == 1);

  CoCoA_ASSERT_ALWAYS(res[m1].myGetPath(res.find(m3)) == -x[1]);
  CoCoA_ASSERT_ALWAYS(res[m2].myGetPath(res.find(m3)) == x[0]);
}

void test_myDirectMorseReduction2(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_myDirectMorseReduction2" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,3), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[2] + 1);
  PolInput.push_back(x[0] * x[2] * x[3]);
  PolInput.push_back(x[0] + x[1] * x[1]);
  PolInput.push_back(x[1] * x[2] + x[3] * x[1]);

  JBMill mill(JBMill::Builder().setInput(PolInput));
  std::vector<RingElem> basis(mill.myReturnJB());

  MorseResolution mg(mill);
  mg.myComputeResolution();

  std::map<MorseElement, MorsePaths> res(mg.myGetBasicResolution());

  std::vector<MorseElement::JBElem > TestBasis;
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[2] + x[1], DynamicBitset(LPP(x[3] * x[2])),               2));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[2] + 1          , DynamicBitset(LPP(x[2] * x[3])),               0));
  TestBasis.push_back(MorseElement::JBElem(x[2] + 1          , DynamicBitset(LPP(x[2] * x[3])),               0));
  TestBasis.push_back(MorseElement::JBElem(x[2] + 1          , DynamicBitset(LPP(x[2] * x[3])),               0));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[2] + x[1], DynamicBitset(LPP(x[2] * x[3])),               2));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[2] + x[1], DynamicBitset(LPP(x[2] * x[3])),               2));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[1]       , DynamicBitset(LPP(x[1] * x[2] * x[3])),        3));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[2] + x[1], DynamicBitset(LPP(x[2] * x[3])),               2));
  TestBasis.push_back(MorseElement::JBElem(x[2] + 1          , DynamicBitset(LPP(x[2] * x[3])),               0));
  TestBasis.push_back(MorseElement::JBElem(x[0]              , DynamicBitset(LPP(x[0] * x[1] * x[2] * x[3])), 4));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[3] - x[1], DynamicBitset(LPP(x[3])),                      1));
  TestBasis.push_back(MorseElement::JBElem(x[1] * x[1]       , DynamicBitset(LPP(x[1] * x[2] * x[3])),        3));

  std::vector<MorseElement::JBElem >::iterator TBIter(TestBasis.begin());
  MorseElement m1(DynamicBitset(LPP(x[0] * x[1] * x[2])), TBIter);
  ++TBIter;

  MorseElement m2(DynamicBitset(LPP(x[0] * x[1])), TBIter);
  ++TBIter;
  MorseElement m3(DynamicBitset(LPP(x[0] * x[1])), TBIter);
  ++TBIter;
  MorseElement m4(DynamicBitset(LPP(x[0] * x[2])), TBIter);
  ++TBIter;
  MorseElement m5(DynamicBitset(LPP(x[1] * x[2])), TBIter);
  ++TBIter;
  MorseElement m6(DynamicBitset(LPP(x[0] * x[1])), TBIter);
  ++TBIter;

  MorseElement m7(DynamicBitset(LPP(x[1])),  TBIter);
  ++TBIter;
  MorseElement m8(DynamicBitset(LPP(x[0])),  TBIter);
  ++TBIter;
  MorseElement m9(DynamicBitset(LPP(x[0])),  TBIter);
  ++TBIter;
  MorseElement m10(DynamicBitset(LPP(x[2])), TBIter);
  ++TBIter;
  MorseElement m11(DynamicBitset(LPP(x[1])), TBIter);
  ++TBIter;
  MorseElement m12(DynamicBitset(LPP(x[1])), TBIter);
  ++TBIter;
  MorseElement m13(DynamicBitset(LPP(x[0])), TBIter);
  ++TBIter;
  MorseElement m14(DynamicBitset(LPP(x[0])), TBIter);
  ++TBIter;

  MorseElement m15(DynamicBitset(LPP(one(PolyRing))), TBIter);
  ++TBIter;
  MorseElement m16(DynamicBitset(LPP(one(PolyRing))), TBIter);
  ++TBIter;
  MorseElement m17(DynamicBitset(LPP(one(PolyRing))), TBIter);
  ++TBIter;
  MorseElement m18(DynamicBitset(LPP(one(PolyRing))), TBIter);
  ++TBIter;
  MorseElement m19(DynamicBitset(LPP(one(PolyRing))), TBIter);

  CoCoA_ASSERT_ALWAYS(res.find(m1) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m2) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m3) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m4) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m5) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m6) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m7) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m8) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m9) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m10) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m11) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m12) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m13) != res.end());
  // outputResolution(res);
  CoCoA_ASSERT_ALWAYS(res.find(m14) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m15) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m16) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m17) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m18) != res.end());
  CoCoA_ASSERT_ALWAYS(res.find(m19) != res.end());


  CoCoA_ASSERT_ALWAYS(res[m1].myGetPaths().size() == 0);
  CoCoA_ASSERT_ALWAYS(res[m2].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m3].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m4].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m5].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m6].myGetPaths().size() == 0);
  CoCoA_ASSERT_ALWAYS(res[m7].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m8].myGetPaths().size() == 1);
  CoCoA_ASSERT_ALWAYS(res[m9].myGetPaths().size() == 3);
  CoCoA_ASSERT_ALWAYS(res[m10].myGetPaths().size() == 2);
  CoCoA_ASSERT_ALWAYS(res[m11].myGetPaths().size() == 2);
  CoCoA_ASSERT_ALWAYS(res[m12].myGetPaths().size() == 2);
  CoCoA_ASSERT_ALWAYS(res[m13].myGetPaths().size() == 2);
  CoCoA_ASSERT_ALWAYS(res[m14].myGetPaths().size() == 2);
  CoCoA_ASSERT_ALWAYS(res[m15].myGetPaths().size() == 4);
  CoCoA_ASSERT_ALWAYS(res[m16].myGetPaths().size() == 2);
  CoCoA_ASSERT_ALWAYS(res[m17].myGetPaths().size() == 4);
  CoCoA_ASSERT_ALWAYS(res[m18].myGetPaths().size() == 3);
  CoCoA_ASSERT_ALWAYS(res[m19].myGetPaths().size() == 3);

  //paths from m1
  CoCoA_ASSERT_ALWAYS(res[m2].myGetPath(res.find(m1)) == -x[3] + 1);
  CoCoA_ASSERT_ALWAYS(res[m3].myGetPath(res.find(m1)) == x[2] + 1);
  CoCoA_ASSERT_ALWAYS(res[m4].myGetPath(res.find(m1)) == -x[1]);
  CoCoA_ASSERT_ALWAYS(res[m5].myGetPath(res.find(m1)) == x[0]);

  //paths from m6
  CoCoA_ASSERT_ALWAYS(res[m7].myGetPath(res.find(m6)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m8].myGetPath(res.find(m6)) == -x[1]);
  CoCoA_ASSERT_ALWAYS(res[m9].myGetPath(res.find(m6)) == one(PolyRing));

  //paths from m5
  CoCoA_ASSERT_ALWAYS(res[m10].myGetPath(res.find(m5)) == x[1]);
  CoCoA_ASSERT_ALWAYS(res[m11].myGetPath(res.find(m5)) == x[3] - 1);
  CoCoA_ASSERT_ALWAYS(res[m12].myGetPath(res.find(m5)) == -x[2] - 1);

  //paths from m4
  CoCoA_ASSERT_ALWAYS(res[m10].myGetPath(res.find(m4)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m13].myGetPath(res.find(m4)) == - x[2] - 1);
  CoCoA_ASSERT_ALWAYS(res[m9].myGetPath(res.find(m4)) == x[3] - 1);

  //paths from m3
  CoCoA_ASSERT_ALWAYS(res[m12].myGetPath(res.find(m3)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m13].myGetPath(res.find(m3)) == - x[1]);
  CoCoA_ASSERT_ALWAYS(res[m14].myGetPath(res.find(m3)) == x[3] - 1);

  //paths from m2
  CoCoA_ASSERT_ALWAYS(res[m11].myGetPath(res.find(m2)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m9].myGetPath(res.find(m2)) == - x[1]);
  CoCoA_ASSERT_ALWAYS(res[m14].myGetPath(res.find(m2)) == x[2] + 1);

  //paths from m7
  CoCoA_ASSERT_ALWAYS(res[m15].myGetPath(res.find(m7)) == - 1);
  CoCoA_ASSERT_ALWAYS(res[m16].myGetPath(res.find(m7)) == x[1]);

  //paths from m8
  CoCoA_ASSERT_ALWAYS(res[m16].myGetPath(res.find(m8)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m17].myGetPath(res.find(m8)) == - x[2] - 1);

  //paths from m9
  CoCoA_ASSERT_ALWAYS(res[m15].myGetPath(res.find(m9)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m17].myGetPath(res.find(m9)) == - x[1] * x[2] - x[1]);

  //paths from m10
  CoCoA_ASSERT_ALWAYS(res[m18].myGetPath(res.find(m10)) == x[2] + 1);
  CoCoA_ASSERT_ALWAYS(res[m15].myGetPath(res.find(m10)) == - x[3] + 1);

  //paths from m11
  CoCoA_ASSERT_ALWAYS(res[m15].myGetPath(res.find(m11)) == x[1]);
  CoCoA_ASSERT_ALWAYS(res[m19].myGetPath(res.find(m11)) == - 1 - x[2]);

  //paths from m12
  CoCoA_ASSERT_ALWAYS(res[m18].myGetPath(res.find(m12)) == x[1]);
  CoCoA_ASSERT_ALWAYS(res[m19].myGetPath(res.find(m12)) == 1 - x[3]);

  //paths from m13
  CoCoA_ASSERT_ALWAYS(res[m18].myGetPath(res.find(m13)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m17].myGetPath(res.find(m13)) == x[1] - x[1] * x[3]);

  //paths from m14
  CoCoA_ASSERT_ALWAYS(res[m19].myGetPath(res.find(m14)) == x[0]);
  CoCoA_ASSERT_ALWAYS(res[m17].myGetPath(res.find(m14)) == - x[1] * x[1]);
}


void test_myMapsAsMatrices(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_myMapsAsMatrices" << std::endl;
  }
  // using result of syzygys
  // output of matrices is transposed content of syzygys
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,3), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[2] + 1);
  PolInput.push_back(x[0] * x[2] * x[3]);
  PolInput.push_back(x[0] + x[1] * x[1]);
  PolInput.push_back(x[1] * x[2] + x[3] * x[1]);

  JBMill mill(JBMill::Builder().setInput(PolInput));
  std::vector<RingElem> basis(mill.myReturnJB());
  MorseResolution mg(mill);
  mg.myComputeResolution();

  std::vector<matrix> syzs(mg.myMapsAsMatrices());
  CoCoA_ASSERT_ALWAYS(syzs.size() == 3);
  matrix syz0(syzs[0]);
  matrix syz1(syzs[1]);
  matrix syz2(syzs[2]);
  CoCoA_ASSERT_ALWAYS(NumRows(syz0) == 5 && NumCols(syz0) == 8);
  CoCoA_ASSERT_ALWAYS(NumRows(syz1) == 8 && NumCols(syz1) == 5);
  CoCoA_ASSERT_ALWAYS(NumRows(syz2) == 5 && NumCols(syz2) == 1);

  matrix tsyz0(NewDenseMat(transpose(syz0)));
  matrix tsyz1(NewDenseMat(transpose(syz1)));
  matrix tsyz2(NewDenseMat(transpose(syz2)));

  CoCoA_ASSERT_ALWAYS(tsyz0(0, 0) == x[1] && tsyz0(0, 1) == 0                   && tsyz0(0, 2) == 0        && tsyz0(0, 3) == -1         && tsyz0(0, 4) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz0(1, 0) == x[0] && tsyz0(1, 1) == -x[2] - 1           && tsyz0(1, 2) == 0        && tsyz0(1, 3) == 0          && tsyz0(1, 4) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz0(2, 0) == 0    && tsyz0(2, 1) == 0                   && tsyz0(2, 2) == x[2] + 1 && tsyz0(2, 3) == -x[3] + 1  && tsyz0(2, 4) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz0(3, 0) == 0    && tsyz0(3, 1) == 0                   && tsyz0(3, 2) == x[1]     && tsyz0(3, 3) == 0          && tsyz0(3, 4) == -x[3] + 1);
  CoCoA_ASSERT_ALWAYS(tsyz0(4, 0) == 0    && tsyz0(4, 1) == -x[1] * x[3] + x[1] && tsyz0(4, 2) == x[0]     && tsyz0(4, 3) == 0          && tsyz0(4, 4) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz0(5, 0) == 0    && tsyz0(5, 1) == 0                   && tsyz0(5, 2) == 0        && tsyz0(5, 3) == x[1]       && tsyz0(5, 4) == -x[2] - 1);
  CoCoA_ASSERT_ALWAYS(tsyz0(6, 0) == 0    && tsyz0(6, 1) == -x[1] * x[2] - x[1] && tsyz0(6, 2) == 0        && tsyz0(6, 3) == x[0]       && tsyz0(6, 4) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz0(7, 0) == 0    && tsyz0(7, 1) == -x[1] * x[1]        && tsyz0(7, 2) == 0        && tsyz0(7, 3) == 0          && tsyz0(7, 4) == x[0]     );

  CoCoA_ASSERT_ALWAYS(tsyz1(0, 0) == x[0]  && tsyz1(0, 1) == -x[1]  && tsyz1(0, 2) == 0    && tsyz1(0, 3) == 0         && tsyz1(0, 4) == 0         && tsyz1(0, 5) == 0         && tsyz1(0, 6) == 1        && tsyz1(0, 7) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz1(1, 0) == 0     && tsyz1(1, 1) == 0      && tsyz1(1, 2) == x[1] && tsyz1(1, 3) == -x[2] - 1 && tsyz1(1, 4) == 0         && tsyz1(1, 5) == x[3] - 1  && tsyz1(1, 6) == 0        && tsyz1(1, 7) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz1(2, 0) == 0     && tsyz1(2, 1) == 0      && tsyz1(2, 2) == x[0] && tsyz1(2, 3) == 0         && tsyz1(2, 4) == -x[2] - 1 && tsyz1(2, 5) == 0         && tsyz1(2, 6) == x[3] - 1 && tsyz1(2, 7) == 0        );
  CoCoA_ASSERT_ALWAYS(tsyz1(3, 0) == 0     && tsyz1(3, 1) == 0      && tsyz1(3, 2) == 0    && tsyz1(3, 3) == x[0]      && tsyz1(3, 4) == -x[1]     && tsyz1(3, 5) == 0         && tsyz1(3, 6) == 0        && tsyz1(3, 7) == x[3] - 1 );
  CoCoA_ASSERT_ALWAYS(tsyz1(4, 0) == 0     && tsyz1(4, 1) == 0      && tsyz1(4, 2) == 0    && tsyz1(4, 3) == 0         && tsyz1(4, 4) == 0         && tsyz1(4, 5) == x[0]      && tsyz1(4, 6) == -x[1]     && tsyz1(4, 7) == x[2] + 1);

  CoCoA_ASSERT_ALWAYS(tsyz2(0, 0) == 0 && tsyz2(0, 1) == x[0] && tsyz2(0, 2) == -x[1] && tsyz2(0, 3) == x[2] + 1 && tsyz2(0, 4) == -x[3] + 1);

}

void test_myBettis(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_myBettis" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[2] + 1);
  PolInput.push_back(x[0] * x[2] * x[3]);
  PolInput.push_back(x[0] + x[1] * x[1]);
  PolInput.push_back(x[1] * x[2] + x[3] * x[1]);
  ideal I(PolyRing, PolInput);
  I = homog(I, x[4]);
  JBMill mill(JBMill::Builder().setInput(gens(I)));
  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);
  matrix bettis(NewDenseMat(RingZZ(), 2, 5));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 0, 1, 2);
  SetEntry(bettis, 0, 2, 1);
  SetEntry(bettis, 0, 3, 0);
  SetEntry(bettis, 0, 4, 0);

  SetEntry(bettis, 1, 0, 0);
  SetEntry(bettis, 1, 1, 2);
  SetEntry(bettis, 1, 2, 5);
  SetEntry(bettis, 1, 3, 4);
  SetEntry(bettis, 1, 4, 1);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}

void test_onlyX(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_onlyX" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,0), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(PolInput));
  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 1, 2));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 0, 1, 1);

  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));

}

void test_onlyOne(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_onlyOne" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing PolyRing = NewPolyRing(Q, SymbolRange("x",0,0), StdDegRevLex);
  const vector<RingElem> x = indets(PolyRing);
  vector<RingElem> PolInput;
  PolInput.push_back(one(PolyRing));
  JBMill mill(JBMill::Builder().setInput(PolInput));
  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 1, 1));

  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));

}


void test_chandra4(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_chandra4" << std::endl;
  }
  ring Q = RingQQ();

  SparsePolyRing polyRingBenchmark_D1 = NewPolyRing_DMPI(Q,5);
  vector<RingElem> x = indets(polyRingBenchmark_D1);
  vector<RingElem> input;
  input.push_back(  - 51234*power(x[0],2) - 34156*x[0]*x[1] - 25617*x[0]*x[2] + 1497532*x[0] - 1600000);
  input.push_back(- 170780*x[0]*x[1] - 128085*power(x[1],2) - 102468*x[1]*x[2] + 7743830*x[1] - 8000000);
  input.push_back(  - 128085*x[0]*x[2] - 102468*x[1]*x[2] - 85390*power(x[2],2) + 7829220*x[2] - 8000000);
  input.push_back(  - 717276*x[0]*x[3] - 597730*x[1]*x[3] - 512340*x[2]*x[3] + 55103405*x[3] - 56000000);
  ideal I(polyRingBenchmark_D1, input);
  I = homog(I, x[4]);
//ring r = 0,(x0, x1, x2, x3, x4),dp;
//ideal i = - 51234*x0^2 - 34156*x0*x1 - 25617*x0*x2 + 1497532*x0 - 1600000, - 170780*x0*x1 - 128085*x1^2 - 102468*x1*x2 + 7743830*x1 - 8000000, - 128085*x0*x2 - 102468*x1*x2 - 85390*x2^2 + 7829220*x2 - 8000000, - 717276*x0*x3 - 597730*x1*x3 - 512340*x2*x3 + 55103405*x3 - 56000000;
//i = groebner(i);
//i = homog(i, x4);
//i = groebner(i);
//resolution rs = mres(i, 0);
//print(betti(rs), "betti");

  JBMill mill(JBMill::Builder().setInput(gens(I)));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 4, 5));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 0, 1, 1);
  SetEntry(bettis, 1, 1, 3);
  SetEntry(bettis, 1, 2, 3);
  SetEntry(bettis, 2, 2, 3);
  SetEntry(bettis, 2, 3, 3);
  SetEntry(bettis, 3, 3, 1);
  SetEntry(bettis, 3, 4, 1);

  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));

}


void test_boon(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_boon" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,7);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;

  input.push_back(power(x[2],2) + power(x[4],2) - 1);
  input.push_back( power(x[3],2) + power(x[5],2) - 1);
  input.push_back(5*x[0]*power(x[2],3) + 5*x[1]*power(x[3],3) - 6);
  input.push_back( 5*x[0]*power(x[4],3) + 5*x[1]*power(x[5],3) - 6);
  input.push_back(      10*x[0]*power(x[2],2)*x[4] + 10*x[1]*power(x[3],2)*x[5] - 7);
  input.push_back(      10*x[0]*x[2]*power(x[4],2) + 10*x[1]*x[3]*power(x[5],2) - 7);

  // ring r = 0,(x0, x1, x2, x3, x4, x5, x6),dp;
  // ideal i = x2^2 + x4^2 - 1, x3^2 + x5^2 - 1, 5*x0*x2^3 + 5*x1*x3^3 - 6, 5*x0*x4^3 + 5*x1*x5^3 - 6, 10*x0*x2^2*x4 + 10*x1*x3^2*x5 - 7, 10*x0*x2*x4^2 + 10*x1*x3*x5^2 - 7;
  // i = groebner(i);
  // i = homog(i, x6);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");

  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[6]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 3, 7));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 0, 1, 2);
  SetEntry(bettis, 0, 2, 1);
  SetEntry(bettis, 1, 1, 7);
  SetEntry(bettis, 1, 2, 22);
  SetEntry(bettis, 1, 3, 26);
  SetEntry(bettis, 1, 4, 14);
  SetEntry(bettis, 1, 5, 3);
  SetEntry(bettis, 2, 2, 6);
  SetEntry(bettis, 2, 3, 20);
  SetEntry(bettis, 2, 4, 25);
  SetEntry(bettis, 2, 5, 14);
  SetEntry(bettis, 2, 6, 3);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}


void test_caprasse(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_caprasse" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,5);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;

  input.push_back(    2*x[3]*x[0]*x[1] - 2*x[0] + power(x[1],2)*x[2] - x[2]);
  input.push_back(  4*x[3]*power(x[0],2)*x[1] + 2*x[3]*power(x[1],3) - 10*x[3]*x[1] - power(x[0],3)*x[2] + 4*power(x[0],2) + 4*x[0]*power(x[1],2)*x[2] + 4*x[0]*x[2] - 10*power(x[1],2) + 2);
  input.push_back(      power(x[3],2)*x[0] + 2*x[3]*x[1]*x[2] - x[0] - 2*x[2]);
  input.push_back(      2*power(x[3],3)*x[1] + 4*power(x[3],2)*x[0]*x[2] - 10*power(x[3],2) + 4*x[3]*x[1]*power(x[2],2) - 10*x[3]*x[1] - x[0]*power(x[2],3) + 4*x[0]*x[2] + 4*power(x[2],2) + 2);

  // ring r = 0,(x0, x1, x2, x3, x4),dp;
  // ideal i = 2*x3*x0*x1 - 2*x0 + x1^2*x2 - x2, 4*x3*x0^2*x1 + 2*x3*x1^3 - 10*x3*x1 - x0^3*x2 + 4*x0^2 + 4*x0*x1^2*x2 + 4*x0*x2 - 10*x1^2 + 2, x3^2*x0 + 2*x3*x1*x2 - x0 - 2*x2, 2*x3^3*x1 + 4*x3^2*x0*x2 - 10*x3^2 + 4*x3*x1*x2^2 - 10*x3*x1 - x0*x2^3 + 4*x0*x2 + 4*x2^2 + 2;
  // i = groebner(i);
  // i = homog(i, x4);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");

  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[4]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 7, 5));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 2, 1, 2);
  SetEntry(bettis, 3, 1, 12);
  SetEntry(bettis, 3, 2, 18);
  SetEntry(bettis, 3, 3, 7);
  SetEntry(bettis, 4, 2, 13);
  SetEntry(bettis, 4, 3, 16);
  SetEntry(bettis, 4, 4, 6);
  SetEntry(bettis, 5, 2, 2);
  SetEntry(bettis, 5, 3, 3);
  SetEntry(bettis, 6, 3, 2);
  SetEntry(bettis, 6, 4, 2);

  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}


void test_cassou(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_cassou" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,5);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(    6*power(x[0],4)*power(x[1],3) + 21*power(x[0],4)*power(x[1],2)*x[2] + 15*power(x[0],4)*x[1]*power(x[2],2) + 9*power(x[0],4)*power(x[2],3) - 8*power(x[0],2)*power(x[1],2)*x[3] - 28*power(x[0],2)*x[1]*x[2]*x[3] - 144*power(x[0],2)*x[1] + 36*power(x[0],2)*power(x[2],2)*x[3] - 648*power(x[0],2)*x[2] - 120);
  input.push_back(  9*power(x[0],4)*power(x[1],4) + 30*power(x[0],4)*power(x[1],3)*x[2] + 39*power(x[0],4)*power(x[1],2)*power(x[2],2) + 18*power(x[0],4)*x[1]*power(x[2],3) - 24*power(x[0],2)*power(x[1],3)*x[3] - 16*power(x[0],2)*power(x[1],2)*x[2]*x[3] - 432*power(x[0],2)*power(x[1],2) + 16*power(x[0],2)*x[1]*power(x[2],2)*x[3] - 720*power(x[0],2)*x[1]*x[2] + 24*power(x[0],2)*power(x[2],3)*x[3] - 432*power(x[0],2)*power(x[2],2) + 16*power(x[1],2)*power(x[3],2) - 32*x[1]*x[2]*power(x[3],2) + 576*x[1]*x[3] - 240*x[1] + 16*power(x[2],2)*power(x[3],2) - 576*x[2]*x[3] + 5184);
  input.push_back(      - 15*power(x[0],2)*power(x[1],3)*x[3] + 15*power(x[0],2)*power(x[1],2)*x[2]*x[3] - 81*power(x[0],2)*power(x[1],2) + 216*power(x[0],2)*x[1]*x[2] - 162*power(x[0],2)*power(x[2],2) + 40*power(x[1],2)*power(x[3],2) - 80*x[1]*x[2]*power(x[3],2) + 1008*x[1]*x[3] + 40*power(x[2],2)*power(x[3],2) - 1008*x[2]*x[3] + 5184);
  input.push_back(- 4*power(x[0],2)*power(x[1],2) + 4*power(x[0],2)*x[1]*x[2] - 3*power(x[0],2)*power(x[2],2) + 22*x[1]*x[3] - 22*x[2]*x[3] + 261);

  // ring r = 0,(x0, x1, x2, x3, x4),dp;
  // ideal i = 6*x0^4*x1^3 + 21*x0^4*x1^2*x2 + 15*x0^4*x1*x2^2 + 9*x0^4*x2^3 - 8*x0^2*x1^2*x3 - 28*x0^2*x1*x2*x3 - 144*x0^2*x1 + 36*x0^2*x2^2*x3 - 648*x0^2*x2 - 120, 9*x0^4*x1^4 + 30*x0^4*x1^3*x2 + 39*x0^4*x1^2*x2^2 + 18*x0^4*x1*x2^3 - 24*x0^2*x1^3*x3 - 16*x0^2*x1^2*x2*x3 - 432*x0^2*x1^2 + 16*x0^2*x1*x2^2*x3 - 720*x0^2*x1*x2 + 24*x0^2*x2^3*x3 - 432*x0^2*x2^2 + 16*x1^2*x3^2 - 32*x1*x2*x3^2 + 576*x1*x3 - 240*x1 + 16*x2^2*x3^2 - 576*x2*x3 + 5184, - 15*x0^2*x1^3*x3 + 15*x0^2*x1^2*x2*x3 - 81*x0^2*x1^2 + 216*x0^2*x1*x2 - 162*x0^2*x2^2 + 40*x1^2*x3^2 - 80*x1*x2*x3^2 + 1008*x1*x3 + 40*x2^2*x3^2 - 1008*x2*x3 + 5184, - 4*x0^2*x1^2 + 4*x0^2*x1*x2 - 3*x0^2*x2^2 + 22*x1*x3 - 22*x2*x3 + 261;
  // i = groebner(i);
  // i = homog(i, x4);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");

  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[4]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 4, 5));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 1, 1, 3);
  SetEntry(bettis, 2, 1, 4);
  SetEntry(bettis, 2, 2, 11);
  SetEntry(bettis, 2, 3, 4);
  SetEntry(bettis, 3, 2, 4);
  SetEntry(bettis, 3, 3, 9);
  SetEntry(bettis, 3, 4, 4);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));

}

void test_conform1(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_conform1" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,4);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;

  input.push_back(  - 3*power(x[1],2)*power(x[2],2) - power(x[1],2) + 8*x[1]*x[2] - power(x[2],2) - 9);
  input.push_back(  - 3*power(x[0],2)*power(x[2],2) - power(x[0],2) + 8*x[0]*x[2] - power(x[2],2) - 9);
  input.push_back(  - 3*power(x[0],2)*power(x[1],2) - power(x[0],2) + 8*x[0]*x[1] - power(x[1],2) - 9);

  // ring r = 0,(x0, x1, x2, x3),dp;
  // ideal i = - 3*x1^2*x2^2 - x1^2 + 8*x1*x2 - x2^2 - 9, - 3*x0^2*x2^2 - x0^2 + 8*x0*x2 - x2^2 - 9, - 3*x0^2*x1^2 - x0^2 + 8*x0*x1 - x1^2 - 9;
  // i = groebner(i);
  // i = homog(i, x3);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");

  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[3]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 5, 4));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 2, 1, 5);
  SetEntry(bettis, 2, 2, 2);
  SetEntry(bettis, 3, 1, 1);
  SetEntry(bettis, 3, 2, 6);
  SetEntry(bettis, 3, 3, 2);
  SetEntry(bettis, 4, 3, 1);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}


void test_lorentz(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_lorentz" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,5);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;

  input.push_back(  x[0]*x[1] - x[0]*x[2] - x[3] + 1);
  input.push_back(  - x[0] + x[1]*x[2] - x[1]*x[3] + 1);
  input.push_back(  - x[0]*x[2] - x[1] + x[2]*x[3] + 1);
  input.push_back(  x[0]*x[3] - x[1]*x[3] - x[2] + 1);

  // ring r = 0,(x0, x1, x2, x3, x4),dp;
  // ideal i = x0*x1 - x0*x2 - x3 + 1, - x0 + x1*x2 - x1*x3 + 1, - x0*x2 - x1 + x2*x3 + 1, x0*x3 - x1*x3 - x2 + 1;
  // i = groebner(i);
  // i = homog(i, x4);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");

  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[4]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 4, 5));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 1, 1, 5);
  SetEntry(bettis, 1, 2, 1);
  SetEntry(bettis, 2, 2, 11);
  SetEntry(bettis, 2, 3, 10);
  SetEntry(bettis, 2, 4, 1);
  SetEntry(bettis, 3, 4, 1);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}

void test_mickey(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_mickey" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,3);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(   power(x[0],2)+ 4*power(x[1],2) - 4);
  input.push_back(   - x[0]+2*power(x[1],2));

  // ring r = 0,(x0, x1, x2),dp;
  // ideal i =x0^2+ 4*x1^2 - 4, - x0+2*x1^2;
  // i = groebner(i);
  // i = homog(i, x2);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");


  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[2]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  matrix bettis(NewDenseMat(RingZZ(), 3, 3));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 1, 1, 2);
  SetEntry(bettis, 2, 2, 1);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}

void test_lorentz_min_res(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_lorentz_min_res" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,5);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(  x[0]*x[1] - x[0]*x[2] - x[3] + 1);
  input.push_back(  - x[0] + x[1]*x[2] - x[1]*x[3] + 1);
  input.push_back(  - x[0]*x[2] - x[1] + x[2]*x[3] + 1);
  input.push_back(  x[0]*x[3] - x[1]*x[3] - x[2] + 1);

  // ring r = 0,(x0, x1, x2, x3, x4),dp;
  // ideal i = x0*x1 - x0*x2 - x3 + 1, - x0 + x1*x2 - x1*x3 + 1, - x0*x2 - x1 + x2*x3 + 1, x0*x3 - x1*x3 - x2 + 1;
  // i = groebner(i);
  // i = homog(i, x4);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");

  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[4]);
  }
  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  std::vector<matrix> maps(MinimalResolution(mill));
  // std::cout <<  "maps.size() = " << maps.size() << std::endl;
  // std::cout <<  "maps[0] = " << maps[0] << std::endl;
  // std::cout <<  "maps[1] = " << maps[1] << std::endl;
  // std::cout <<  "maps[2] = " << maps[2] << std::endl;
  // std::cout <<  "maps[3] = " << maps[3] << std::endl;
  matrix m = maps[0] * maps[1];
  for (long i = 0; i < NumRows(m); ++i)
  {
    for (long j = 0; j < NumCols(m); ++j)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(m(i,j)));
    }
  }
  m = maps[1] * maps[2];
  for (long i = 0; i < NumRows(m); ++i)
  {
    for (long j = 0; j < NumCols(m); ++j)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(m(i,j)));
    }
  }
  m = maps[2] * maps[3];
  for (long i = 0; i < NumRows(m); ++i)
  {
    for (long j = 0; j < NumCols(m); ++j)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(m(i,j)));
    }
  }
}


void test_lorentz_res(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_lorentz_res" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,5);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(  x[0]*x[1] - x[0]*x[2] - x[3] + 1);
  input.push_back(  - x[0] + x[1]*x[2] - x[1]*x[3] + 1);
  input.push_back(  - x[0]*x[2] - x[1] + x[2]*x[3] + 1);
  input.push_back(  x[0]*x[3] - x[1]*x[3] - x[2] + 1);

  // ring r = 0,(x0, x1, x2, x3, x4),dp;
  // ideal i = x0*x1 - x0*x2 - x3 + 1, - x0 + x1*x2 - x1*x3 + 1, - x0*x2 - x1 + x2*x3 + 1, x0*x3 - x1*x3 - x2 + 1;
  // i = groebner(i);
  // i = homog(i, x4);
  // i = groebner(i);
  // resolution rs = mres(i, 0);
  // print(betti(rs), "betti");
  //
  // R = QQ[x0, x1, x2, x3, x4]
  // I = ideal(x0*x1 - x0*x2 - x3 + 1, - x0 + x1*x2 - x1*x3 + 1, - x0*x2 - x1 + x2*x3 + 1, x0*x3 - x1*x3 - x2 + 1)
  // I = homogenize(I, x4)
  // C = res(I)
  // betti C
  input = JBMill(JBMill::Builder().setInput(input)).myReturnJB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[4]);
  }
  JBMill mill(JBMill::Builder().setInput(input));
  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  std::vector<matrix> maps(Resolution(mill));
  // std::cout <<  "maps.size() = " << maps.size() << std::endl;
  // std::cout <<  "maps[0] = " << maps[0] << std::endl;
  // std::cout <<  "maps[1] = " << maps[1] << std::endl;
  // std::cout <<  "maps[2] = " << maps[2] << std::endl;
  // std::cout <<  "maps[3] = " << maps[3] << std::endl;
  CoCoA_ASSERT_ALWAYS(maps.size() == 4);
  matrix m = maps[0] * maps[1];
  for (long i = 0; i < NumRows(m); ++i)
  {
    for (long j = 0; j < NumCols(m); ++j)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(m(i,j)));
    }
  }
  m = maps[1] * maps[2];
  for (long i = 0; i < NumRows(m); ++i)
  {
    for (long j = 0; j < NumCols(m); ++j)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(m(i,j)));
    }
  }
  m = maps[2] * maps[3];
  for (long i = 0; i < NumRows(m); ++i)
  {
    for (long j = 0; j < NumCols(m); ++j)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(m(i,j)));
    }
  }
}


void test_anna(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_anna" << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing_DMPI(Q,3);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> input;
  input.push_back(x[0]);
  input.push_back(x[1]);
  input.push_back(x[2] * x[2]);

  JBMill mill(JBMill::Builder().setInput(input));

  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);

  std::vector<matrix> maps(MinimalResolution(mill));
  CoCoA_ASSERT_ALWAYS(maps[0](0,0) == x[1]);
  CoCoA_ASSERT_ALWAYS(maps[0](0,1) == x[0]);
  CoCoA_ASSERT_ALWAYS(maps[0](0,2) == x[2] * x[2]);

  CoCoA_ASSERT_ALWAYS(maps[1](0,0) == x[0]);
  CoCoA_ASSERT_ALWAYS(maps[1](0,1) == -x[2] * x[2]);
  CoCoA_ASSERT_ALWAYS(maps[1](0,2) == 0);
  CoCoA_ASSERT_ALWAYS(maps[1](1,0) == -x[1]);
  CoCoA_ASSERT_ALWAYS(maps[1](1,1) == 0);
  CoCoA_ASSERT_ALWAYS(maps[1](1,2) == - x[2] * x[2]);
  CoCoA_ASSERT_ALWAYS(maps[1](2,0) == 0);
  CoCoA_ASSERT_ALWAYS(maps[1](2,1) == x[1]);
  CoCoA_ASSERT_ALWAYS(maps[1](2,2) == x[0]);

  CoCoA_ASSERT_ALWAYS(maps[2](0,0) == x[2] * x[2]);
  CoCoA_ASSERT_ALWAYS(maps[2](1,0) == x[0]);
  CoCoA_ASSERT_ALWAYS(maps[2](2,0) == -x[1]);
}

void test_BinomialIdealCycle(bool WriteName)
{
  if(WriteName)
  {
    std::cout << "test_BinomialIdealCycle" << std::endl;
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

  matrix bettis(NewDenseMat(RingZZ(), 3, 5));
  SetEntry(bettis, 0, 0, 1);
  SetEntry(bettis, 0, 1, 0);
  SetEntry(bettis, 0, 2, 0);
  SetEntry(bettis, 0, 3, 0);
  SetEntry(bettis, 0, 4, 0);

  SetEntry(bettis, 1, 0, 0);
  SetEntry(bettis, 1, 1, 4);
  SetEntry(bettis, 1, 2, 0);
  SetEntry(bettis, 1, 3, 0);
  SetEntry(bettis, 1, 4, 0);

  SetEntry(bettis, 2, 0, 0);
  SetEntry(bettis, 2, 1, 0);
  SetEntry(bettis, 2, 2, 9);
  SetEntry(bettis, 2, 3, 8);
  SetEntry(bettis, 2, 4, 2);
  CoCoA_ASSERT_ALWAYS(bettis == BettiDiagram(mill));
}



void program()
{
  GlobalManager CoCoAFoundations;
  bool WriteName(false);

  // test_myVectorBoolToVectorLong(WriteName);
  // test_myVariationWithoutRepetition(WriteName);
  // test_myPossibleWedges(WriteName);
  // test_myComputeGeneralBasis(WriteName);
  // test_myComputeBasicGraph(WriteName);
  test_myDirectMorseReduction1(WriteName);
  test_myDirectMorseReduction2(WriteName);
  test_myMapsAsMatrices(WriteName);
  test_myBettis(WriteName);
  test_onlyX(WriteName);
  test_onlyOne(WriteName);
  test_chandra4(WriteName);
  test_cassou(WriteName);
  test_boon(WriteName);
  test_caprasse(WriteName);
  test_conform1(WriteName);
  test_lorentz(WriteName);
  test_mickey(WriteName);
  test_lorentz_min_res(WriteName);
  test_lorentz_res(WriteName);
  test_anna(WriteName);
  test_BinomialIdealCycle(WriteName);
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
