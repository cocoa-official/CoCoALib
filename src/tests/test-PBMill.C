//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Orig author:  2015 Mario Albert

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
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/error.H"
#include "CoCoA/RingQQ.H"

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


void test_JBProjDim(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  CoCoA_ASSERT_ALWAYS(1 - 1 == mill.myProjDim());
}

void test_JBDepth(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  CoCoA_ASSERT_ALWAYS(1 == mill.myDepth());
}

void test_JBCls(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  RingElem elem = x[0] * x[0] + 1;
  //Test -- mill.myCls PPMonoidElem
  PPMonoidElem ElemPPM(LPP(elem));  // == x^2
  CoCoA_ASSERT_ALWAYS(0 == mill.myCls(ElemPPM));
  //Test -- mill.myCls for 1
  elem = one(polyRing);
  CoCoA_ASSERT_ALWAYS(-1 == mill.myCls(LPP(elem)));
}

void test_JBMinCls(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  CoCoA_ASSERT_ALWAYS(0 == mill.myMinCls());
}

void test_JBElementsWithClass(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  std::vector<RingElem> ClassZero(mill.myElementsWithClass(0));
  std::vector<RingElem> ClassOne(mill.myElementsWithClass(1));
  std::vector<RingElem> ClassTwo(mill.myElementsWithClass(2));
  CoCoA_ASSERT_ALWAYS(ClassZero.size() == 1);
  CoCoA_ASSERT_ALWAYS(ClassZero[0] == x[0]);
  CoCoA_ASSERT_ALWAYS(ClassOne.empty());
  CoCoA_ASSERT_ALWAYS(ClassTwo.empty());
}

void test_JBExtremalBettiNumbers(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  std::map<std::pair<long, long>, long> bettis(mill.myExtremalBettiNumbers());
  /* Expected Result:
   Deg
   1       1  --- x
   |       |      |
   0       0  --- 0
   |       |      |
   Module P0 --- I

   */
  CoCoA_ASSERT_ALWAYS(bettis.size() == 1);
  CoCoA_ASSERT_ALWAYS((*bettis.begin()).first.first == 0);
  CoCoA_ASSERT_ALWAYS((*bettis.begin()).first.second == 1);
  CoCoA_ASSERT_ALWAYS((*bettis.begin()).second == 1);
}


void test_JBRegSeq(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  vector<RingElem> RegSeq(mill.myRegSeq());
  CoCoA_ASSERT_ALWAYS(RegSeq.size() == 1);
  CoCoA_ASSERT_ALWAYS(RegSeq[0] == x[0]);
}

void test_JBMaxStronglyIndependentSet(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  vector<RingElem> MaxIndepent(mill.myMaxStronglyIndependentSet());
  CoCoA_ASSERT_ALWAYS(MaxIndepent.size() == 0);
}

void test_JBIsCohenMacualay(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  CoCoA_ASSERT_ALWAYS(false == mill.IamCohenMacaulay());
}

void test_JBRegularity(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  CoCoA_ASSERT_ALWAYS(1 == mill.myRegularity());
}

void test_JBSatiety(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  CoCoA_ASSERT_ALWAYS(1 == mill.mySatiety());
}

void test_JBSaturation(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill jbmill(JBMill::Builder().setInput(pols));
  PBMill mill(PBMill::Converter().setJBMill(jbmill));

  vector<RingElem> saturation(mill.mySaturation());
  CoCoA_ASSERT_ALWAYS(1 == saturation.size());
  CoCoA_ASSERT_ALWAYS(IsOne(saturation[0]));
}



void program()
{
  GlobalManager CoCoAFoundations;
//  test_stuff();

  test_JBProjDim();
  test_JBDepth();
  test_JBCls();
  test_JBMinCls();
  test_JBElementsWithClass();
  test_JBExtremalBettiNumbers();
  test_JBRegSeq();
  test_JBMaxStronglyIndependentSet();
  test_JBIsCohenMacualay();
  test_JBRegularity();
  test_JBSatiety();
  test_JBSaturation();

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
