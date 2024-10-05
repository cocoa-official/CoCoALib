//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Original author:  2015 Mario Albert

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
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/error.H"
#include "CoCoA/RingQQ.H"

using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

const bool WRITE_VECS = false;
const bool WRITE_NAME = false;


  bool operator()(const RingElem& r1, const RingElem& r2) {
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

void test_computation_empty_input() {
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
   // test empty input
  try {
    JBMill(JBMill::Builder().setInput(pols));
  } catch (const CoCoA::ErrorInfo& err) {
    CoCoA_ASSERT_ALWAYS(err == ERR::nonstandard);
  }
}

void test_computation_zero_input() {
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  // test zero input
  pols.push_back(zero(polyRing));
  try {
    JBMill(JBMill::Builder().setStrategy(TQDegree).setInput(pols));
  } catch (const CoCoA::ErrorInfo& err) {
    CoCoA_ASSERT_ALWAYS(err == ERR::nonstandard);
  }
  try {
    JBMill(JBMill::Builder().setStrategy(TQBlockHigh).setInput(pols));
  } catch (const CoCoA::ErrorInfo& err) {
    CoCoA_ASSERT_ALWAYS(err == ERR::nonstandard);
  }
  try {
    JBMill(JBMill::Builder().setStrategy(TQBlockLow).setInput(pols));
  } catch (const CoCoA::ErrorInfo& err) {
    CoCoA_ASSERT_ALWAYS(err == ERR::nonstandard);
  }
  try {
    JBMill(JBMill::Builder().setStrategy(GBCompletion).setInput(pols));
  } catch (const CoCoA::ErrorInfo& err) {
    CoCoA_ASSERT_ALWAYS(err == ERR::nonstandard);
  }
}

void test_computation_only_one() {
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(one(polyRing));
  result = JBMill(JBMill::Builder().setStrategy(TQDegree).setInput(pols)).myReturnJB();
  test_one(result);
  result = JBMill(JBMill::Builder().setStrategy(TQBlockLow).setInput(pols)).myReturnJB();
  test_one(result);
  result = JBMill(JBMill::Builder().setStrategy(TQBlockHigh).setInput(pols)).myReturnJB();
  test_one(result);
  result = JBMill(JBMill::Builder().setStrategy(GBCompletion).setInput(pols)).myReturnJB();
  test_one(result);
}

void test_computation_only_x() {
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  result = JBMill(JBMill::Builder().setStrategy(TQDegree).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
  result = JBMill(JBMill::Builder().setStrategy(TQBlockLow).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
  result = JBMill(JBMill::Builder().setStrategy(TQBlockHigh).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
  result = JBMill(JBMill::Builder().setStrategy(GBCompletion).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
}

void test_computation_ex1(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  pols.push_back(x[0] * x[0]);
  pols.push_back(x[0] * x[0]);
  result = JBMill(JBMill::Builder().setStrategy(TQDegree).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
  result = JBMill(JBMill::Builder().setStrategy(TQBlockLow).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
  result = JBMill(JBMill::Builder().setStrategy(TQBlockHigh).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
  result = JBMill(JBMill::Builder().setStrategy(GBCompletion).setInput(pols)).myReturnJB();
  test_principal_ideal(result, x[0]);
}

void test_myNonMultVar(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  std::map<PPMonoidElem, std::vector<bool> > ResJBNonMultVar = mill.myNonMultVars();
  CoCoA_ASSERT_ALWAYS(1 == ResJBNonMultVar.size());
  std::vector<bool> NonMultVarVector = (*ResJBNonMultVar.begin()).second;
  CoCoA_ASSERT_ALWAYS(1 == NonMultVarVector.size());
  CoCoA_ASSERT_ALWAYS(NonMultVarVector[0] == false);
}

void test_myMultVar(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  std::map<PPMonoidElem, std::vector<bool> > ResJBMultVar = mill.myMultVars();
  CoCoA_ASSERT_ALWAYS(1 == ResJBMultVar.size());
  std::vector<bool> MultVarVector = (*ResJBMultVar.begin()).second;
  CoCoA_ASSERT_ALWAYS(1 == MultVarVector.size());
  CoCoA_ASSERT_ALWAYS(MultVarVector[0] == true);
}

void test_IAmPommaretBasis(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  CoCoA_ASSERT_ALWAYS(mill.IamPommaretBasis() == true);
}

void test_IAmHomogeneous(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  CoCoA_ASSERT_ALWAYS(mill.IamHomogenous() == true);
}

void test_IAmMonomialIdeal(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  CoCoA_ASSERT_ALWAYS(mill.IamMonomialIdeal() == true);
}

void test_myHilbertPol(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill1(JBMill::Builder().setInput(pols));

  SparsePolyRing HilbPolPolyRing = NewPolyRing(Q, SymbolRange("s", 0, 0), NewStdDegLexOrdering(1));
  const std::vector<RingElem> s = indets(HilbPolPolyRing);

  CoCoA_ASSERT_ALWAYS(IsZero(mill1.myHilbertPol(s[0])));

  pols.clear();
  SparsePolyRing P2 = NewPolyRing(Q, 2);
  x = indets(P2);

  pols.push_back(x[0] * x[0]);
  pols.push_back(x[0] * x[1]);

  JBMill mill2(JBMill::Builder().setInput(pols));

  CoCoA_ASSERT_ALWAYS(1 == mill2.myHilbertPol(s[0]));

  pols.clear();
  SparsePolyRing P4 = NewPolyRing(Q, 4);
  x = indets(P4);

  pols.push_back(power(x[0],3)*power(x[1],2)+4*power(x[0],2)*power(x[1],2)*x[2]+288*power(x[0],2)*power(x[1],2)-power(x[0],2)*x[1]*power(x[2],2)+207*power(x[0],2)*x[1]*x[2]-3456*power(x[0],2)*x[1]+1152*x[0]*power(x[1],2)*x[2]+20736*x[0]*power(x[1],2)+156*x[0]*x[1]*power(x[2],2)+19008*x[0]*x[1]*x[2]-497664*x[0]*x[1]+x[0]*power(x[2],3)+432*x[0]*power(x[2],2)+62208*x[0]*x[2]+2985984*x[0]+82944*power(x[1],2)*x[2]);
  pols.push_back(power(x[3],3)*power(x[1],3)+4*power(x[3],3)*power(x[1],2)+4*power(x[3],2)*power(x[1],3)-power(x[3],2)*power(x[1],2)*x[2]-48*power(x[3],2)*power(x[1],2)-5*power(x[3],2)*x[1]*x[2]+108*x[3]*x[1]*x[2]+x[3]*power(x[2],2)+144*x[3]*x[2]-1728*x[2]);
  pols.push_back(4*power(x[3],2)*x[0]*power(x[2],2)+1152*power(x[3],2)*x[0]*x[2]+82944*power(x[3],2)*x[0]+power(x[3],2)*power(x[2],3)+288*power(x[3],2)*power(x[2],2)+20736*power(x[3],2)*x[2]-x[3]*power(x[0],2)*power(x[2],2)+156*x[3]*power(x[0],2)*x[2]+207*x[3]*x[0]*power(x[2],2)+19008*x[3]*x[0]*x[2]-3456*x[3]*power(x[2],2)-497664*x[3]*x[2]+power(x[0],3)*x[2]+432*power(x[0],2)*x[2]+62208*x[0]*x[2]+2985984*x[2]);
  pols.push_back(power(x[3],3)*power(x[1],3)+4*power(x[3],3)*power(x[1],2)-power(x[3],2)*x[0]*power(x[1],2)+4*power(x[3],2)*power(x[1],3)-48*power(x[3],2)*power(x[1],2)-5*x[3]*x[0]*power(x[1],2)+108*x[3]*x[0]*x[1]+power(x[0],2)*x[1]+144*x[0]*x[1]-1728*x[0]);

  JBMill mill4(JBMill::Builder().setInput(pols));
  CoCoA_ASSERT_ALWAYS( (8 * s[0] + 2) == mill4.myHilbertPol(s[0]));
}


void test_myHilbertFunc(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  CoCoA_ASSERT_ALWAYS(mill.myHilbertFunc(BigInt(0)) == BigInt(1));
  CoCoA_ASSERT_ALWAYS(mill.myHilbertFunc(BigInt(1)) == BigInt(0));
  CoCoA_ASSERT_ALWAYS(mill.myHilbertFunc(BigInt(2)) == BigInt(0));
  CoCoA_ASSERT_ALWAYS(mill.myHilbertFunc(BigInt(3)) == BigInt(0));
}

void test_myHilbertSeries(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));
  FractionField QR = NewFractionField(polyRing);
  RingHom embeddedHom(EmbeddingHom(QR));
  RingElem qx = embeddedHom(x[0]);
  CoCoA_ASSERT_ALWAYS(mill.myHilbertSeries(qx) == qx / (1 - qx));
}

void test_myStandardRepresentation(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));

  RingElem elem(x[0] * x[0] + x[0] + 1);
  std::pair<std::map<PPMonoidElem, RingElem>, RingElem> res(mill.myStandardRepresentation(elem));
  CoCoA_ASSERT_ALWAYS(res.second == one(polyRing));
  CoCoA_ASSERT_ALWAYS((*res.first.begin()).second == x[0] + 1);

  RingElem JBSRZero(zero(polyRing));
  res = mill.myStandardRepresentation(JBSRZero);
  CoCoA_ASSERT_ALWAYS(res.second == zero(polyRing));
  CoCoA_ASSERT_ALWAYS((*res.first.begin()).second == zero(polyRing));
}

void test_myNormalForm(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));

  RingElem elem(x[0] * x[0] + x[0] + 1);
  CoCoA_ASSERT_ALWAYS(one(polyRing) == mill.myJNormalForm(elem));

  RingElem JBSRZero(zero(polyRing));
  CoCoA_ASSERT_ALWAYS(zero(polyRing) == mill.myJNormalForm(JBSRZero));
}

void test_mySyzygy(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));

  FGModule firstSyz = mill.mySyzygy();
  vector<ModuleElem> gensFirstSyz = gens(firstSyz);
  CoCoA_ASSERT_ALWAYS(0 == gensFirstSyz.size());
}

void test_myDim(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing P1 = NewPolyRing(Q, 1);
  vector<RingElem> x = indets(P1);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill mill1(JBMill::Builder().setInput(pols));

  CoCoA_ASSERT_ALWAYS(0 == mill1.myDim());

  pols.clear();
  SparsePolyRing P2 = NewPolyRing(Q, 2);
  x = indets(P2);

  pols.push_back(x[0] * x[0]);
  pols.push_back(x[0] * x[1]);

  JBMill mill2(JBMill::Builder().setInput(pols));

  CoCoA_ASSERT_ALWAYS(0 == mill2.myDim());

  // Cohn 2
  pols.clear();
  SparsePolyRing P4 = NewPolyRing(Q, 4);
  x = indets(P4);

  pols.push_back(power(x[0],3)*power(x[1],2)+4*power(x[0],2)*power(x[1],2)*x[2]+288*power(x[0],2)*power(x[1],2)-power(x[0],2)*x[1]*power(x[2],2)+207*power(x[0],2)*x[1]*x[2]-3456*power(x[0],2)*x[1]+1152*x[0]*power(x[1],2)*x[2]+20736*x[0]*power(x[1],2)+156*x[0]*x[1]*power(x[2],2)+19008*x[0]*x[1]*x[2]-497664*x[0]*x[1]+x[0]*power(x[2],3)+432*x[0]*power(x[2],2)+62208*x[0]*x[2]+2985984*x[0]+82944*power(x[1],2)*x[2]);
  pols.push_back(power(x[3],3)*power(x[1],3)+4*power(x[3],3)*power(x[1],2)+4*power(x[3],2)*power(x[1],3)-power(x[3],2)*power(x[1],2)*x[2]-48*power(x[3],2)*power(x[1],2)-5*power(x[3],2)*x[1]*x[2]+108*x[3]*x[1]*x[2]+x[3]*power(x[2],2)+144*x[3]*x[2]-1728*x[2]);
  pols.push_back(4*power(x[3],2)*x[0]*power(x[2],2)+1152*power(x[3],2)*x[0]*x[2]+82944*power(x[3],2)*x[0]+power(x[3],2)*power(x[2],3)+288*power(x[3],2)*power(x[2],2)+20736*power(x[3],2)*x[2]-x[3]*power(x[0],2)*power(x[2],2)+156*x[3]*power(x[0],2)*x[2]+207*x[3]*x[0]*power(x[2],2)+19008*x[3]*x[0]*x[2]-3456*x[3]*power(x[2],2)-497664*x[3]*x[2]+power(x[0],3)*x[2]+432*power(x[0],2)*x[2]+62208*x[0]*x[2]+2985984*x[2]);
  pols.push_back(power(x[3],3)*power(x[1],3)+4*power(x[3],3)*power(x[1],2)-power(x[3],2)*x[0]*power(x[1],2)+4*power(x[3],2)*power(x[1],3)-48*power(x[3],2)*power(x[1],2)-5*x[3]*x[0]*power(x[1],2)+108*x[3]*x[0]*x[1]+power(x[0],2)*x[1]+144*x[0]*x[1]-1728*x[0]);

  JBMill mill4(JBMill::Builder().setInput(pols));

  CoCoA_ASSERT_ALWAYS(1 == mill4.myDim());

}

void test_myComplementaryDecomposition(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));

  vector<std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp(mill.myComplementaryDecomposition());
  CoCoA_ASSERT_ALWAYS(1 == CompDecomp.size());
  CoCoA_ASSERT_ALWAYS(IsOne(CompDecomp[0].first));
  CoCoA_ASSERT_ALWAYS(1 == CompDecomp[0].second.size());
  CoCoA_ASSERT_ALWAYS(false == CompDecomp[0].second[0]);
}

void test_myStandardPairs(){
  if(WRITE_NAME)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(x[0]);
  JBMill mill(JBMill::Builder().setInput(pols));

  vector<std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp(mill.myComplementaryDecomposition());
  CoCoA_ASSERT_ALWAYS(1 == CompDecomp.size());
  CoCoA_ASSERT_ALWAYS(IsOne(CompDecomp[0].first));
  CoCoA_ASSERT_ALWAYS(1 == CompDecomp[0].second.size());
  CoCoA_ASSERT_ALWAYS(false == CompDecomp[0].second[0]);
}

void test_TrivialJB(){
  if(WRITE_NAME){
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 1);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  pols.push_back(one(polyRing));
  JBMill OneMill(JBMill::Builder().setInput(pols));
  // Test -- Correct Janet Basis
  std::vector<RingElem> ReturnedJB(OneMill.myReturnJB());
  CoCoA_ASSERT_ALWAYS(1 == ReturnedJB.size());
  CoCoA_ASSERT_ALWAYS(IsOne(ReturnedJB[0]));
  std::map<PPMonoidElem, std::vector<bool> > ResJBNonMultVar = OneMill.myNonMultVars();
  CoCoA_ASSERT_ALWAYS(1 == ResJBNonMultVar.size());
  std::vector<bool> NonMultVarVector = (*ResJBNonMultVar.begin()).second;
  CoCoA_ASSERT_ALWAYS(1 == NonMultVarVector.size());
  CoCoA_ASSERT_ALWAYS(NonMultVarVector[0] == false);
   // Test -- JBMultVar
  std::map<PPMonoidElem, std::vector<bool> > ResJBMultVar = OneMill.myMultVars();
  CoCoA_ASSERT_ALWAYS(1 == ResJBMultVar.size());
  std::vector<bool> MultVarVector = (*ResJBMultVar.begin()).second;
  CoCoA_ASSERT_ALWAYS(1 == MultVarVector.size());
  CoCoA_ASSERT_ALWAYS(MultVarVector[0] == true);
  //Test -- JBIsPommaretBasis
  CoCoA_ASSERT_ALWAYS(OneMill.IamPommaretBasis() == true);
  //Test -- JBIsHomogenous
  CoCoA_ASSERT_ALWAYS(OneMill.IamHomogenous() == true);
  //Test -- JBIsMonomialIdeal
  CoCoA_ASSERT_ALWAYS(OneMill.IamMonomialIdeal() == true);
  //Test -- JBHilbertPol
  SparsePolyRing HilbPolPolyRing = NewPolyRing(Q, SymbolRange("s", 0, 0), NewStdDegLexOrdering(1));
  const std::vector<RingElem> s = indets(HilbPolPolyRing);
  CoCoA_ASSERT_ALWAYS(IsZero(OneMill.myHilbertPol(s[0])));
  //Test -- JBHilbertFunc
  CoCoA_ASSERT_ALWAYS(0 == OneMill.myHilbertFunc(BigInt(0)));
  CoCoA_ASSERT_ALWAYS(0 == OneMill.myHilbertFunc(BigInt(1)));
  CoCoA_ASSERT_ALWAYS(0 == OneMill.myHilbertFunc(BigInt(2)));
  CoCoA_ASSERT_ALWAYS(0 == OneMill.myHilbertFunc(BigInt(3)));
  //Test -- JBHilbertSeries
  FractionField QR = NewFractionField(polyRing);
  RingHom embeddedHom(EmbeddingHom(QR));
  RingElem qx = embeddedHom(x[0]);
  OneMill.myHilbertSeries(qx);
  CoCoA_ASSERT_ALWAYS(OneMill.myHilbertSeries(qx) == 1 / (1 - qx));
  //Test -- JBStandardRepresentation
  RingElem elem = x[0] * x[0] + x[0] + 1;
  std::pair<std::map<PPMonoidElem, RingElem>, RingElem> res = OneMill.myStandardRepresentation(elem);
  CoCoA_ASSERT_ALWAYS(IsZero(res.second));
  CoCoA_ASSERT_ALWAYS((*res.first.begin()).second == elem);
  //Test -- JBNormalForm
  CoCoA_ASSERT_ALWAYS(IsZero(OneMill.myJNormalForm(elem)));
  //Test -- JBSyzygy
  FGModule firstSyz = OneMill.mySyzygy();
  std::vector<ModuleElem> gensFirstSyz = gens(firstSyz);
  CoCoA_ASSERT_ALWAYS(0 == gensFirstSyz.size());
  //Test -- JBDim
  CoCoA_ASSERT_ALWAYS(0 == OneMill.myDim());
  //Test -- JBMinCls
  CoCoA_ASSERT_ALWAYS(-1 == OneMill.myMinCls());
  //Test -- JBElementsWithClass
  CoCoA_ASSERT_ALWAYS(OneMill.myElementsWithClass(0).empty());
  CoCoA_ASSERT_ALWAYS(OneMill.myElementsWithClass(1).empty());
  CoCoA_ASSERT_ALWAYS(OneMill.myElementsWithClass(2).empty());
  //Test -- JBComplementaryDecomposition
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > CompDecomp = OneMill.myComplementaryDecomposition();
  CoCoA_ASSERT_ALWAYS(0 == CompDecomp.size());

  //Test -- JBStandardPairs
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > StdPairs = OneMill.myStandardPairs();
  CoCoA_ASSERT_ALWAYS(0 == StdPairs.size());
}


void test_myGBCompletionStrategy()
{
  if(WRITE_NAME){
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 2);
  const vector<RingElem> x = indets(polyRing);
  vector<RingElem> pols;
  vector<RingElem> result;
  pols.push_back(x[0] * x[0]);
  pols.push_back(x[1]);

  result = JBMill(JBMill::Builder().setInput(pols).setStrategy(GBCompletion)).myReturnGB();


  vector<RingElem> expectedResult;
  expectedResult.push_back(x[1]);
  expectedResult.push_back(x[0] * x[0]);

  test_equal_vectors(expectedResult, result);

  expectedResult.clear();

  expectedResult.push_back(x[1]);
  expectedResult.push_back(x[0] * x[1]);
  expectedResult.push_back(x[0] * x[0]);

  result = JBMill(JBMill::Builder().setInput(pols).setStrategy(GBCompletion)).myReturnJB();

  test_equal_vectors(expectedResult, result);
}

void test_myComputation_solotarev()
{
  if(WRITE_NAME){
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 8);
  vector<RingElem> x(indets(polyRing));
  vector<RingElem> input;

  input.push_back(  - 2*x[3]*x[7] - x[4]*x[6]);
  input.push_back(  4*x[1] + 9*x[4]);
  input.push_back(  - 4*x[2]*x[7] - 3*x[3]*x[6] - 2*x[4]*x[5]);
  input.push_back(  9*x[0] - 7*x[2] - 8*x[5]);
  input.push_back(  - 5*x[2]*x[6] - 4*x[3]*x[5] - 3*x[4] - 6*x[7]);
  input.push_back(  9*x[1] - 6*x[2]*x[5] - 5*x[3] - 7*x[6]);
  input.push_back(  6*x[0] - 5*x[1] + 9*x[3]);
  input.push_back(  - 7*x[0] + 9*x[2] + 8);

  // test whether all Janet Basis are equal

  vector<RingElem> tqDegree(JBMill(JBMill::Builder().setInput(input).setStrategy(TQDegree)).myReturnJB());
  vector<RingElem> tqBlockLow(JBMill(JBMill::Builder().setInput(input).setStrategy(TQBlockLow)).myReturnJB());
  vector<RingElem> tqBlockHigh(JBMill(JBMill::Builder().setInput(input).setStrategy(TQBlockHigh)).myReturnJB());

  test_equal_vectors(tqDegree, tqBlockLow);
  test_equal_vectors(tqDegree, tqBlockHigh);

  // test whether all GB are qual

  tqDegree    = JBMill(JBMill::Builder().setInput(input).setStrategy(TQDegree)).myReturnGB();
  tqBlockLow  = JBMill(JBMill::Builder().setInput(input).setStrategy(TQBlockLow)).myReturnGB();
  tqBlockHigh = JBMill(JBMill::Builder().setInput(input).setStrategy(TQBlockHigh)).myReturnGB();
  std::vector<RingElem> gb;
  std::vector<RingElem> minGens;

  ComputeGBasis(gb, minGens, input);

  if (WRITE_VECS){
    cout << "gb <-> tqdegree" << endl;
  }
  test_equal_vectors(gb, tqDegree);
  if (WRITE_VECS){
    cout << "tqdegree<-> tqblocklow" << endl;
  }
  test_equal_vectors(tqDegree, tqBlockLow);
  if (WRITE_VECS){
    cout << "tqdegree <-> tqblockhigh" << endl;
  }
  test_equal_vectors(tqDegree, tqBlockHigh);
}

void test_myCompletion_braid411()
{
  if(WRITE_NAME){
    std::cout << __FUNCTION__ << std::endl;
  }
  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 8);
  vector<RingElem> x(indets(polyRing));
  vector<RingElem> input;
  input.push_back(x[0] * x[1] * x[1] - x[1] * x[2] * x[2]);
  input.push_back(x[0] * x[2] * x[2] - x[1] * x[2] * x[3]);
  input.push_back(power(x[0], 3) + power(x[1], 3) + x[0] * x[1] * x[2] + x[2] * x[2] * x[2]);

  // computing GBbasis
  std::vector<RingElem> gb;
  std::vector<RingElem> minGens;
  ComputeGBasis(gb, minGens, input);

  vector<RingElem> completion(JBMill(JBMill::Builder().setInput(input).setStrategy(GBCompletion)).myReturnGB());

  test_equal_vectors(gb, completion);
}

void test_jdivisorofidealone()
{
  if(WRITE_NAME){
    std::cout << __FUNCTION__ << std::endl;
  }

  ring Q = RingQQ();
  SparsePolyRing polyRing = NewPolyRing(Q, 8);
  vector<RingElem> x(indets(polyRing));
  vector<RingElem> input;
  input.push_back(one(polyRing));

  // computing GBbasis
  std::vector<RingElem> gb;
  std::vector<RingElem> minGens;
  ComputeGBasis(gb, minGens, input);

  JBMill mill(JBMill(JBMill::Builder().setInput(input)));

  CoCoA_ASSERT_ALWAYS(IsOne(mill.myJDivisor(one(polyRing))));
  CoCoA_ASSERT_ALWAYS(IsOne(mill.myJDivisor(x[0])));
  CoCoA_ASSERT_ALWAYS(IsZero(mill.myJDivisor(zero(polyRing))));
}

void program()
{
  GlobalManager CoCoAFoundations;
  test_computation_empty_input();
  test_computation_zero_input();
  test_computation_only_one();
  test_computation_only_x();
  test_computation_only_x();
  test_computation_ex1();
  test_myNonMultVar();
  test_myMultVar();
  test_IAmPommaretBasis();
  test_IAmHomogeneous();
  test_IAmMonomialIdeal();
  test_myHilbertPol();
  test_myHilbertFunc();
  test_myHilbertSeries();
  test_myStandardRepresentation();
  test_myNormalForm();
  test_mySyzygy();
  test_myDim();
  test_myComplementaryDecomposition();
  test_myStandardPairs();
  test_TrivialJB();
  test_myGBCompletionStrategy();
  test_myComputation_solotarev();
  test_myCompletion_braid411();
  test_jdivisorofidealone();
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
