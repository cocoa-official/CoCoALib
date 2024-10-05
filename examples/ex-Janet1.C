// Copyright (c) 2013-2015  John Abbott,  Anna M Bigatti
// Orig author: 2013-2015  Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
#include <bitset>
using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

//----------------------------------------------------------------------
const string ShortDescription =
  "In this file we explain how to compute the Janet Basis in CoCoA.  \n";

const string LongDescription =
  " We explain the function JanetBasis and show whicht options we can choose. \n";
//----------------------------------------------------------------------

void output(vector<RingElem> vec)
{
  for(vector<RingElem>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
    {
      cout << *iter << endl;
    }
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  ring Q = RingQQ();

  cout <<"/////////////////////////////////////Cyclic 5//////////////////////////////////////////" << endl;
  SparsePolyRing polyRingCyc5 = NewPolyRing(Q, SymbolRange("x",0,4), StdDegLex);
  const vector<RingElem> x = indets(polyRingCyc5);
  vector<RingElem>cyclic5;
  cyclic5.push_back(x[0]+x[1]+x[2]+x[3]+x[4]);
  cyclic5.push_back(x[0]*x[1] + x[2]*x[1] + x[2]*x[3] + x[3]*x[4] + x[0]*x[4]);
  cyclic5.push_back(x[1]*x[2]*x[0] + x[2]*x[3]*x[1] + x[4]*x[2]*x[3] + x[0]*x[4]*x[1] + x[0]*x[4]*x[3]);
  cyclic5.push_back(x[0]*x[1]*x[2]*x[3]+x[0]*x[1]*x[2]*x[4]+x[0]*x[1]*x[3]*x[4]+x[0]*x[2]*x[3]*x[4]+x[2]*x[3]*x[1]*x[4]);
  cyclic5.push_back(x[0]*x[1]*x[2]*x[3]*x[4] -1 );
  ideal I = ideal(cyclic5);
  double timing = CpuTime();
  vector<RingElem> GroebnerBasis = TidyGens(I);
  cout << "Time taken (TinyGens) is " << CpuTime() - timing << endl;
  cout << "size of GroebnerBasis = " << GroebnerBasis.size() << endl << endl;


  cout << "---------computing the Janet Base with a list of polynomials(TQBlockLow):------------------------" << endl;
  cout << "-----------------------------first two criteria ---------------------------------" << endl;
  // the following line computes the Janet basis
  Involutive::JBMill mill1(Involutive::JBMill::Builder().setInput(cyclic5)
                                                        .setStrategy(Involutive::TQBlockLow)
                                                        .setInvolutiveCriteria(bitset<3>(string("110"))));
  cout << "-----------------------------Output GB ---------------------------------" << endl;
  vector<RingElem> result(mill1.myReturnGB());
  std::cout << "len(result) = " << len(result) << std::endl;

  cout << "-----------------------------Output JB---------------------------------" << endl;
  result = mill1.myReturnJB();
  std::cout << "len(result) = " << len(result) << std::endl;

  cout << "-----------------------------first + third criteria -----------------------------" << endl;
  Involutive::JBMill::Builder builder;
  builder.setInput(cyclic5);
  builder.setStrategy(Involutive::TQBlockLow).setInvolutiveCriteria(bitset<3>(string("101")));

  Involutive::JBMill mill2(builder);
  std::cout <<  "len(mill2.myReturnJB()) = " << len(mill2.myReturnJB()) << std::endl;

  cout << "---------computing the Janet Base with a list of polynomials (TQBlockHigh):------------------------" << endl;
  // reusing builder object above
  builder.setStrategy(Involutive::TQBlockHigh);
  builder.setInvolutiveCriteria(bitset<3>(string("110")));
  // the following line computes the Janet basis
  Involutive::JBMill mill3(builder);

  std::cout <<  "len(mill3.myReturnGB()) = " << len(mill3.myReturnGB()) << std::endl;

  cout << "//////////////////////Example over Fp7//////////////////////////////////////////////" << endl;
  ring Fp7 = NewZZmod(7);
  SparsePolyRing polyRing = NewPolyRing(Q, SymbolRange("x",0,2), StdDegLex);
  const vector<RingElem> y = indets(polyRing);
  vector<RingElem> gEx;
  gEx.push_back(2*y[0]*y[1]*power(y[2],3) - power(y[0],4));
  gEx.push_back(7*power(y[1],5)*y[2] - 6*power(y[1],4));
  gEx.push_back(4*power(y[1],6));
  gEx.push_back(2*power(y[0],4)*y[1]*y[1]);

  cout << "---------computing the Janet Base with a list of polynomials(TQBlockLow):------------------------" << endl;
  cout << "-----------------------------first two criteria ---------------------------------" << endl;
  // the following line computes the Janet basis
  Involutive::JBMill mill4(Involutive::JBMill::Builder().setInput(gEx)
                                                        .setStrategy(Involutive::TQBlockLow)
                                                        .setInvolutiveCriteria(bitset<3>(string("110"))));
  cout << "-----------------------------Length GB ---------------------------------" << endl;
  std::cout << "len(mill4.myReturnGB()) = " << len(mill4.myReturnGB()) << std::endl;

  cout << "-----------------------------Length JB---------------------------------" << endl;
  std::cout << "len(mill4.myReturnJB()) = " << len(mill4.myReturnJB()) << std::endl;

  cout << "-----------------------------Output myReturnJB---------------------------------" << endl;
  output(mill4.myReturnJB());
  cout << "-----------------------------Output myReturnGB---------------------------------" << endl;
  output(mill4.myReturnGB());

  cout << "------------------------------TidyGens-----------------------------------------------" << endl;
  ideal gIdeal(gEx);
  vector<RingElem> result2 = TidyGens(gIdeal);
  std::cout << "len(result2) = " << len(result2) << std::endl << endl;

  cout << "-----------------------------Output Tidy Gens----------------------------------------" << endl;
  output(result2);
}

//----------------------------------------------------------------------
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
  return 1;
}
