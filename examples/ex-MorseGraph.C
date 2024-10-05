// Copyright (c) 2013-2015  John Abbott,  Anna Bigatti
// Orig author: 2013-2015  Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

//----------------------------------------------------------------------
const string ShortDescription =
  "In this file we explain how to compute a free Resolution via Pommaret Basis in CoCoA.  \n";

const string LongDescription =
  " JBResolution, JBBettiDiagramm, JBMinimalResolution. \n";
//----------------------------------------------------------------------



void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  const ring& QQ = RingQQ();
  SparsePolyRing polyRing2 = NewPolyRing_DMPI(QQ, 4);
  vector<RingElem> x = indets(polyRing2);
  vector<RingElem> input;
  input.push_back(   power(x[0],2) + 4*power(x[1],2) - 4);
  input.push_back(   - x[0] + 2*power(x[1],2));


  //////////////////////////////////////////////////////////////////
  std::cout << "Resolution of ideal" << std::endl;

  JBMill mill1(Involutive::JBMill::Builder().setInput(input));
  std::vector<matrix> res1(Involutive::Resolution(mill1));
  long LenRes = len(res1);
  for (long i=0; i < LenRes; ++i)
  {
    std::cout << "F[" << i << "] == " << res1[i] << std::endl;
  }


  //////////////////////////////////////////////////////////////////
  std::cout << "Betti Diagram of homogenized Ideal" << std::endl;

  input = mill1.myReturnGB();
  for (std::vector<RingElem>::iterator i = input.begin(); i != input.end(); ++i)
  {
    *i = homog(*i, x[2]);
  }
  JBMill mill2(Involutive::JBMill::Builder().setInput(input));
  std::cout << Involutive::BettiDiagram(mill2) << std::endl;

  //////////////////////////////////////////////////////////////////
  std::cout << "Second Betti Column of homogenized Ideal" << std::endl;

  std::cout << Involutive::BettiColumn(mill2, 2) << std::endl;

  //////////////////////////////////////////////////////////////////
  std::cout << "Betti Number (2,2) of homogenized Ideal" << std::endl;

  std::cout << Involutive::BettiNumber(mill2, 2, 2) << std::endl;

  //////////////////////////////////////////////////////////////////
  std::cout << "Minimal Free Resolution of homogenized Ideal" << std::endl;

  std::vector<matrix> res2(Involutive::MinimalResolution(mill2));
  for (long i = 0; i < len(res2); ++i)
  {
    std::cout << "F[" << i << "] == " << res2[i] << std::endl;
  }

  ///////////// Janet example
  input.clear();
  input.push_back(x[0] * x[0] * x[2]);
  input.push_back(x[0] * x[0] * x[3]);
  input.push_back(x[0] * x[2] - x[2] * x[3]);
  input.push_back(x[0] * x[3]);
  input.push_back(x[1] * x[3]);


  //////////////////////////////////////////////////////////////////
  std::cout << "Resolution of janet ideal" << std::endl;

  JBMill mill3(Involutive::JBMill::Builder().setInput(input));
  std::cout << mill3.myReturnJB() << std::endl;
  std::vector<matrix> res3(Involutive::Resolution(mill3));
  LenRes = len(res3);
  for (long i=0; i < LenRes; ++i)
  {
    std::cout << "F[" << i << "] == " << res3[i] << std::endl;
  }
  std::cout <<  "res3[0] * res3[1] = " << res3[0] * res3[1] << std::endl;
  std::cout <<  "res3[1] * res3[2] = " << res3[1] * res3[2] << std::endl;

  std::cout << "Minimal Resolution of janet ideal" << std::endl;

  std::vector<matrix> res4(Involutive::MinimalResolution(mill3));
  LenRes = len(res4);
  for (long i=0; i < LenRes; ++i)
  {
    std::cout << "F[" << i << "] == " << res4[i] << std::endl;
  }
  std::cout <<  "res4[0] * res4[1] = " << res4[0] * res4[1] << std::endl;
  std::cout <<  "res4[1] * res4[2] = " << res4[1] * res4[2] << std::endl;
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
