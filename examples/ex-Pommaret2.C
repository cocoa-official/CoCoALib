// Copyright (c) 2013-2015  John Abbott,  Anna M. Bigatti
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
  "In this file we explain how to using some methods which are using the Pommaret Basis.  \n";

const string LongDescription =
  "In this file we explain how to using some methods which are using the Pommaret Basis.  \n";
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
  cout << boolalpha << endl;
  ring Q = RingQQ();
  cout <<"///////////////A homogenous ideal with a Pommaret basis//////////////////////////////////////////" << endl;
  SparsePolyRing polyRing6 = NewPolyRing(Q, SymbolRange("x",0,5));
  vector<RingElem> x = indets(polyRing6);
  vector<RingElem> polys;

  polys.push_back(  power(x[0],2) + power(x[1],2) + power(x[2],2) );
  polys.push_back(  x[0]*x[3] + x[1]*x[4] + x[2]*x[5]);
  polys.push_back(  50*power(x[0],2) - 2*x[0]*x[4] + 14*x[0]*x[5] + power(x[1],2) - 14*x[1]*x[2] + 2*x[1]*x[3] + 49*power(x[2],2) - 14*x[2]*x[3] + power(x[3],2) + power(x[4],2) + power(x[5],2) );
  polys.push_back(  29*power(x[0],2) - 10*x[0]*x[1] - 4*x[0]*x[2] - 4*x[0]*x[4] + 10*x[0]*x[5] + 5*power(x[1],2) - 20*x[1]*x[2] + 4*x[1]*x[3] - 2*x[1]*x[5] + 26*power(x[2],2) - 10*x[2]*x[3] + 2*x[2]*x[4] + power(x[3],2) + power(x[4],2) + power(x[5],2) );
  polys.push_back(  9*power(x[0],2) - 6*x[0]*x[5] + 9*power(x[2],2) + 6*x[2]*x[3] + power(x[3],2) + power(x[4],2) + power(x[5],2) );
  polys.push_back(  9*power(x[0],2) + 6*x[0]*x[5] + 9*power(x[2],2) - 6*x[2]*x[3] + power(x[3],2) + power(x[4],2) + power(x[5],2) );


  cout << "---------computing the Janet basis/Pommaret basis:------------------------" << endl;
  Involutive::JBMill jbMill(Involutive::JBMill::Builder().setInput(polys));
  Involutive::PBMill mill(Involutive::PBMill::Converter().setJBMill(jbMill));

  cout << "Janet basis" << endl;
  vector<RingElem> basis = mill.myReturnPB();
  for(std::vector<RingElem>::const_iterator i = basis.begin(); i != basis.end(); ++i)
  {
    cout << *i << endl;
  }

  cout << "-----computing the depth of the ideal" << endl;
  cout << "depth = " << mill.myDepth() << endl << endl;

  cout << "-----computing the projective dimension" << endl;
  cout << "projDim = " << mill.myProjDim() << endl;
  cout << "-----test if I is Cohen Macaulay" << endl;
  cout << "Ideal is a Cohen-Macaulay ideal = " << mill.IamCohenMacaulay() << endl;


  cout << "-----regular sequence" << endl;
  std::vector<RingElem> regularSeq = mill.myRegSeq();
  for (std::vector<RingElem>::iterator i = regularSeq.begin(); i != regularSeq.end(); ++i)
  {
    if(i != regularSeq.begin())
    {
      cout << ", " << *i;
    }
    else
    {
      cout << *i;
    }
  }
  cout << endl;


  cout << "-----maxStronglyIndependentSet" << endl;
  std::vector<RingElem> stronglySet = mill.myMaxStronglyIndependentSet();
  for (std::vector<RingElem>::iterator i = stronglySet.begin(); i != stronglySet.end(); ++i)
  {
    if(i != stronglySet.begin())
    {
      cout << ", " << *i;
    }
    else
    {
      cout << *i;
    }
  }
  cout << endl;

  cout << "---------regularity" << endl;
  cout << "regularity = " << mill.myRegularity() << endl;

  cout << "----------satiety" << endl;
  cout << "satiety = " << mill.mySatiety() << endl;

  cout << "----------saturation" << endl;
  std::vector<RingElem> satur = mill.mySaturation();
  long counter = 0;
  for (std::vector<RingElem>::iterator i = satur.begin(); i != satur.end(); ++i)
  {
    ++counter;
    cout << *i << endl;
  }

  cout << "----------classes stuff" << endl;

  RingElem clsRingElem = x[0]*x[1];

  cout << "cls of " << LPP(clsRingElem) << " = " << mill.myCls(LPP(clsRingElem)) <<  endl;

  cout << "minimal class of result = " << mill.myMinCls() << endl;

  cout << "all elements with cls 3 in the current basis" << endl;
  std::vector<RingElem> clsElems = mill.myElementsWithClass(3);
  for (std::vector<RingElem>::iterator i = clsElems.begin(); i != clsElems.end(); ++i)
  {
    cout << *i << endl;
  }

  cout << endl << endl;

  cout << "extremalBettiNumbers" << endl;

  map<pair<long, long>, long> bettis = mill.myExtremalBettiNumbers();
  for (map<pair<long, long>, long>::iterator i = bettis.begin(); i != bettis.end(); ++i)
  {
    cout << "(" << i->first.first << ", " << i->first.second << ")" << "=" << i->second << endl;
  }
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
