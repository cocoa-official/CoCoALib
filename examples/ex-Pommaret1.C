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
  "In this file we explain how to compute several stability positions.  \n";

const string LongDescription =
  "In this file we explain how to compute several stability positions.  \n";
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
  // This example shows does not use coordinate transformations
  // The example shows how to convert a Janet basis which is delta regular to a Pommaret basis
  cout << "First Example: Convert a Janet basis to a Pommaret basis" << endl;
  SparsePolyRing polyRing1 = NewPolyRing(Q, symbols("x,y"));
  vector<RingElem> x = indets(polyRing1);
  vector<RingElem> gens;

  gens.push_back(x[0] * x[0]);
  gens.push_back(x[0] * x[1]);

  // introduce JBMill
  Involutive::JBMill jbMill(Involutive::JBMill::Builder().setInput(gens));
  // convert to PBMill
  // if jbMill.IamPommaretBasis() is false the method would raise an error
  Involutive::PBMill mill1(Involutive::PBMill::Converter().setJBMill(jbMill));
  cout << "Pommaret Basis" << std::endl;
  output(mill1.myReturnPB());



  cout << "Second Example: Convert a non-delta-regular input set to delta regular coordinates and compute the Pommaret basis" << endl;
  SparsePolyRing polyRing2 = NewPolyRing(Q, symbols("x,y"));
  x = indets(polyRing2);
  gens.clear();

  gens.push_back(x[1] * x[1]);
  gens.push_back(x[0] * x[1]);

  Involutive::PBMill::DeltaRegularTransformator builder1;
  builder1.setInput(gens);
  // default: SingleWithoutPermutation
  // or AllWithPermutation or AllWithoutPermutation
  builder1.setStrategy(Involutive::PBMill::DeltaRegularTransformator::SingleWithPermutation);
  // during the computation we compute several Janet basis. This method defines the algorithm which
  // should be used during this computation
  builder1.setJanetStrategy(Involutive::GBCompletion);
  // finally building the Pommaret basis
  std::cout << __FUNCTION__ << "1" << std::endl;
  Involutive::PBMill mill2(builder1);
  cout << "Pommaret Basis" << std::endl;
  output(mill2.myReturnPB());



  cout << "Third Example: Convert an input set to stable position and compute the Pommaret basis" << endl;
  // SparsePolyRing polyRing2 = NewPolyRing(Q, symbols("x,y"));
  // x = indets(polyRing2);
  gens.clear();

  gens.push_back(x[0] * x[0] * x[0]);
  gens.push_back(x[1]);

  Involutive::PBMill::StableLTITransformator builder2;
  builder2.setInput(gens);
  // default: Single
  builder2.setStrategy(Involutive::PBMill::StableLTITransformator::All);
  // during the computation we compute several Janet basis. This method defines the algorithm which
  // should be used during this computation
  builder2.setJanetStrategy(Involutive::GBCompletion);
  // finally building the Pommaret basis
  Involutive::PBMill mill3(builder2);
  cout << "Pommaret Basis" << std::endl;
  output(mill3.myReturnPB());



  cout << "Fourth Example: Convert an input set to strongly stable position and compute the Pommaret basis" << endl;
  SparsePolyRing polyRing4 = NewPolyRing(Q, symbols("x,y,z"));
  x = indets(polyRing4);
  gens.clear();

  gens.push_back(x[0] * x[1] * x[2]);
  gens.push_back(x[0] * x[1] * x[1]);
  gens.push_back(x[0] * x[0] * x[1]);
  gens.push_back(x[0] * x[0] * x[0]);

  Involutive::PBMill::StronglyStableLTITransformator builder3;
  builder3.setInput(gens);
  // default: Single
  builder3.setStrategy(Involutive::PBMill::StronglyStableLTITransformator::All);
  // during the computation we compute several Janet basis. This method defines the algorithm which
  // should be used during this computation
  builder3.setJanetStrategy(Involutive::GBCompletion);
  // finally building the Pommaret basis
  Involutive::PBMill mill4(builder3);
  cout << "Pommaret Basis" << std::endl;
  output(mill4.myReturnPB());
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
