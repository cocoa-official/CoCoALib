// Copyright (c) 2012  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to create and use a PPMonoidHom.  \n";

const string LongDescription =
  "This program shows how to create and use a PPMonoidHom.  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    PPMonoid PPM1 = NewPPMonoid(symbols("x,y,z"), StdDegRevLex);
    PPMonoid PPM2 = NewPPMonoid(symbols("alpha,beta,gamma,delta"), lex);

    const vector<PPMonoidElem>& x = indets(PPM1);
    const vector<PPMonoidElem>& alpha = indets(PPM2);

    vector<PPMonoidElem> images;
    images.push_back(alpha[0]*alpha[1]);
    images.push_back(alpha[1]*alpha[2]);
    images.push_back(alpha[2]*power(alpha[3],4));

    PPMonoidHom phi = GeneralHom(PPM1, images);
    cout << phi << endl;

    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
        {
          PPMonoidElem t = power(x[0],i) * power(x[1],j) * power(x[2],k);
          cout << t << " |--> " << phi(t) << endl;
        }
  }
  
} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
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

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
