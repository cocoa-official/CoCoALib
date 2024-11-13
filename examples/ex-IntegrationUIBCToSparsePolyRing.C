// Copyright (c) 2015  John Abbott,  Anna M Bigatti
// Orig author: 2015  Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a template file for example programs.  \n"
  "The program itself does nothing whatsoever.    \n";

const string LongDescription =
  "Make a copy of this file (called \"foo.C\", say) and put your code \n"
  "inside the procedure \"program\".                                  \n"
  "To compile your file in the examples directory just do this:       \n"
  "  make foo                                                         \n";
//----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  cout << boolalpha; // so that bools print out as true/false
  ring Q = RingQQ();

  SparsePolyRing P = NewPolyRing(Q, symbols("x[0],x[1]"));
  const vector<RingElem>& x = indets(P);
  ideal I = ideal(x[0]*x[0], x[0] * x[1]);
  std::cout <<  "JanetBasis(I)    = " << Involutive::JanetBasis(I) << std::endl;
  std::cout << "IsDeltaRegular(I) = " << Involutive::IsDeltaRegular(I) << std::endl;
  std::cout << "IsMonomial(I)     = " << Involutive::IsMonomial(I) << std::endl;
  std::cout << "IsHomogeneous(I)  = " << Involutive::IsHomogeneous(I) << std::endl;
  std::map<PPMonoidElem, std::vector<bool> > multVars(Involutive::MultVars(I));

  std::cout << "MultVars" << std::endl;
  for (std::map<PPMonoidElem, std::vector<bool> >::iterator i = multVars.begin(); i != multVars.end(); ++i)
  {
    std::cout << i->first << ": " << i->second<< std::endl;
  }
  std::cout << std::endl;
  std::map<PPMonoidElem, std::vector<bool> > nonMultVars(Involutive::NonMultVars(I));

  std::cout << "NonMultVars" << std::endl;
  for (std::map<PPMonoidElem, std::vector<bool> >::iterator i = nonMultVars.begin(); i != nonMultVars.end(); ++i)
  {
    std::cout << i->first << ": " << i->second<< std::endl;
  }
  std::cout << std::endl;
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
