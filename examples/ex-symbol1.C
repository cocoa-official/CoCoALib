// Copyright (c) 2008  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Creation of symbols, and some simple operations on them.\n";

const string LongDescription =
  "Symbols are used to give print names to indeterminates.  Their main   \n"
  "use is as an argument to a pseudo-ctor for a PPMonoid or a polynomial \n"
  "ring -- see examples for those types.                                 \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Symbols without subscripts" << endl;
    symbol x("x");
    symbol alpha("alpha");
    cout << x << endl;
    cout << alpha << endl << endl;  

    cout << "Symbols with 1 subscript" << endl;
    symbol y2("y",2);
    symbol zeta("zeta",0);
    cout << y2 << endl;
    cout << zeta << endl << endl;  

    cout << "Symbols with 2 subscripts" << endl;
    symbol amn("a",4,7);
    symbol s("Y_phi", -3,8);
    cout << amn << endl;
    cout << s << endl << endl;

    cout <<"Symbols with many subscripts -- flexible but cumbersome!" << endl;
    vector<long> subs = {1,3,5};
    symbol XX("X", subs);
    cout << XX << endl << endl;

    vector<symbol> v;
    v.push_back(alpha);
    v.push_back(XX);
    v.push_back(s);
    cout << "Constructors of polynomial ring requires a vector of symbols\n"
         << "to give PRINT NAMES to indeterminates, for example\n"
         << v << endl << endl;

    cout << "There are convenience functions for making vectors of symbols"
         << endl;
    cout << "- vector of length 1:  " << symbols("omega") << endl;
    cout << "- vector of length 2:  " << symbols("mu, nu") << endl;
    cout << "- vector with range of subscripts:  " 
         << SymbolRange("B", -2,2) << endl;
    cout << "- vector with a \"rectangle\" of double subscripts:\n  "
         << SymbolRange(symbol("y",0,1), symbol("y",2,2)) << endl << endl;

    cout << "Querying a symbol to get its head, number and values of subscripts "
         << endl;
    cout << "s\t is " << s
         << "\t head(s)\t is  " << head(s) << endl;
    cout << "XX\t is "<< XX
         << "\t NumSubscripts(XX)\t is  " << NumSubscripts(XX) << endl;
    cout << "y2\t is " << y2
         << "\t subscript(y2,0)\t is  " << subscript(y2, 0) << endl;

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
