// Copyright (c) 2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows simple use of the functions for creating new   \n"
  "anonymous symbols.                                                \n";

const string LongDescription =
  "Often algorithms need one to work with an addition new indeterminate\n"
  "which must be different from all indeterminates already being use.  \n"
  "In CoCoALib we can create new \"anonymous\" indeterminates guaranteed\n"
  "to be different from all others.  This examples shows how to create \n"
  "them singly, and also several all together.                         \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Function to count how many distinct irred factors x^n-1 has in k[x].
  // This function will work regardless of which symbols appear in k.
  long NumFactors(ring k, int n)
  {
    PolyRing P = NewPolyRing(k, NewSymbols(1));
    const RingElem x = indet(P,0);
    const factorization<RingElem> facs = factor(power(x,n)-1);
    return len(facs.myFactors());
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // To create a single new symbol use the function NewSymbol().
    // Each time you call it you will get a different symbol.
    cout << "First new symbol: " << NewSymbol() << endl;
    cout << "Second new symbol: " << NewSymbol() << endl;

    // If you want several new symbols, use the function NewSymbols(n).
    vector<symbol> v = NewSymbols(5);
    cout << "Here are " << len(v) << " newly minted symbols: " << v << endl;

    ring k1 = RingQQ();
    long n1 = NumFactors(k1, 8);
    cout << "x^8-1  has  " << n1 << " factors over " << k1 << endl;
    ring k2 = NewZZmod(257);
    long n2 = NumFactors(k2, 8);
    cout << "x^8-1  has  " << n2 << " factors over " << k2 << endl;

    PolyRing Qx = NewPolyRing(k1, symbols("x"));
    ring k3 = NewQuotientRing(Qx, "x^2+1"); // QQ(i)
    // same as  NewQuotientRing(Qx, ideal(RingElem(Qx, "x^2+1")));
    long n3 = NumFactors(k1, 8);
    cout << "x^8-1  has  " << n3 << " factors over " << k3 << endl;
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
