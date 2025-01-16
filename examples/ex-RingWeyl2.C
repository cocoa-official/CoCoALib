// Copyright (c) 2003-2006  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "An example about RingWeyl, the interface is not quite settled yet.\n";

const string LongDescription =
  "This shows a computation of a Groebner Basis.\n"
  "All these examples about RingWeyl will probably be merged into one.\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // praticamente identico a 3 e 4: unificare?
  // ideali: da dove vengono?

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    vector<symbol> names = symbols("u,v,x,y");
    vector<long> ElimIndets; // empty
    SparsePolyRing WA = NewWeylAlgebra(RingQQ(), names, ElimIndets);

    RingElem x(WA, symbol("x"));
    RingElem y(WA, symbol("y"));
    RingElem u(WA, symbol("u"));
    RingElem v(WA, symbol("v"));

    RingElem dx(WA, symbol("dx"));
    RingElem dy(WA, symbol("dy"));

    ideal I = ideal(3*x*dx + 2*y*dy +6,  3*y*dx + 2*x*dy);
    cout << "gens(I) = " << gens(I) << endl;
    cout << "TidyGens(I) = " << TidyGens(I) << endl;
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

//output
// gens(I) = [3*x*dx +2*y*dy +6,  3*y*dx +2*x*dy]
// TidyGens(I) = [x*dx +2/3*y*dy +2,  1]

// ANNA: so is it [(1)] ???
