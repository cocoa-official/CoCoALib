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

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    vector<symbol> names = symbols("x,y,z,u");
    vector<long> ElimIndets; // empty
    SparsePolyRing WA = NewWeylAlgebra(RingQQ(), names, ElimIndets);

    RingElem x(WA, symbol("x"));
    RingElem y(WA, symbol("y"));
    RingElem z(WA, symbol("z"));
    RingElem u(WA, symbol("u"));

    RingElem dx(WA, symbol("dx"));
    RingElem dy(WA, symbol("dy"));
    RingElem dz(WA, symbol("dz"));
    RingElem du(WA, symbol("du"));

    ideal I = ideal(dy*dz -dx*du, x*dx -u*du +1, y*dy +u*du +1, z*dz +u*du +2);
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

// gens(I)=[d[1]*d[2] -d[0]*d[3],  x[0]*d[0] -x[3]*d[3] +1,  x[1]*d[1] +x[3]*d[3] +1,  x[2]*d[2] +x[3]*d[3] +2]
// TidyGens(I)=[d[1]*d[2] -d[0]*d[3],  x[2]*d[2] +x[3]*d[3] +2,  x[1]*d[1] +x[3]*d[3] +1,  x[0]*d[0] -x[3]*d[3] +1,  x[2]*d[0]*d[3] +x[3]*d[1]*d[3] +2*d[1],  x[1]*d[0]*d[3] +x[3]*d[2]*d[3] +d[2],  x[0]*x[3]*d[1]*d[3] +x[2]*x[3]*d[3]^2 +2*x[0]*d[1],  x[0]*x[3]*d[2]*d[3] +x[1]*x[3]*d[3]^2 +x[0]*d[2],  x[1]*x[2]*x[3]*d[3]^2 -x[0]*x[3]^2*d[3]^2 -4*x[0]*x[3]*d[3] -2*x[0]]
