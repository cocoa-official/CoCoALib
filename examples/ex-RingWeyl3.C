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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingWeyl3.C,v 1.9 2015/07/27 11:50:50 bigatti Exp $
// $Log: ex-RingWeyl3.C,v $
// Revision 1.9  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.8  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.7  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.6  2012/02/08 17:48:04  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.5  2011/03/08 18:03:52  bigatti
// -- changed size_t into long
// -- using RingElem(ring, symbol) ctor instead of "derivation"
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/05/14 09:45:29  bigatti
// -- improved syntax/style
//
// Revision 1.2  2007/09/24 14:12:37  abbott
// Added missing newline in a string.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 17:18:22  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 16:15:37  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1  2006/08/17 10:08:18  cocoa
// -- modernized version of WeylAlgebraN example
//
