// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to compute a factor base using an object  \n"
  "of type GCDFreeBasis_RingElem.                                   \n";

const string LongDescription =
  "This example shows how to compute a factor base using an object of\n"
  "type GCDFreeBasis_RingElem.  The generators of the factor base    \n"
  "may be supplied individually, or all together in a vector.        \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));

    // --------------------------------------------
    // Ex.1 Supplying generators in a vector:
    vector<RingElem> v;
    v.push_back(RingElem(P, "x*y^2*z^3"));
    v.push_back(RingElem(P, "x*y^4*z^9"));

    cout << "Generators in a vector = " << v << endl;
    CoprimeFactorBasis_RingElem FacB;
    FacB.myAddInfo(v);
    cout << "First CoprimeFactorBasis = " << FactorBase(FacB) << endl << endl;

    // --------------------------------------------
    // Ex.2 Supplying generators one at a time:
    CoprimeFactorBasis_RingElem FacB2;
    cout << "Generators supplied one at a time..." << endl;
    FacB2.myAddInfo(RingElem(P, "x*y^4*z^9"));
    FacB2.myAddInfo(RingElem(P, "x*y^2*z^3"));
    cout << "Second CoprimeFactorBasis = " << FactorBase(FacB2) << endl;
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
