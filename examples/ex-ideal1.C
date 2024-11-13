// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to create ideals from ring elements. \n";

const string LongDescription =
  "This short example shows how to create ideals from one or more \n"
  "ring elements.  If there are at most 4 generators this can be  \n"
  "done directly; otherwise put the generators in a vector, and   \n"
  "create the ideal from the vector.                              \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // First we create a PolyRing in which to compute -- see also ex-PolyRing1
    PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    // Two polynomials...
    RingElem f = RingElem(P, "x*y-1");
    RingElem g = RingElem(P, "x^2+y^2-1");

    // Now create some ideals
    ideal I1 = ideal(f);
    ideal I2 = ideal(f,g);
    cout << "Ideal I1 = " << I1 << "  and I2 = " << I2 << endl;

    // To create an ideal with many generators, put them in a vector:
    vector<RingElem> GenList;
    GenList.push_back(f);
    GenList.push_back(g);
    GenList.push_back(f*g);
    GenList.push_back(indet(P,2));
    ideal I3 = ideal(GenList);

    // Alternative form (which works even if GenList is empty)
    ideal I4 = ideal(P, GenList);
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
