// Copyright (c) 2017  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "The example in this file shows how to create a RingHom from \n"
  "one quotient ring to another.                               \n";

const string LongDescription =
  "The example in this file shows how to create a RingHom from \n"
  "one quotient ring to another.  First we make a RingHom from \n"
  "a polynomial ring to the codomain (a quotient ring), then we\n"
  "use the function `InducedHom' to extend it to the full domain.\n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // We shall build a RingHom from the quotient ring
    // Q1 := QQ[x,y,z]/ideal(x^2+y^2+z^2-1)
    // to the quotient ring
    // Q2 := QQ[a,b]/ideal(a^2+4*b^2-1)
    
    const ring P1 = NewPolyRing(RingQQ(), symbols("x,y,z"));
    const ideal I1 = ideal(RingElem(P1,"x^2+y^2+z^2-1"));
    const QuotientRing Q1 = NewQuotientRing(P1,I1);

    const ring P2 = NewPolyRing(RingQQ(), symbols("a,b"));
    const ideal I2 = ideal(RingElem(P2, "a^2+4*b^2-1"));
    const QuotientRing Q2 = NewQuotientRing(P2,I2);

    // First build PolyAlgebraHom from P1 to Q2
    vector<RingElem> images;
    images.push_back(RingElem(Q2, "2*b"));
    images.push_back(RingElem(Q2, "0"));
    images.push_back(RingElem(Q2, "a"));
    const RingHom phi = PolyAlgebraHom(P1, Q2, images);

    // Now use InducedHom to extend to hom from Q1 to Q2:
    const RingHom psi = InducedHom(Q1, phi);

    // Print out the RingHom -- printed form is ugly/hard-to-read :-(
    cout << "psi is " << psi << endl << endl;
    
    cout << "psi(z^2) = " << psi(RingElem(Q1, "z^2")) << endl;
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

