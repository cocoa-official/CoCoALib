// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how we define a ring homomorphism       \n"
  "to perform a change of coordinates in a polynomial ring.   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create coeff ring, then poly ring.
    ring Fp = NewZZmod(32003);
    PolyRing Fpx = NewPolyRing(Fp, symbols("a,b,c,d")); // Fpx is Z/(32003)[a,b,c,d]

    cout << "We shall work in the PolyRing " << Fpx << endl;

    // Line below gives easy access to the indets: they are just x[0], x[1],...
    const vector<RingElem>& x = indets(Fpx);

    // Careful with rationals in C++ because C++ treats 2/3 as an integer division
    const RingElem TwoThirds = RingElem(Fpx,2)/3;

    RingElem f = 2*x[0] + 4*x[1]*x[1] - TwoThirds*x[0]*x[2];
    cout << "Original poly is f = " << f << endl;

    // Now create the poly algebra homomorphism; start with vector of images:
    vector<RingElem> images;
    images.push_back(zero(Fpx));       // x[0] |--> 0
    images.push_back(x[0] + 2*x[1]);   // x[1] |--> x[0] + 2*x[1]
    images.push_back(one(Fpx));        // x[2] |--> 1
    images.push_back(TwoThirds*x[3]);  // x[3] |--> (2/3)*x[3]

    RingHom phi = PolyAlgebraHom(Fpx, Fpx, images);

    cout << "We have built the following homomorphism"
         << " which we have called phi:" << endl
         << phi << endl
         << "Here are some values calculated using phi:" << endl
         << endl;

    cout << x[0] << " |--> " << phi(x[0]) << endl;
    cout << x[1] << " |--> " << phi(x[1]) << endl;
    cout << x[2] << " |--> " << phi(x[2]) << endl;
    cout << x[3] << " |--> " << phi(x[3]) << endl;
    cout << "f    |--> " << phi(f) << endl;
    cout << endl;

    ideal I(x[0]+x[1], x[2]-x[3]);

    const vector<RingElem>& h=gens(I);
    long NGens = len(h);
    vector<RingElem> v;
    v.reserve(NGens); // not necessary but good style: reserves memory (C++ STL)
    for (long i=0; i < NGens; ++i)  v.push_back(phi(h[i])); 
    cout << "v = " << v << endl;
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
    cerr << "***ERROR***  UNCAUGHT CoCoAError";
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
