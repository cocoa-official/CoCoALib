// Copyright (c) 2007  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "The example in this file shows how to create and use some \n"
  "homomorphisms between rings.                              \n";

const string LongDescription =
  "CanonicalHom is an easy way to make these homomorphisms:  \n"
  "R --> R/I,    R --> R[x],   R --> FractionField(R),       \n"
  "R --> R,      QQ --> R,     ZZ --> R,                     \n"
  "PolyAlgebraHom makes the R-algebra homomorphisms:         \n"
  "R[x] --> R,   R[x] --> R[y]                               \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create some rings
    ring QQ = RingQQ();
    ring Fp = NewZZmod(101);
    PolyRing P = NewPolyRing(Fp, symbols("a,b"));   // P is Fp[a,b]

    cout << "-- CanonicalHom into P = " << P << endl;

    RingHom FpToP = CanonicalHom(Fp, P);  // same as CoeffEmbeddingHom(P)
    RingHom QQToP = CanonicalHom(QQ, P);  // same as QQEmbeddingHom(P)
    // NB!! QQToP is a partial homomorphism:
    // e.g. we cannot compute QQToP(one(QQ)/101)

    RingElem a = RingElem(Fp, 13);
    cout << "FpToP(" << a << ") = " << FpToP(a) << endl;

    RingElem q = RingElem(QQ, 5)/2;
    cout << "QQToP(" << q << ") = " << QQToP(q) << endl;
    cout << "  same as (RingElem calls CanonicalHom)" << endl;
    cout << "  RingElem(P,q) = " << RingElem(P,q) << endl;
    cout << endl;
 
    vector<RingElem> IndetImages;
    IndetImages.push_back(RingElem(Fp,2));
    IndetImages.push_back(RingElem(Fp,5));
    cout << "-- PolyAlgebraHom:    "
         << indet(P,0) << " |--> " << IndetImages[0] << "   "
         << indet(P,1) << " |--> " << IndetImages[1] << endl;

    RingHom PToFp = PolyAlgebraHom(P, Fp, IndetImages);
  
    RingElem f = 100*indet(P,0) + indet(P,1);
    cout << "PToFp(" << f << ") = " << PToFp(f) << endl;
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
