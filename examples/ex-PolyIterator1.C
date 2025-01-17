// Copyright (c) 2005 John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing how to iterate through \"sparse\" polynomials.\n";

const string LongDescription =
  "Unlike in CoCoA-5, there no functions in CoCoALib to extract the      \n"
  "support, the coeff list or the monomials from a polynomial.  Instead  \n"
  "CoCoALib offers an \"iterator\" over the terms in a sparse polynomial:\n"
  "you can access the coeff and PP directly from the iterator.           \n";

// ----------------------------------------------------------------------

namespace CoCoA
{

  void SimpleDemo(const SparsePolyRing& P)
  {
    if (NumIndets(P) < 4) return;       // Do nothing if there are fewer than 4 indets...
    if (!IsInvertible(RingElem(P,3))) return; // ... or if 3 is not invertible.

    const vector<RingElem>& x = indets(P);
    RingElem f = 2*x[0] + 4*x[1]*x[1] - 2*x[0]*x[3]/3;
    cout << "Our poly f = " << f << "   element of " << owner(f) << endl;

    // ---- ITERATION over f: we just print out the coeffs and the PPs ----
    cout << "The terms in f are as follows\n";
    for (SparsePolyIter iter=BeginIter(f); !IsEnded(iter); ++iter)
    {
      cout << "coeff: " << coeff(iter) << "\t   element of " << owner(coeff(iter)) << endl
           << "PP: " << PP(iter) << "\t   element of " << owner(PP(iter)) << endl
           << endl;
    }
    cout << "Polynomial iterators are read-only, so f is unchanged: " << f << endl
         << "---------------------------------" << endl << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring Fp = NewZZmod(32003);          // Coefficient ring.
    // Now create 3 poly rings over Fp; they are all isomorphic to ZZ/(32003)[a,b,c],
    // but specify different internal implementations (and use different indet names).
    SparsePolyRing Fpx = NewPolyRing(Fp, SymbolRange("x",0,3));
    SparsePolyRing Fpy = NewPolyRing_DMPI(Fp, SymbolRange("y",0,3));
    SparsePolyRing Fpz = NewPolyRing_DMPII(Fp, SymbolRange("z",0,3));

    // Now run the demo in each of these rings:
    SimpleDemo(Fpx);
    SimpleDemo(Fpy);
    SimpleDemo(Fpz);
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
