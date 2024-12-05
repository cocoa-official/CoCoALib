// Copyright (c) 2013  John Abbott, Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows how to extract the \"coefficients\" \n"
  "of a polynomial, and also how to compute some other operations \n"
  "on the coefficients                                            \n";

const string LongDescription =
  "This example program shows how to extract the \"coefficients\" \n"
  "of a polynomial, and also how to compute some other operations \n"
  "on the coefficients (e.g. content).  CoCoALib offers two notions\n"
  "of coefficient: one is the natural one dictated by the ring,   \n"
  "and the other comes from identifying e.g. k[x,y,z] and k[x,z][y]\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    SparsePolyRing ZZxy = NewPolyRing(RingZZ(), symbols("x,y"));
    const RingElem x = indet(ZZxy,0);
    const RingElem y = indet(ZZxy,1);

    RingElem f = 2*x*x*y - 4*y*y + 6*x*x + 36;

    cout << "In the following we have f = " << f << endl << endl;

    cout << "Product of indets actually in f is " << IndetsProd(f) << endl << endl;

    // Accessing coeffs via SparsePolyIter:
    cout << "Using a SparsePolyIter we decompose f as follows:" << endl;
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      cout << "coeff = " << coeff(it) << "  in ring " << owner(coeff(it)) << endl;
      cout << "PP    = " << PP(it)    << "  in " << owner(PP(it)) << endl;
      cout << endl;
    }

    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << endl;

    // Regard f as a poly in just "x" or just "y" we obtain:
    cout << "Coefficients with respect to certain indeterminates:" << endl << endl;
    cout << "Case (1) considering f as a *univariate* polynomial in..." << endl;
    cout << "...the indeterminate x the coeffs are " << CoeffVecWRT(f, x) << endl;
    cout << "...the indeterminate y the coeffs are " << CoeffVecWRT(f, y) << endl;
    cout << endl;

    cout << "Case (2) considering f as a sparse multivariate polynomial in..." << endl;
    cout << "...the indet x its structure is " << CoefficientsWRT(f, x) << endl;
    cout << "...the indet y its structure is " << CoefficientsWRT(f, y) << endl;
    vector<long> XandY; XandY.push_back(0); XandY.push_back(1);
    cout << "...the indets x & y its structure is " << CoefficientsWRT(f, XandY) << endl;
    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "The content of a polynomial" << endl << endl;

    // Content of f
    RingElem ContF = content(f);
    cout << "The \"numerical\" content of f is " << ContF << "  -- an element of ring " << owner(ContF) << endl;
    cout << endl;
    RingElem ContWRTx = ContentWRT(f, x);
    cout << "Content WRT x is " << ContWRTx << "  -- element of " << owner(ContWRTx) << endl;

    RingElem ContWRTy = ContentWRT(f, y);
    cout << "Content WRT y is " << ContWRTy << "  -- element of " << owner(ContWRTy) << endl;

    RingElem ContWRTxy = ContentWRT(f, XandY);
    cout << "Content WRT x & y is " << ContWRTxy << "  -- element of " << owner(ContWRTxy) << endl;

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
