// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how we define a ring homomorphism       \n"
  "to evaluate polynomials.                                   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  RingElem EvalAtZero(RingElem f) //--> the constant term
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    for (SparsePolyIter it= BeginIter(f); !IsEnded(it); ++it)
      if (IsOne(PP(it))) return coeff(it);
    return zero(CoeffRing(owner(f)));
  }
  

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create some coefficient rings
    ring Fp = NewZZmod(32003);

    PolyRing Fpx = NewPolyRing(Fp, SymbolRange("x", 0,3)); // Fpx is Z/(32003)[x[0..3]]
    cout << "  --- PolyRing is " << Fpx << endl;

    const vector<RingElem>& x = indets(Fpx);
    RingElem f = 4*x[1]*x[1] + x[0]*x[2] - 2*x[0] +5;
    cout << "  --- f = " << f << endl;

    cout << "  --- EvalAtZero(f) = " << EvalAtZero(f) << endl;

    cout << "  --- Definition of PolyAlgebraHom  eval" << endl;
    RingHom eval = PolyAlgebraHom(Fpx, Fp, "1, 0, 2/3, 100");

    for (long i=0; i<4; ++i)
      cout << x[i] << " |-> " << eval(x[i]) << endl;
    cout << "f    |-> " << eval(f) << endl;

    cout << " --- Explicit construction of a C++ vector  images" << endl;
    vector<RingElem> images;
    images.push_back(one(Fp));
    images.push_back(zero(Fp));
    images.push_back(RingElem(Fp,"2/3")); // x[2] |-> 2/3
    // using string syntax, because C++ reads 2/3 as an integer division!
    images.push_back(RingElem(Fp,100)); // x[3] |-> 100

    eval = PolyAlgebraHom(Fpx, Fp, images);
    // same as PolyRingHom(Fpx, Fp, IdentityHom(Fp), images);
    cout << "f    |-> " << eval(f) << endl;

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
