// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to use geobuckets for long sums.          \n"
  "It compares timings between normal sum \"+=\" and geobucket sum. \n";

const string LongDescription =
  "We simulate a long sum:\n"
  "we add the summands of a long polynomial f (and f*x) one at a time,   \n"
  "and do it twice to consider also the arithmetics on the coefficients. \n"
  "The advantage for geobuckets comes when f is quite long.              \n"
  "Geobuckets are to be used as a temporary value;                       \n"
  "the final result is then copied into a RingElem.";

//----------------------------------------------------------------------

namespace CoCoA
{

  void CompareSums(SparsePolyRing P)
  {
    RingElem f = power(sum(indets(P)) + 1, 9);
    RingElem x = indet(P,0);
    RingElem h_poly(P);
    RingElem h_gbk(P); // just the copy of gbk final value
    geobucket gbk(P);
    cout << "----------------------------------------------------" << endl;
    cout << "---- P is          " << P << endl;
    cout << "---- NumTerms(f) = " << NumTerms(f) << endl;
    cout << "---- sum of monomials in f and x*f:" << endl;

    double t0 = CpuTime();
    for (long i=0; i<2; ++i)
    {
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
        h_poly += monomial(P, coeff(it), PP(it))*x;
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
        h_poly += monomial(P, coeff(it), PP(it));
    }
    cout << "sum using \"+=\"      time = " << CpuTime()-t0 << endl;  

    t0 = CpuTime();
    for (long i=0; i<2; ++i)
    {
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
        gbk.myAddMulLM(monomial(P, coeff(it), PP(it))*x, one(P), 1);
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
        gbk.myAddMulLM(monomial(P, coeff(it), PP(it)), one(P), 1);
    }
    AddClear(h_gbk, gbk); // copies gbk to f_gbk and clears gbk
    cout << "sum using geobucket time = " << CpuTime()-t0 << endl;

    cout <<  "check:  2*(1+x)*f == h_poly is " << (2*(1+x)*f == h_poly);
    cout <<  ";  2*(1+x)*f == h_gbk is " << (2*(1+x)*f == h_gbk) << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << LongDescription << endl;
    cout << boolalpha; // prints true/false for bool

    CompareSums(NewPolyRing(RingQQ(), symbols("x,y,z")));
    CompareSums(NewPolyRing(RingQQ(), symbols("a,b,c,d")));
    CompareSums(NewPolyRing(RingQQ(), SymbolRange("x",1,6)));
    CompareSums(NewPolyRing_DMPI(RingQQ(), SymbolRange("y",1,6)));
    CompareSums(NewPolyRing_DMPII(NewZZmod(101), SymbolRange("z",1,6)));
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
