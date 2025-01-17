// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <algorithm>
using std::sort;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing how to compute a quotient basis of a 0-dimensional ideal.\n";

const string LongDescription =
  "The function \"QuotientBasis\" is now included in CoCoALib,          \n"
  "so this example is a lot shorter than it was before version 0.9943.  \n"
  "It returns a vector of PPMonoidElem.\n";
// ----------------------------------------------------------------------

// Includes from the standard C++ library
#include <functional> // using std::bind2nd; using std::ptr_fun;
// #include <iostream> // using std::endl;
#include <iterator> // using std::back_insert_iterator;
#include <list> // using std::list;

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring Fp = NewZZmod(32003);                      // coefficient ring
    ring Fpwyz = NewPolyRing_DMPII(Fp, symbols("w,y,z"), lex);
    RingElem w = RingElem(Fpwyz, symbol("w"));
    RingElem y = RingElem(Fpwyz, symbol("y"));
    RingElem z = RingElem(Fpwyz, symbol("z"));
    ideal J = ideal(power(w,2), power(y,3), power(z,3), w*y*z*z);
    cout << "J  = " << J << endl;
    cout << "GBasis(J)  = " << GBasis(J) << endl;  
    cout << "QuotientBasis(J) = " << QuotientBasis(J) << endl;
    cout << endl;
  
    SparsePolyRing Fpx = NewPolyRing_DMPII(Fp, 8); // Fp[x[0..7]]
    const vector<RingElem>& x = indets(Fpx);
    vector<RingElem> g;
    for (long i=0 ; i < NumIndets(Fpx) ; ++i )
      g.push_back(power(x[i],6));
    ideal I = ideal(g) + ideal(power(x[2],2)*power(x[5],4),
                               power(x[1],3)*power(x[4],4));
    cout << "I  = " << I << endl;  
    vector<PPMonoidElem> QB = QuotientBasis(I);
    double t0 = CpuTime();
    sort(QB.begin(),QB.end());
    double t1 = CpuTime();
    cout << "Time to sort: " << t1 -t0 << endl;
    cout << "len(QB) = " << len(QB) << endl;
    cout << "QB[0] = " << QB[0] << endl;
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
