// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <cstdlib>
// using exit

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows a very crude way to read a polynomial from the input;\n"
  "this could be useful for importing some kind of raw data.\n"
  "See ex-PolyInput2.C for much better, human usable, functions.\n";
//----------------------------------------------------------------------

// Includes from the standard C++ library
// #include<iostream> // using std::endl; using std::flush;
#include<fstream>  // using std::ifstream; using std::ofstream;

namespace CoCoA
{

  RingElem ReadNewPoly(ostream& prompt, istream& in, const SparsePolyRing& P)
  {
    RingElem f(P);
    size_t InputNumTerms, NumInds = NumIndets(P);
    BigInt IntCoeff;
    vector<long> expv(NumInds);

    // look in program() for changing Global*put from default values
    prompt << "how many terms? " << flush;
    in     >> InputNumTerms;
    if ( !in )
    {
      cerr << "*** ERROR *** Input must be a positive integer" << endl;
      exit(1);
    }
    for (size_t t=1 ; t <= InputNumTerms ; ++t) 
    {
      prompt << "[" << t << "]\tcoeff (integer)? " << flush;
      in     >> IntCoeff;
      prompt << "   \t" << NumInds << " exponents? (separated by spaces) " << flush;
      for (size_t i=0; i<NumInds; ++i)  in >> expv[i];
      f += monomial(P, IntCoeff, expv);
    }
    return f;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    SparsePolyRing Qxy = NewPolyRing(RingQQ(), symbols("x,y"));
    SparsePolyRing Fpabc = NewPolyRing(NewZZmod(32003), symbols("a,b,c"));
    SparsePolyRing Pxy = NewPolyRing(Fpabc, symbols("x,y"));

    cout << "We start by reading a polynomial from stdin: " << endl;
    cout << "-- reading h ..." << endl;
    RingElem h = ReadNewPoly(cout, cin, Qxy);
    cout << "-- h in " << owner(h) << endl;
    cout << "-- h     = " << h << endl;

    cout << endl;
    cout << "... and now we read f and g from file \"ex-PolyInput1.in\": " << endl;
    //----------------------------------------------------------------------
    // we can, at any time, change the default streams in, out, log, err like this:
    const char* FileName = "ex-PolyInput1.in";
    ifstream in(FileName);
    if (!in)
    {
      cerr << "Cannot find input file `" << FileName << "'.  Aborting." << endl;
      abort();
    }

    // now we read from file: so we do not want "ReadNewPoly" to print:
    ofstream NoPrompt("/dev/null");  // output-file-stream /dev/null

    cout << "-- reading f ..." << endl;
    RingElem f = ReadNewPoly(NoPrompt, in, Fpabc); // NoPrompt instead of cout
    cout << "-- reading g ..." << endl;
    RingElem g = ReadNewPoly(NoPrompt, in, Pxy); // NoPrompt instead of cout
    cout << endl;

    cout << "-- f in " << owner(f) << endl;
    cout << "-- f     = " << f << endl;
    cout << endl;

    cout << "The printed form of g will have parentheses around each" << endl;
    cout << "coefficient to remind the user that coefficients" << endl;
    cout << "are not \"simple numbers\"" << endl;
    cout << "-- g in " << owner(g) << endl;
    cout << "-- g     = " << g << endl;
    cout << endl;

    RingHom phi = CoeffEmbeddingHom(Pxy);
    cout << "-- phi(f) * g = " << phi(f) * g << endl;
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
