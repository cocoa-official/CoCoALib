// Copyright (c) 2005  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <cstdlib>
// using exit

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how read a polynomial from a string or a file.\n";

const string LongDescription =
  "This program shows the easiest way to read a RingElem (not just polys) \n"
  "from an expression written in a string or in a file.                   \n"
  "It is much more refined than what was available/shown in ex-PolyInput1.\n"
  "NB As always, reading from string is convenient but NOT efficient.     \n";

//----------------------------------------------------------------------

// Includes from the standard C++ library
#include<fstream>  // using std::ifstream; using std::ofstream;

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    SparsePolyRing QQxy = NewPolyRing(RingQQ(), symbols("x,y"));
    SparsePolyRing Fpabc = NewPolyRing(NewZZmod(32003), symbols("a,b,c"));
    SparsePolyRing QQx = NewPolyRing(RingQQ(), SymbolRange("x",0,2));

    const char* FileName = "ex-PolyInput2.in";
    ifstream in(FileName);
    if (!in)
    {
      cerr << "Cannot find input file `" << FileName << "'.  Aborting." << endl;
      abort();
    }

    cout << "------ reading form file -----------------" << endl;
    cout << "-- reading f ..." << endl;
    RingElem f = ReadExprSemicolon(Fpabc, in); // NoPrompt instead of cout
    cout << "-- f is " << f << endl << endl;

    cout << "-- reading g ..." << endl;
    RingElem g = ReadExprSemicolon(QQxy, in); // NoPrompt instead of cout
    cout << "-- g is " << g << endl;
    cout << endl;

    cout << "-- reading h ..." << endl;
    RingElem h = ReadExprSemicolon(QQx, in); // NoPrompt instead of cout
    cout << "-- h is " << h << endl;
    cout << endl;

    cout << "------ reading form string -----------------" << endl;
    string s = "(a^200*b*(-2) + 1) * (- a + b)";
    cout << "-- s = \"" << s << "\"" << endl;
    cout << "-- RingElem(Fpabc, s) gives " << RingElem(Fpabc, s) << endl<< endl;

    cout << "-- reading into a vector ..." << endl;
    s = "";
    cout << "-- s = \"" << s << "\"" << endl;
    cout << "-- RingElems(QQxy, s) gives " << RingElems(QQxy, s) << endl<< endl;
    cout << endl;

    s = "x-y,    (x^2 +1) * (x +y),    (x-1)^3";
    cout << "-- s = \"" << s << "\"" << endl;
    cout << "-- RingElems(QQxy, s) gives " << RingElems(QQxy, s) << endl<< endl;
    cout << endl;

    cout << "-- and also in other rings ..." << endl;
    cout << "-- RingElems(NewQuotientRing(QQxy, \"x^2+1, y-4\"), s) gives\n   "
         << RingElems(NewQuotientRing(QQxy, "x^3-1, y-4"), s) << endl<< endl;
    cout << endl;

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
