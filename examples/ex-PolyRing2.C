// Copyright (c) 2010  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example showing how to write polynomials using \"monomial\" \n"
  "See also ex-PolyRing1.C                                     \n";

const string LongDescription =
  "Because of C/C++ precedence on operator \"^\" we cannot overload\n"
  "it to define powers, as a consequence writing power-products can\n"
  "be quite tedious or difficult to read.                          \n"
  "Here we show how to write a simple function to create exponent  \n"
  "vectors to be passed to \"monomial\".                           \n"
  "I wish there were a better way to initialise a C++ vector...    \n";
//-----------------------------------------------------------------------

namespace CoCoA
{

  vector<long> exps(long a1, long a2, long a3, long a4)
  {
    vector<long> v;
    v.push_back(a1); v.push_back(a2); v.push_back(a3); v.push_back(a4);
    return v;
  }


  void PolyDemo(SparsePolyRing P)
  {
    cout << "----------------------------------------" << endl;
    cout << "Polynomial ring is " << P << endl;
    if (NumIndets(P) != 4)
    {
      cout << "our function exps works only for 4 indets" << endl;
      return;
    }

    // Use the vector x to access easily the indeterminates:
    const vector<RingElem>& x = indets(P);

    RingElem p(P), f(P), g(P), h(P);

    try
    {
      // Writing as in CoCoA-5 with RingElem
      // This works only if the symbols "x[0]".."x[3]" are in P
      p = RingElem(P, "13*x[0]^4*x[1]^3*x[2] -(1/2)*x[0]^5*x[2]*x[3]^2"
                   "+ 1234567890*x[0]^2*x[1]*x[2]*x[3] -1*x[1]");
      cout << "p = " << p << endl;
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      if (message(err) != "symbol not in ring")
      cout << "Some symbol not in ring: p could not be defined" << endl;
        //        ANNOUNCE(cout, err);
    }


    // Writing the terms of a polynomial using monomial(...)
    f = monomial(P, 13,                             exps(4,3,1,0)) 
      + monomial(P, BigRat(-1,2),                   exps(5,0,1,2))
      + monomial(P, BigIntFromString("1234567890"), exps(2,1,1,1))
      + monomial(P,  -1,                            exps(0,1,0,0));
    cout << "f = " << f << endl;

    // or equivalently writing each term as a product of powers
    g = 13 * power(x[0],4) * power(x[1],3) * x[2]
      + BigRat(-1,2) * power(x[0],5) * x[2] * power(x[3],2)
      + BigIntFromString("1234567890") * power(x[0],2) * x[1] * x[2] * x[3]
      - x[1];
    cout << "g = " << g << endl;

    // You can even mix the two ways
    h = monomial(P, 13,           exps(4,3,1,0)) 
      + monomial(P, BigRat(-1,2), exps(5,0,1,2))
      + BigIntFromString("1234567890") * power(x[0],2) * x[1] * x[2] * x[3]
      - x[1];
    // or use strings if you know the indet "names" (see ex-PolyInput2.C)

    cout << "h = " << h << endl;
    cout << "deg(h) = " << deg(h) << endl;
    cout << "wdeg(h) = " << wdeg(h) << endl;
    cout << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring Z11 = NewZZmod(11);
    PolyDemo(NewPolyRing(Z11, NewSymbols(4)));
    PolyDemo(NewPolyRing(Z11, SymbolRange("x",0,3), StdDegLex));
    PolyDemo(NewPolyRing(RingQQ(), symbols("a,b,c,d"), lex));
    PolyDemo(NewPolyRingWeights(RingQQ(), symbols("x[0],x[1],x[2],x[3]"),
                                RowMat(RingElems(RingZZ(),"1,1,2,2"))));
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
