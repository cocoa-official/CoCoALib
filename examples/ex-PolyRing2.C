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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PolyRing2.C,v 1.21 2024/03/18 16:06:05 bigatti Exp $
// $Log: ex-PolyRing2.C,v $
// Revision 1.21  2024/03/18 16:06:05  bigatti
// Summary: improved readability of output;  catching (wanted) error;
// added call to ring with weights and deg/wdeg
//
// Revision 1.20  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.19  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.18  2018/04/20 18:51:24  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.17  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.16  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.15  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.14  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.13  2015/06/26 14:54:55  abbott
// Summary: Revised after update to exceptions/errors
// Author: JAA
//
// Revision 1.12  2015/04/27 12:51:08  bigatti
// Summary: added example of ReadExpr
//
// Revision 1.11  2014/03/21 16:44:53  bigatti
// -- added comment about ReadExpr
//
// Revision 1.10  2013/05/27 17:13:50  abbott
// Minor improvement to long description.
//
// Revision 1.9  2012/02/08 17:42:26  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.8  2011/09/06 13:36:29  abbott
// Layout change.
//
// Revision 1.7  2011/08/24 10:46:32  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.6  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2010/03/23 08:32:13  bigatti
// -- improved: 3 poly rings, more interesting coefficients
//
// Revision 1.3  2010/03/22 10:16:40  abbott
// Removed a useless pair of brackets.
//
// Revision 1.2  2010/03/18 18:01:50  bigatti
// -- added mixed sum
//
// Revision 1.1  2010/03/18 17:32:55  bigatti
// -- first import
//
// Revision 1.10  2010/03/15 12:03:04  bigatti
// -- rearranged margins for descriptions
// -- added one constant to g
//
// Revision 1.9  2010/03/11 11:01:35  abbott
// Removed pointless print command.
// Added some helpful(?) comments.
// Added a new example.
//
// Revision 1.8  2010/03/08 08:27:03  bigatti
// -- added constant to polynomial g
//
// Revision 1.7  2010/03/05 21:33:55  abbott
// Added long description.
// Changed an example to use the new PPOrderingCtor arg, and a few
// other changes to improve readability.
//
// Revision 1.6  2009/06/05 12:14:56  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from BigInt.H/C
//   removed gcd from MachineInteger.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
// Revision 1.5  2008/07/21 07:58:51  abbott
// Several cosmetic changes to examples.
//
// Revision 1.4  2008/03/12 16:35:19  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.3  2007/12/05 12:27:54  bigatti
// -- minor fix
//
// Revision 1.2  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.11  2007/03/08 22:26:27  cocoa
// Removed try..catch constructs from some "simple" examples.
//
// Revision 1.10  2007/03/08 17:26:32  bigatti
// -- improved NewPolyRing
//
// Revision 1.9  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.8  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.7  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/02/28 14:51:19  bigatti
// -- simplified call for Z/(BigPrime)
//
// Revision 1.5  2007/02/26 15:50:25  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/12 15:45:15  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.3  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2006/07/19 07:10:27  cocoa
// -- added Q[x,y,z]
// -- use of "sum" and other minor changes
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/03/01 14:26:03  cocoa
// -- added Z / (big prime)
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.5  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.4  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.3  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.2  2005/06/30 15:53:37  cocoa
// -- some comments were wrong (indexed indets)
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.1  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from BigInt).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
