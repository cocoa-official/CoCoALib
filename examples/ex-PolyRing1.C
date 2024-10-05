// Copyright (c) 2005-2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example showing how to create some simple polynomial rings.             \n"
  "It also shows some of the operations specific to elements of PolyRings. \n"
  "See also ex-ring2.C                                                     \n";

const string LongDescription =
  "This example program exhibits two things: various ways of creating \n"
  "polynomial rings, and several operations specific to elements of a \n"
  "polynomial ring.                                                   \n"
  "In the procedure \"program\" there are examples of creating        \n"
  "polynomial rings.                                                  \n"
  "In the procedure \"SomeComputations\" there are brief examples of  \n"
  "the operations specific to elements of PolyRings (e.g. deg).       \n"
  "See also ex-ring2.C.                                               \n";
//-----------------------------------------------------------------------

namespace CoCoA
{

  void SomeComputations(const PolyRing& P)
  {
    const long n = NumIndets(P);
    const vector<RingElem>& x = indets(P);

    // Put f = x[0] + x[1] + .. x[n-1]
    RingElem f(sum(indets(P)));

    // Put g = x[0] + 2*x[1] + ... + n*x[n-1] + 1234
    RingElem g(P, 1234);
    for (long i=0; i < n; ++i)  g += (i+1)*x[i];

    cout << "Some computations in the ring " << P << endl;
    cout << "Let f = " << f << endl
         << "and g = " << g << endl << endl;

    cout << "NumTerms(f) = " << NumTerms(f)
         << "   and NumTerms(g) = " << NumTerms(g) << endl;
    cout << "f+g = " << f+g << endl;
    cout << "f-g = " << f-g << endl;
    cout << "f*g = " << f*g << endl;
    if (IsDivisible(f, g))  cout << "f/g = " << f/g << endl;

    cout << "deg(f) = "  << deg(f) << endl;      // of type long
    cout << "StdDeg(f) = " << StdDeg(f) << endl; // same as deg(f)
    cout << "wdeg(f) = " << wdeg(f) << endl;     // of type CoCoA::degree
    cout << "deg(f, 0) = " << deg(f, 0) << endl; // of type long
    cout << "CmpWDeg(f, g) = " << CmpWDeg(f, g) << endl;
    cout << "CmpWDeg(power(f,2), g) = " << CmpWDeg(power(f,2), g) << endl;
    cout << "CmpWDeg(f, power(g,2)) = " << CmpWDeg(f, power(g,2)) << endl;
    if (!IsZero(f-g))
      cout << "LPP(f-g) = " << LPP(f-g) << "   in " << owner(LPP(f-g)) << endl;
    if (!IsZero(f+g))
      cout << "LC(f+g) = " << LC(f+g) << "   in " << owner(LC(f+g)) << endl;
    cout << "deriv(g, n-1) = " << deriv(g, n-1) << endl; // note second arg!
    if (GradingDim(owner(f)) > 0)
      cout << "IsHomog(f) = " << IsHomog(f) << endl;
    else
      cout << "IsHomog: works only if the grading dimension is positive!" << endl;
    cout << "--------------------------------------------" << endl << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << boolalpha; // so that bools print out as true/false

    cout << ShortDescription << endl;

    // Define some coefficient rings.
    ring QQ   = RingQQ();
    ring FF_7 = NewZZmod(7);

    // Now we show some different ways of creating polynomial rings.

    // [1] If we want, we can specify the names of the indets:
    PolyRing P1 = NewPolyRing(FF_7, symbols("beta,Psi")); // FF_7[beta,Psi]
    SomeComputations(P1);

    // [2] We can a range of names for the indets:
    PolyRing P2 = NewPolyRing(FF_7, SymbolRange("x", 1,3)); // FF_7[x[1],x[2],x[3]]
    SomeComputations(P2);

    // [3] We can let CoCoA choose "new" names for the indets:
    PolyRing P3 = NewPolyRing(QQ, NewSymbols(5));         // indets have (strange) "new" names
    SomeComputations(P3);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PolyRing1.C,v 1.27 2022/02/14 11:23:11 bigatti Exp $
// $Log: ex-PolyRing1.C,v $
// Revision 1.27  2022/02/14 11:23:11  bigatti
// Summary: added boolalpha
//
// Revision 1.26  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.25  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.24  2017/02/08 16:58:32  abbott
// Summary: Minor changes to improve readability
//
// Revision 1.23  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.22  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.21  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.20  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.19  2014/05/08 15:26:55  bigatti
// Summary: deg instead of MaxExponent
//
// Revision 1.18  2014/04/17 13:22:36  bigatti
// -- simplified
//
// Revision 1.17  2014/04/17 08:40:03  bigatti
// -- simplified (focussed on RingElem operations)
//
// Revision 1.16  2013/03/25 17:23:23  abbott
// Added check that GradingDim is positive before testing for homogeneity.
//
// Revision 1.15  2012/07/04 12:17:18  abbott
// Changed size_t into long; several other minor "cleaning up" changes.
//
// Revision 1.14  2012/02/08 17:42:26  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.13  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.12  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.11  2010/03/22 10:17:02  abbott
// Added a comment about BigPrime.
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
