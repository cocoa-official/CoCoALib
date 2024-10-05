// Copyright (c) 2005,2007  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example of use of power products and PPMonoids.     \n"
  "Program exhibiting most functions on power products.\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x",0,3), lex);
    // Indets will be printed are x[0], x[1], x[2] & x[3].

    // We can ask the PPM to create powers of its own indets...
    PPMonoidElem t1 = indet(PPM,0) * IndetPower(PPM,1,3) * IndetPower(PPM,2,7);

    // or, for handy access to the indeterminates in PPM...
    const vector<PPMonoidElem>& x = indets(PPM);

    PPMonoidElem t2 = x[0] * power(x[1],7) * power(x[2],3);

    cout << "initial power products are\nt1 = " << t1 << "\nt2 = " << t2 << endl;
    cout << endl;

    cout << "t1 == t2 gives " << (t1 == t2) << endl;
    cout << "t1 != t2 gives " << (t1 != t2) << endl;
    cout << "t1 < t2  gives " << (t1 < t2) << endl;
    cout << "t1 <= t2 gives " << (t1 <= t2) << endl;
    cout << "t1 > t2  gives " << (t1 > t2) << endl;
    cout << "t1 >= t2 gives " << (t1 >= t2) << endl;
    cout << endl;

    cout << "t1 * t2  gives " << t1*t2 << endl;
    cout << "IsDivisible(t1, t2) gives " << IsDivisible(t1, t2) << endl;
    cout << "We CANNOT compute t2 / t1" << endl;
    cout << endl;

    cout << "colon(t1, t2) gives " << colon(t1, t2) << endl;
    cout << "colon(t2, t1) gives " << colon(t2, t1) << endl;
    cout << "gcd(t1, t2)   gives " << gcd(t1, t2) << endl;
    cout << "lcm(t1, t2)   gives " << lcm(t1, t2) << endl;
    cout << "power(t1, 5)  gives " << power(t1, 5) << endl;
    cout << "IsCoprime(t1, t2) gives " << IsCoprime(t1, t2) << endl;
    cout << endl;

    cout << "StdDeg(t1) gives " << StdDeg(t1) << endl;
    cout << "wdeg(t1) gives " << wdeg(t1) << endl;
    cout << "[note: Lex is automatically ungraded (i.e. graded over Z^0)]" << endl;
    cout << "[see: ex-OrderingGrading1.C]" << endl;
    cout << endl;

    cout << "We can find the power to which each indeterminate appears:" << endl;
    cout << "exponent(t1, 0) gives " << exponent(t1, 0) << endl;
    cout << "exponent(t1, 1) gives " << exponent(t1, 1) << endl;
    cout << "exponent(t1, 2) gives " << exponent(t1, 2) << endl;
    cout << "exponent(t1, 3) gives " << exponent(t1, 3) << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PPMonoidElem1.C,v 1.7 2022/02/13 09:56:57 abbott Exp $
// $Log: ex-PPMonoidElem1.C,v $
// Revision 1.7  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2015/06/29 15:32:38  bigatti
// *** empty log message ***
//
// Revision 1.5  2015/06/29 12:38:13  bigatti
// -- code in namespace cocoa
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/03/18 18:05:05  bigatti
// -- improved syntax for PPOrdering argument
//
// Revision 1.2  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 22:26:27  cocoa
// Removed try..catch constructs from some "simple" examples.
//
// Revision 1.7  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.6  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/28 14:00:13  bigatti
// -- minor: just a comment
//
// Revision 1.3  2007/02/12 15:31:57  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/12/09 15:08:42  cocoa
// -- added log info
//
