// Copyright (c) 2008  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Creation of symbols, and some simple operations on them.\n";

const string LongDescription =
  "Symbols are used to give print names to indeterminates.  Their main   \n"
  "use is as an argument to a pseudo-ctor for a PPMonoid or a polynomial \n"
  "ring -- see examples for those types.                                 \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Symbols without subscripts" << endl;
    symbol x("x");
    symbol alpha("alpha");
    cout << x << endl;
    cout << alpha << endl << endl;  

    cout << "Symbols with 1 subscript" << endl;
    symbol y2("y",2);
    symbol zeta("zeta",0);
    cout << y2 << endl;
    cout << zeta << endl << endl;  

    cout << "Symbols with 2 subscripts" << endl;
    symbol amn("a",4,7);
    symbol s("Y_phi", -3,8);
    cout << amn << endl;
    cout << s << endl << endl;

    cout <<"Symbols with many subscripts -- flexible but cumbersome!" << endl;
    vector<long> subs = {1,3,5};
    symbol XX("X", subs);
    cout << XX << endl << endl;

    vector<symbol> v;
    v.push_back(alpha);
    v.push_back(XX);
    v.push_back(s);
    cout << "Constructors of polynomial ring requires a vector of symbols\n"
         << "to give PRINT NAMES to indeterminates, for example\n"
         << v << endl << endl;

    cout << "There are convenience functions for making vectors of symbols"
         << endl;
    cout << "- vector of length 1:  " << symbols("omega") << endl;
    cout << "- vector of length 2:  " << symbols("mu, nu") << endl;
    cout << "- vector with range of subscripts:  " 
         << SymbolRange("B", -2,2) << endl;
    cout << "- vector with a \"rectangle\" of double subscripts:\n  "
         << SymbolRange(symbol("y",0,1), symbol("y",2,2)) << endl << endl;

    cout << "Querying a symbol to get its head, number and values of subscripts "
         << endl;
    cout << "s\t is " << s
         << "\t head(s)\t is  " << head(s) << endl;
    cout << "XX\t is "<< XX
         << "\t NumSubscripts(XX)\t is  " << NumSubscripts(XX) << endl;
    cout << "y2\t is " << y2
         << "\t subscript(y2,0)\t is  " << subscript(y2, 0) << endl;

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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-symbol1.C,v 1.10 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-symbol1.C,v $
// Revision 1.10  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.9  2020/01/13 17:09:23  bigatti
// -- using C++11 { } for making a vector
//
// Revision 1.8  2015/07/27 11:50:49  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.7  2015/06/29 15:52:45  bigatti
// *** empty log message ***
//
// Revision 1.6  2015/06/29 13:17:29  bigatti
// -- code inside namespace CoCoA
//
// Revision 1.5  2012/05/30 16:08:45  bigatti
// -- improvements suggested by A.Caleo
//
// Revision 1.4  2012/05/24 14:53:35  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2008/12/12 16:32:09  abbott
// Improved comments.
//
// Revision 1.1  2008/12/12 11:29:48  abbott
// Fixed a bug in SymbolRange.  Added example and test for symbols.
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
