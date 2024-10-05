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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PolyInput2.C,v 1.14 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-PolyInput2.C,v $
// Revision 1.14  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.13  2020/10/23 07:40:27  bigatti
// Summary: added example for RingElems with empty string
//
// Revision 1.12  2019/10/18 14:12:41  bigatti
// -- updated after having "RingElems"
//
// Revision 1.11  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.10  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.9  2016/06/22 14:22:29  abbott
// Summary: Corrected condition for checking if ifstream was created OK
//
// Revision 1.8  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.7  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.6  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.5  2014/03/21 15:57:31  bigatti
// -- Ring is now first argument of ReadExpr
//
// Revision 1.4  2014/01/30 16:18:25  bigatti
// -- fixed LongDescription for creating index.html
//
// Revision 1.3  2014/01/30 16:15:19  bigatti
// -- changed file name for input (as provided)
//
// Revision 1.2  2014/01/30 15:45:43  bigatti
// -- removed code (now part of CoCoALib: RingElemInput.[CH])
// -- added description
//
// Revision 1.1  2014/01/30 14:22:16  bigatti
// -- first import
//
// Revision 1.6  2012/11/30 14:32:50  bigatti
// -- added some comments about use of /dev/null
//
// Revision 1.5  2012/02/08 17:41:29  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2009/10/29 18:24:49  abbott
// Added missing include of cstdlib for exit.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/26 15:49:08  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/12 15:31:57  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.3  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2006/08/17 10:05:27  cocoa
// -- changed: SmallExponent_t --> long
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
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
// Revision 1.2  2005/05/04 16:32:03  cocoa
// -- simplified syntax for NewPolyRing
// -- input read both from stdin and file
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from BigInt).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
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
// Revision 1.4  2004/12/09 15:07:33  cocoa
// -- simplified code using monomial instead of PushBack
//
// Revision 1.3  2004/12/09 13:31:30  cocoa
// -- choice of I/O moved into program()
//
