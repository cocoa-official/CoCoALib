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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-QuotientBasis.C,v 1.13 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-QuotientBasis.C,v $
// Revision 1.13  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.12  2019/11/14 17:44:50  abbott
// Summary: Reduced degree to make it faster
//
// Revision 1.11  2017/07/08 19:07:02  abbott
// Summary: Removed comment out (dodgy) code for reporting unhandled interrupts
//
// Revision 1.10  2015/11/21 19:15:19  abbott
// Summary: Cleaned a comment; folded a long line
//
// Revision 1.9  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.8  2015/06/29 15:42:11  bigatti
// *** empty log message ***
//
// Revision 1.7  2015/06/29 13:06:53  bigatti
// -- code in namespace CoCoA
//
// Revision 1.6  2012/02/08 17:42:26  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.5  2011/05/20 16:17:35  bigatti
// -- added QuotientBasis, so the example is now simpler
//
// Revision 1.4  2011/03/10 17:15:50  abbott
// Replaced size_t by long.
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.11  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.10  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.9  2007/02/26 15:50:25  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.8  2007/02/12 15:45:15  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.7  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/11/27 13:22:23  cocoa
// -- deteted unused #include file
//
// Revision 1.5  2006/11/21 14:35:13  cocoa
// -- fixed/hacked a problem with bind2nd: what's the correct way to fix it?
//
// Revision 1.4  2006/11/20 15:52:54  cocoa
// Fixed some buglets.
//
// Revision 1.3  2006/11/17 18:15:28  cocoa
// -- added const & to QuotientBasis argument
//
// Revision 1.2  2006/08/17 10:06:47  cocoa
// -- changed: list as argument is now reference (very embarrassing..)
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.4  2006/03/09 17:11:21  cocoa
// -- changed: "insert" instead of "copy"
//
// Revision 1.3  2006/03/09 17:03:19  cocoa
// -- modified using "list" instead of "vector" to use "remove_if" member
//    function on lists
//
// Revision 1.2  2005/11/29 13:01:36  cocoa
// -- added example for passing overloaded function to bind2nd
//
// Revision 1.1  2005/11/16 14:07:25  cocoa
// -- first import
//
