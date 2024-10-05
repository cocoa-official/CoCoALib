// Copyright (c) 2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows simple use of the functions for creating new   \n"
  "anonymous symbols.                                                \n";

const string LongDescription =
  "Often algorithms need one to work with an addition new indeterminate\n"
  "which must be different from all indeterminates already being use.  \n"
  "In CoCoALib we can create new \"anonymous\" indeterminates guaranteed\n"
  "to be different from all others.  This examples shows how to create \n"
  "them singly, and also several all together.                         \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Function to count how many distinct irred factors x^n-1 has in k[x].
  // This function will work regardless of which symbols appear in k.
  long NumFactors(ring k, int n)
  {
    PolyRing P = NewPolyRing(k, NewSymbols(1));
    const RingElem x = indet(P,0);
    const factorization<RingElem> facs = factor(power(x,n)-1);
    return len(facs.myFactors());
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // To create a single new symbol use the function NewSymbol().
    // Each time you call it you will get a different symbol.
    cout << "First new symbol: " << NewSymbol() << endl;
    cout << "Second new symbol: " << NewSymbol() << endl;

    // If you want several new symbols, use the function NewSymbols(n).
    vector<symbol> v = NewSymbols(5);
    cout << "Here are " << len(v) << " newly minted symbols: " << v << endl;

    ring k1 = RingQQ();
    long n1 = NumFactors(k1, 8);
    cout << "x^8-1  has  " << n1 << " factors over " << k1 << endl;
    ring k2 = NewZZmod(257);
    long n2 = NumFactors(k2, 8);
    cout << "x^8-1  has  " << n2 << " factors over " << k2 << endl;

    PolyRing Qx = NewPolyRing(k1, symbols("x"));
    ring k3 = NewQuotientRing(Qx, "x^2+1"); // QQ(i)
    // same as  NewQuotientRing(Qx, ideal(RingElem(Qx, "x^2+1")));
    long n3 = NumFactors(k1, 8);
    cout << "x^8-1  has  " << n3 << " factors over " << k3 << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-symbol2.C,v 1.7 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-symbol2.C,v $
// Revision 1.7  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2019/10/10 16:36:19  bigatti
// -- syntax with string for NewQuotientRing
//
// Revision 1.5  2019/10/03 14:53:19  bigatti
// -- uncommented example: factorization over QQ(i) works.
//
// Revision 1.4  2015/06/29 15:52:45  bigatti
// *** empty log message ***
//
// Revision 1.3  2015/06/29 13:17:29  bigatti
// -- code inside namespace CoCoA
//
// Revision 1.2  2014/03/24 12:09:20  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.1  2012/05/10 14:41:52  abbott
// New example for anonymous symbols.
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
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
