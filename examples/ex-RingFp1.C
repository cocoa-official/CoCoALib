// Copyright (c) 2005  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <cstdlib>
// using exit

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Inefficient program to compute sqrt(2) modulo a given number.    \n"
  "Simple example using finite fields or integers modulo N.\n";

const string LongDescription =
  "The program asks the user for the value of N, it creates the     \n"
  "ring of integers mod N, and finally uses \"brute force\" to find \n"
  "all square-roots of 2 modulo N.                                  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Enter a small number: ";
    int p;
    cin >> p;
    if ( !cin )
    {
      cerr << "*** ERROR *** Input must be a positive integer" << endl;
      exit(1);
    }
    if (p == 0)
    {
      cout << "There is no integer square-root of 2 modulo 0." << endl;
      return;
    }

    ring Fp = NewZZmod(p);
    if (!IsField(Fp))
      cout << "The number you entered, p=" << p << ", is not prime:" << endl
           << "so the integers mod p are implemented as a general QuotientRing."
           << endl;
    cout << "Fp is " << Fp << endl;

    // Just blindly try all elements of Fp (except -1)
    for (RingElem x(Fp); !IsZero(x+1); x += 1)
      if (x*x-2 == 0)
        cout << "A square-root of 2 modulo " << p
             << " is " << x << endl;
    cout << "Search finished" << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingFp1.C,v 1.9 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingFp1.C,v $
// Revision 1.9  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2015/07/01 16:31:35  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.7  2015/06/29 15:47:58  bigatti
// -- code in namespace CoCoA
//
// Revision 1.6  2012/02/08 17:43:14  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.5  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2009/10/29 18:25:08  abbott
// Added missing include of cstdlib for exit.
//
// Revision 1.2  2007/03/22 22:45:31  abbott
// Removed spaces at ends of lines.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 17:41:50  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 15:59:00  bigatti
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
// Revision 1.2  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
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
