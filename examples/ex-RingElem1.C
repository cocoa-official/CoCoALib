// Copyright (c) 2012,2021  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "First introduction to polynomial rings and polynomials (RingElem) \n";

const string LongDescription =
  "Create your first polynomial ring and polynomial.           \n"
  "Basic operations an printing. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Here is an intro to creating/using RingElems in C++.
  // You should consult a good introductory book.
  // Later, the books by Scott Meyers contain many useful
  // design hints (but they assume you know basic C++ already).

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // First step is to create the ring you want to use:
    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    
    // An easy way to create specific elements is with RingElem ctor
    RingElem f = RingElem(P, "x^2+y^2-3");

    cout << "f  is  " << f << endl;
    cout << "2*f +3  is  " << 2*f +3 << endl;
    cout << "Square of  f  is  " << power(f,2) << endl;
    // Recall: you CANNOT use ^ for exponentiation!!

    cout << "-------------------------------" << endl;
    // Suppose we do not know the ring of f:
    // "owner" gives the ring to which a RingElem belongs
    cout << "f is an element of the ring " << owner(f) << endl;

    // We can check whether a ring is a polynomial ring:
    if (IsPolyRing(owner(f)))
    {
      cout << "f is an element of a polynomial ring" << endl;
      cout << "the indeterminates of this ring are  " << indets(owner(f)) << endl;
      cout << "and its coefficient ring is  " << CoeffRing(owner(f)) << endl;
    }

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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingElem1.C,v 1.16 2021/09/16 11:36:12 bigatti Exp $
// $Log: ex-RingElem1.C,v $
// Revision 1.16  2021/09/16 11:36:12  bigatti
// Summary: polished
//
// Revision 1.15  2021/09/13 13:55:08  abbott
// Summary: Completely new version: better for complete beginners (hopefully)
//
// Revision 1.1  2021/09/03 13:04:58  abbott
// Summary: Basic Example showing how to use rings and ringelems
//
// Revision 1.8  2020/02/14 09:42:37  abbott
// Summary: Moved Anna's example about translating first/last into C++ into new file ex-c++-vector2.C
//
// Revision 1.7  2020/02/11 10:50:56  bigatti
// -- added examples for assigning vectors (RingElems and c++14 {})
//
// Revision 1.6  2020/01/13 17:08:07  bigatti
// -- minor improvements
// -- removed use of CoCoAVector
//
// Revision 1.5  2016/11/18 07:36:56  bigatti
// -- more examples
//
// Revision 1.4  2015/06/29 15:51:44  bigatti
// -- code in namespace CoCoA
//
// Revision 1.3  2013/05/31 15:10:27  bigatti
// -- just a little more...
//
// Revision 1.2  2012/05/11 10:13:21  bigatti
// -- moved semicolon for bug in my html-index code
//
// Revision 1.1  2012/05/11 10:07:57  bigatti
// first import
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
