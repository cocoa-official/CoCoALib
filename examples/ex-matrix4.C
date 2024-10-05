// Copyright (c) 2014  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to use MatByRows & MatByCols.  \n";

const string LongDescription =
  "A simple example showing how a C++ vector can be viewed as a matrix.\n"
  "You can even view the same vector in several different ways.        \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    ring R = RingZZ();
    vector<RingElem> v(6);
    for (int i=1; i <= 6; ++i) v[i-1] = i;

    cout << "Vector v = " << v << endl;
    MatrixView M2x3ByRows = MatByRows(2,3, v);
    MatrixView M2x3ByCols = MatByCols(2,3, v);
    MatrixView M1x6ByRows = MatByRows(1,6, v);

    cout << "MatByRows(2,3, v) sees v as " << M2x3ByRows << endl;
    cout << "MatByCols(2,3, v) sees v as " << M2x3ByCols << endl;
    cout << "MatByRows(1,6, v) sees v as " << M1x6ByRows << endl;
    cout << endl << endl;

    // As usual, setting a matrix entry will change the vector v.
    SetEntry(M2x3ByRows,1,1, 0);
    cout << "Setting the (1,1) entry of M2x3ByRows changes v; and this" << endl
         << "change is seen by all views into v:" << endl;
    cout << "The underlying vector is now v = " << v << endl;
    cout << "MatByRows(2,3, v) sees v as " << M2x3ByRows << endl;
    cout << "MatByCols(2,3, v) sees v as " << M2x3ByCols << endl;
    cout << "MatByRows(1,6, v) sees v as " << M1x6ByRows << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-matrix4.C,v 1.3 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-matrix4.C,v $
// Revision 1.3  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.2  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.1  2014/07/17 14:50:19  abbott
// Summary: Example for MatByRows/MatByCols
// Author: JAA
//
// Revision 1.7  2013/05/28 07:07:04  bigatti
// -- added "cout << boolalpha": useful for testing
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
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
