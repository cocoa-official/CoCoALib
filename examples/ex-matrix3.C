// Copyright (c) 2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to solve linear systems (matrix equations). \n";

const string LongDescription =
  "This program shows how to solve linear systems (matrix equations). \n"
  "At the moment linear system solving is only partially implemented. \n"
  "It will work if the matrix is over a field; otherwise it may fail. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring QQ = RingQQ();

    // Here is a normal 3x3 matrix (with entries 1,2,3, 4,5,6, 7,8,9)
    // we shall use it later on.
    matrix M = NewDenseMat(QQ,3,3);
    SetEntry(M,0,0, 1);  SetEntry(M,0,1, 2);  SetEntry(M,0,2, 3);
    SetEntry(M,1,0, 4);  SetEntry(M,1,1, 5);  SetEntry(M,1,2, 6);
    SetEntry(M,2,0, 7);  SetEntry(M,2,1, 8);  SetEntry(M,2,2, 9);

    matrix V = NewDenseMat(QQ,3,1);
    SetEntry(V,0,0, 6);  SetEntry(V,1,0, 15);  SetEntry(V,2,0, 24);
    cout << "A solution to the matrix equation is " << LinSolve(M,V) << endl;

    // If no solution exists, LinSolve returns as answer a 0-by-0 matrix...
    SetEntry(V,0,0, 5);
    const matrix soln = LinSolve(M,V); // this one has no solution...
    if (IsMat0x0(soln))
      cout << "The second matrix equation was correctly identified as unsolvable." << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-matrix3.C,v 1.6 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-matrix3.C,v $
// Revision 1.6  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.4  2012/12/04 19:57:03  abbott
// Corrected check for unsolvable linear system (after change in spec of LinSolve).
//
// Revision 1.3  2012/07/04 12:25:08  abbott
// Improved a printed message.
//
// Revision 1.2  2012/06/19 15:45:10  abbott
// Added a comment.
//
// Revision 1.1  2012/04/27 14:47:12  abbott
// New example for LinSolve
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
