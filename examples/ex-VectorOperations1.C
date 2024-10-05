// Copyright (c) 2006-2007,2009,2014  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program to demonstrate printing of vectors (and lists). \n";

const string LongDescription =
  "This example shows how to print out a C++ vector of values.          \n"
  "It also shows that you can just as easily print a vector of vectors. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    // sum & product of elements in a vector
    vector<int> v1;
    v1.push_back(3);
    v1.push_back(5);
    v1.push_back(5);
    cout << "Vector v1 is " << v1 << endl;
    cout << "sum(v1) = " << sum(v1) << endl;
    cout << "product(v1) = " << product(v1) << endl;

    vector<int> v2;
    v2.push_back(1);
    v2.push_back(1);
    v2.push_back(3);
    cout << "Vector v2 is " << v2 << endl;

    // You can easily print out even a vector<vector<...>>.
    vector< vector<int> > vv;
    vv.push_back(v1);
    vv.push_back(v2);
    cout << "Vector of vector vv = " << vv << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-VectorOperations1.C,v 1.5 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-VectorOperations1.C,v $
// Revision 1.5  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2015/06/29 15:49:01  bigatti
// *** empty log message ***
//
// Revision 1.3  2015/06/29 13:25:54  bigatti
// -- code in namespace CoCoA
//
// Revision 1.2  2014/08/16 13:22:21  abbott
// Summary: Added some explanatory comments
// Author: JAA
//
// Revision 1.1  2014/07/31 16:01:25  abbott
// Summary: renamed from ex-io.C
// Author: JAA
//
// Revision 1.5  2014/05/15 12:27:22  abbott
// Summary: Major update: removed all refs to GlobalOutput etc; now just does std::vector example
// Author: JAA
//
// Revision 1.4  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2009/10/29 18:25:40  abbott
// Added missing include of cstdlib for exit.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
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
// Revision 1.1  2006/05/22 15:52:16  cocoa
// Added preprocess-disg algorithm to ApproxPts.
// Sundry minor improvements.
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
