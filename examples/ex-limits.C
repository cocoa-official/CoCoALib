// Copyright (c) 2011  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows the numeric limits for CoCoALib.  \n";

const string LongDescription =
  "The numeric limits of CoCoALib depend on many factors: \n"
  "some choices from the authors of CoCoALib,                  \n"
  "the compilation flags of CoCoALib, which in turn depend on  \n"
  "the architecture of your machine and the compilations flags \n"
  "used by the gmp and boost libraries.                        \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "CoCoALib returns no unsigned types" << endl;
  
    cout << "max of int:  " << std::numeric_limits<int>::max() << endl;
    cout << "digits int:  " << std::numeric_limits<int>::digits << endl;
    cout << "max of long: " << std::numeric_limits<long>::max() << endl;
    cout << "digits long: " << std::numeric_limits<long>::digits << endl;

    double t;
    cout << endl;
    cout << "Here is a test trying to show which is faster (integer division):" << endl;
    cout << "I'm not sure it is reliable:" << endl;
    {
      t=CpuTime();
      long sum = 0;
      for (long i = 0; i<10000000; ++i) sum += (i%1000)/(i);
      cout << "long " << sum << " computed in " << CpuTime()-t << "s" << endl;
    }
    {
      t=CpuTime();
      int sum = 0;
      for (int i = 0; i<10000000; ++i) sum += (i%1000)/(i);
      cout << "int  " << sum << " computed in " << CpuTime()-t << "s" << endl;
    }
    {
      t=CpuTime();
      long sum = 0;
      for (long i = 0; i<10000000; ++i) sum += (i%1000)/(i);
      cout << "long " << sum << " computed in " << CpuTime()-t << "s" << endl;
    }
    {
      t=CpuTime();
      int sum = 0;
      for (int i = 0; i<10000000; ++i) sum += (i%1000)/(i);
      cout << "int  " << sum << " computed in " << CpuTime()-t << "s" << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-limits.C,v 1.3 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-limits.C,v $
// Revision 1.3  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.2  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.1  2011/03/04 17:30:16  bigatti
// -- first import
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
