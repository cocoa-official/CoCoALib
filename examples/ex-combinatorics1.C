// Copyright (c) 2024  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example illustrates use of the iterators for subsets and for tuples.  \n";

const string LongDescription =
  "This example illustrates use of the iterators for subsets and for tuples.  \n"
  "The iterators actually work on \"indexes\" in the range 0 to n-1.          \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // ---------------------------------
    // Subset iterator example
    cout << "Subsets of {0,1,2}" << endl;
    SubsetIter it(3);  // iterates through all subsets of {0,1,2}
    long counter = 0;
    while (!IsEnded(it))
    {
      if (counter == 4)
        cout << "The 4-th subset is " << *it << endl;
      ++it;
      ++counter;
    }
    cout << "A set with 3 elements has " << counter << " distinct subsets." << endl;

    // ---------------------------------
    // Subset iterator example (fixed cardinality)
    cout << endl << "2-Subsets of {0,1,2,3}" << endl;
    it = SubsetIter(4,2);  // iterates through all 2-subsets of {0,1,2,3}
    counter = 0;
    while (!IsEnded(it))
    {
      if (counter == 4)
        cout << "The 4-th 2-subset is " << *it << endl;
      ++it;
      ++counter;
    }
    cout << "A set with 4 elements has " << counter << " distinct 2-subsets." << endl;

    // ---------------------------------
    // Tuple iterator example:
    cout << endl << "3-tuples of elements from {0,1,2,...,7}" << endl;
    TupleIter iter(8,3); // iterate [0,0,0], [0,0,1], [0,0,2],...,[7,7,7]
    counter = 0;
    while (!IsEnded(iter))
    {
      if (counter == 444)
        cout << "The 444-th tuple is " << *iter << endl;
      ++iter;
      ++counter;
    }
    cout << "There are " << counter << " distinct 3-tuples of elements from {0,1,...,7}." << endl;
    
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

  catch (const CoCoA::InterruptReceived& intr)
  {
    cerr << endl
         << "------------------------------" << endl
         << ">>>  CoCoALib interrupted  <<<" << endl
         << "------------------------------" << endl
         << "-->>  " << intr << "  <<--" << endl;
    return 2;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-combinatorics1.C,v 1.3 2024/03/13 17:28:23 bigatti Exp $
// $Log: ex-combinatorics1.C,v $
// Revision 1.3  2024/03/13 17:28:23  bigatti
// Summary: added example for k-subsets
//
// Revision 1.2  2024/03/08 20:03:00  abbott
// Summary: Added comments
//
// Revision 1.1  2024/02/08 20:46:30  abbott
// Summary: Added new class TupleIter (redmine 379)
//
// Revision 1.16  2021/12/14 08:35:43  abbott
// Summary: Uncommented code for printing out "interrupted" message
//
// Revision 1.15  2020/01/09 13:32:48  abbott
// Summary: Added comment
//
// Revision 1.14  2019/11/14 17:45:59  abbott
// Summary: Added SignalWatcher (in case you want to make your code interruptible)
//
// Revision 1.13  2017/12/01 21:30:10  abbott
// Summary: Added Anna to copyright
//
// Revision 1.12  2017/07/08 19:07:02  abbott
// Summary: Removed comment out (dodgy) code for reporting unhandled interrupts
//
// Revision 1.11  2016/11/18 18:05:15  abbott
// Summary: Added commented out code to catch InterruptReceived
//
// Revision 1.10  2015/06/29 14:23:19  abbott
// Summary: added missing CoCoA:: prefix
// Author: JAA
//
// Revision 1.9  2015/06/29 13:25:54  bigatti
// -- code in namespace CoCoA
//
// Revision 1.8  2015/06/25 14:19:02  abbott
// Summary: Added call to CoCoA::BuildInfo::Printall
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
//
