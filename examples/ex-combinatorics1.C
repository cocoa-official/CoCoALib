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
  "The iterators actually work on \"indices\" in the range 0 to n-1.          \n";
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
