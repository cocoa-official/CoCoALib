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
