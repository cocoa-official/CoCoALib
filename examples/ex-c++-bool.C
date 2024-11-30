// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing basic use of booleans in C++.  \n"
  "See also ex-c++-basic.C.                                  \n";

const string LongDescription =
  "This is an example showing basic use of booleans in C++.  \n"
  "Note especially the `smart' operators for `and' and `or'. \n"
  "See also ex-c++-basic.C.                                  \n";


//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;
    
    bool flag1 = true;  // boolean value: either "true" or "false"
    bool flag2 = false; // boolean value: either "true" or "false"
    
    long a = 1;
    long b = 2;
    // Careful!  TEST EQUALITY with DOUBLE EQUALS!
    // operator == tests whether a is equal to b
    cout << "Are a and b equal?  " <<  (a == b)  << endl;

    // assignment of booleans
    flag1 = (a == b);

    // NOT operator (prefix exclamation mark)
    // AND operator (infix  &&)
    flag2 = !(a<0 && b>0);

    cout << "flag1 = " << flag1 << "   and   flag2 = " << flag2 << endl;
    // This command makes bools print out as "true" or "false"
    // (by default they are printed as 1 and 0 respectively).
    cout << boolalpha;
    cout << "  After \"boolalpha\" are printed in this way:" << endl;
    cout << "flag1 = " << flag1 << "   and flag2 = " << flag2 << endl;

    // SMART AND operator -- check LHS first, check RHS only if necessary
    if (a > 2 && b > 99)  cout << "MISTAKE!" << endl;
    // The first condition is false, so the computer will
    // not check the second condition (as there is no need).


    // SMART OR operator -- check LHS first, check RHS only if necessary
    if (a > 0 || b > 0) cout << "a or b (or both) is positive." << endl;
    // The first condition is true, so the computer will
    // not check the second condition (as there is no need).
  }

} // end of namespace CoCoA

// IGNORE THE STUFF BELOW (at least for now)

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
