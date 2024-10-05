// Copyright (c) 2021  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing the new style \"for\" loop in C++.   \n"
  "See ex-c++-loop-for1.C for \"continue\" and \"break\" commands. \n";

const string LongDescription =
  "This is an example showing the new style \"for\" loop in C++.   \n"
  "See ex-c++-loop-for1.C for \"continue\" and \"break\" commands. \n"
  "We show how to state whether the elements are copied.           \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // SYNTAX of new "for" loop:
  //
  // for  ( <entry-type> <var-name>  :  <list/vector> )
  // {
  //   <commands>
  // }
  // In reality the "for" loop is more flexible than this.

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    vector<int> v{3,1,4,1,5,9,2,6,5};

    // Print out the elements of v:
    // [NB each element of v is copied into digit]
    for (int digit: v)
    {
      cout << digit << "  ";
    }
    cout << endl;

    // Reduce each digit by 1:
    // In C++ "&" means reference to a value
    for (int& digit: v)
    {
      --digit;
    }
    cout << endl;

    // Print out the elements of v:
    // [without copying them]
    for (const int& digit: v)
    {
      cout << digit << "  ";
    }
    cout << endl;
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-loop-for2.C,v 1.2 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-c++-loop-for2.C,v $
// Revision 1.2  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.1  2021/09/03 13:02:02  abbott
// Summary: Added example using new for loop syntax (for containers)
//
// Revision 1.3  2019/03/27 14:57:37  bigatti
// -- better layout for example of C++ syntax
//
// Revision 1.2  2017/02/15 12:22:24  abbott
// Summary: Added countdown loop; several minor changes
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
