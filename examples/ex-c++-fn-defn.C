// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some basic function definitions in C++. \n"
  "The command `return' ends computation in the function/procedure;   \n"
  "in a function it also says which value to return to the caller.    \n";

const string LongDescription =
  "This is an example showing some basic function definitions in C++. \n"
  "The command `return' ends computation in the function/procedure;   \n"
  "in a function it also says which value to return to the caller.    \n"
  "C++ barely distinguishes between \"functions\" (which return a     \n"
  "value) and \"procedures\" (which return no value).                 \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // function definition syntax:
  //
  // <rtn-type>   <fn-name>(<arg-type> <arg-name>)
  // {
  //   <commands>
  // }


  // a simple FUNCTION:
  long square(long n)
  {
    return n*n;
  }


  // a slightly longer function (with 2 args)
  bool IsInDisc100(long x, long y)
  {
    if (x < -100 || x > 100 || y < -100 || y > 100) return false;
    long SqDist = square(x) + square(y);
    return (SqDist <= square(100));
  }


  // a PROCEDURE has return type "void"
  // Here first two params are REFERENCES:
  // this procedure will change their values.  See the call below!
  void QuotientAndRemainder(long& quot, long& rem, long a, long b)
  {
    quot = a/b; // integer division!
    rem = a%b;  // operator % computes the remainder
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // Call "square"
    long a = 7;
    long b = 3;
    cout << "square of " << a+b << " is " << square(a+b) << endl;

    // Call "QuotientAndRemainder"
    long q = 0; // initial value does not matter
    long r = 0; // (ditto)
    QuotientAndRemainder(q, r, a, b); // NB changes values of q and r!
    cout << "Quotient and remainder for " << a << " and " << b
         << " are  quot=" << q << " and rem=" << r << endl;

    // Call "IsInDisc100"
    cout << "Does the point (" << a << ", " << b << ") lie in disc of radius 100? "
         << IsInDisc100(a,b) << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-fn-defn.C,v 1.4 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-c++-fn-defn.C,v $
// Revision 1.4  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.3  2021/09/03 12:59:37  abbott
// Summary: Improved comments/descr
//
// Revision 1.2  2019/03/27 14:57:37  bigatti
// -- better layout for example of C++ syntax
//
// Revision 1.1  2017/02/15 12:22:41  abbott
// Summary: New C++ examples
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
