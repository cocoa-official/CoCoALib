// Copyright (c) 2017,2020  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing how to translate CoCoA-5 list  \n"
  "functions APPEND, FIRST & LAST into C++.                  \n";

const string LongDescription =
  "This is an example showing how to translate CoCoA-5 list  \n"
  "functions APPEND, FIRST & LAST into C++.  APPEND translates\n"
  "nicely into push_back. Also first(L) and last(L) have nice\n"
  "translations; BUT first(L,k) and last(L,k) do not have    \n"
  "efficient translations -- there are usually better, but   \n"
  "different techniques in C++.                              \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // First: look at ex-c++-vector1.C

  // If you want to translate code from CoCoA-5 to C++ (using
  // features from CoCoALib too :-) then a LIST in CoCoA-5
  // should most likely bo converted to a C++ "vector".
  // But do remember that for vectors INDEXES START AT 0.
  // This example gives some guidance.

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // ---------------------------------------------------------
    // Translating CoCoA-5 functions APPEND, FIRST, LAST, CONCAT
    // ---------------------------------------------------------

    // Note: direct use of "first" and "last" with 2 parameters does
    // not translate well into C++ (because C++ offers more efficient
    // alternatives in many cases).
    // Nevertheless, it can be helpful to know how to translate them
    // directly... we show this here.

    // Construct a "list" with given elements:
    // L := [3,4,5];
    std::vector<long> L = {3,4,5};
    cout << "L is " << L << endl;

    // append(ref L, 6);
    L.push_back(6);
    cout << "After push_back  L is  " << L << endl;

    // first(L);
    cout << "L.front() is  " << L.front() << endl;

    // last(L);
    cout << "L.back() is  " << L.back() << endl;

    // first(L,2);  <-- This is a bit fiddly.
    cout << "first(L,2) gives  " << vector<long>(L.begin(), L.begin()+2) << endl;

    // last(L,3);  <-- This is even more fiddly.
    const long n = len(L);
    cout << "last(L,3) gives  " << vector<long>(L.begin()+(n-3), L.end()) << endl;

    // CONCAT is trickier in C++
    // This example modifies L by putting new elements at the end
    const std::vector<long> L2 = {6,7,8,9};
    L.insert( L.end(),  L2.begin(),  L2.end() );  // L := concat(L,L2);
    cout << "Concatenated with \"L.insert\"  modifies L  --> " << L << endl;
    // L2 remains unchanged.
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-vector2.C,v 1.5 2021/10/07 14:30:37 abbott Exp $
// $Log: ex-c++-vector2.C,v $
// Revision 1.5  2021/10/07 14:30:37  abbott
// Summary: Minor clarifications
//
// Revision 1.4  2021/09/16 11:40:25  bigatti
// Summary: polished
//
// Revision 1.3  2021/09/03 13:04:08  abbott
// Summary: Added example for "concat"
//
// Revision 1.2  2020/03/06 17:37:08  bigatti
// -- long (instead of int)
//
// Revision 1.1  2020/02/14 09:42:37  abbott
// Summary: Moved Anna's example about translating first/last into C++ into new file ex-c++-vector2.C
//
// Revision 1.1  2017/02/15 12:22:42  abbott
// Summary: New C++ examples
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
