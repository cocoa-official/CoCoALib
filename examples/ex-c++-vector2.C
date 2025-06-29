// Copyright (c) 2017,2020  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "C++ vectors behave much like LISTs in CoCoA-5, but faster.\n"
  "This is an example showing how to translate CoCoA-5 list  \n"
  "functions APPEND, CONCAT, FIRST & LAST into C++.          \n";

const string LongDescription =
  "C++ vectors behave much like LISTs in CoCoA-5, but faster. \n"
  "This is an example showing how to translate CoCoA-5 list   \n"
  "functions APPEND, CONCAT FIRST & LAST into C++.            \n"
  "APPEND translates nicely into push_back. Also FIRST(L) and \n"
  "and LAST(L) have simple translations.  However FIRST(L,k)  \n"
  "and LAST(L,k) do not have such simple translations; instead\n"
  "there are usually better, but different techniques in C++. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // First: look at ex-c++-vector1.C

  // If you want to translate code from CoCoA-5 to C++ (using
  // features from CoCoALib too :-) then a LIST in CoCoA-5
  // should most likely bo converted to a C++ "vector".
  // But do remember that for C++ vectors ***INDICES START AT 0***.
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

    // Translate APPEND into the member function push_back:
    // append(ref L, 6);
    L.push_back(6);
    cout << "After push_back  L is  " << L << endl;

    // Translate FIRST into the member function front:
    // first(L);
    cout << "L.front() is  " << L.front() << endl;

    // Translate LAST into the member function back
    // last(L);
    cout << "L.back() is  " << L.back() << endl;

    // This will become easier when CoCoALib jumps to C++20 standard.
    // first(L,2);  <-- This is a bit fiddly.
    cout << "first(L,2) gives  " << vector<long>(L.begin(), L.begin()+2) << endl;

    // This will become easier when CoCoALib jumps to C++20 standard.
    // last(L,3);  <-- This is even more fiddly.
    const long n = len(L);
    cout << "last(L,3) gives  " << vector<long>(L.begin()+(n-3), L.end()) << endl;

    // Translating CONCAT into C++ is cumbersome (but also faster).
    // This example modifies L by putting new elements at the end
    const std::vector<long> L2 = {6,7,8,9};
    // L := concat(L,L2);  becomes the line below
    L.insert( L.end(),  L2.begin(),  L2.end() );
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
