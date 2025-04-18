// Copyright (c) 2017,2020  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some basic use of C++ vectors    \n"
  "(which are a bit like lists in CoCoA-5).                    \n"
  "See also ex-c++-loop-for2.C for iterating over the elements.\n";

const string LongDescription =
  "This is an example showing some basic use of C++ vectors.    \n"
  "C++ has few built-in types, but the STL (Standard Template   \n"
  "Library) contains many useful extensions.  std::vector is    \n"
  "is one of these extensions.  It is a homogeneous array: a    \n"
  "collection of many objects of the same type, and each one    \n"
  "may be accessed directly by an integer index (from 0 to n-1) \n"
  "See also ex-c++-loop-for2.C for iterating over the elements. \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // If you want to translate code from CoCoA-5 to C++ (using
  // features from CoCoALib too :-) then a LIST in CoCoA-5
  // should most likely be converted to a C++ "vector".
  // But do remember that for vectors INDICES START AT 0.
  // This example gives some guidance; see also ex-c++-vector2.C

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // An array of 10 "long" integers; initially all 0.
    vector<long> v(10);

    // Use function len to count how many values are in a vector
    const long n = len(v);  // number of entries in v

    // >>> VECTOR INDICES START AT 0 AND GO TO n-1 <<<
    for (long i=0; i < n; ++i)
    {
      v[i] = i; // index goes inside square brackets.
    }
    cout << "v = " << v << endl; // CoCoALib can print a vector


    // --------------------------------------------
    // C++ offers many useful functions on vectors.
    // We find the function push_back quite useful:
    // it appends values one at a time.
    // Also look up the C++ function "reserve"!

    // An array of "long" integers; initially empty.
    vector<long> w;
    w.reserve(10);  // reserve not necessary, but a good idea!

    // Now we "fill" w by appending values repeatedly:
    for (int i=0; i < 10; ++i)
    {
      w.push_back(i);
    }
    cout << "w = " << w << endl; // print out w
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
