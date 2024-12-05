// Copyright (c) 2006-2007,2009,2014  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program to demonstrate printing of vectors (and lists). \n";

const string LongDescription =
  "This example shows how to print out a C++ vector of values.          \n"
  "It also shows that you can just as easily print a vector of vectors. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    // sum & product of elements in a vector
    vector<int> v1;
    v1.push_back(3);
    v1.push_back(5);
    v1.push_back(5);
    cout << "Vector v1 is " << v1 << endl;
    cout << "sum(v1) = " << sum(v1) << endl;
    cout << "product(v1) = " << product(v1) << endl;

    vector<int> v2;
    v2.push_back(1);
    v2.push_back(1);
    v2.push_back(3);
    cout << "Vector v2 is " << v2 << endl;

    // You can easily print out even a vector<vector<...>>.
    vector< vector<int> > vv;
    vv.push_back(v1);
    vv.push_back(v2);
    cout << "Vector of vector vv = " << vv << endl;
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
