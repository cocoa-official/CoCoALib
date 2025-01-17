// Copyright (c) 2014  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to use MatByRows & MatByCols.  \n";

const string LongDescription =
  "A simple example showing how a C++ vector can be viewed as a matrix.\n"
  "You can even view the same vector in several different ways.        \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    ring R = RingZZ();
    vector<RingElem> v(6);
    for (int i=1; i <= 6; ++i) v[i-1] = i;

    cout << "Vector v = " << v << endl;
    MatrixView M2x3ByRows = MatByRows(2,3, v);
    MatrixView M2x3ByCols = MatByCols(2,3, v);
    MatrixView M1x6ByRows = MatByRows(1,6, v);

    cout << "MatByRows(2,3, v) sees v as " << M2x3ByRows << endl;
    cout << "MatByCols(2,3, v) sees v as " << M2x3ByCols << endl;
    cout << "MatByRows(1,6, v) sees v as " << M1x6ByRows << endl;
    cout << endl << endl;

    // As usual, setting a matrix entry will change the vector v.
    SetEntry(M2x3ByRows,1,1, 0);
    cout << "Setting the (1,1) entry of M2x3ByRows changes v; and this" << endl
         << "change is seen by all views into v:" << endl;
    cout << "The underlying vector is now v = " << v << endl;
    cout << "MatByRows(2,3, v) sees v as " << M2x3ByRows << endl;
    cout << "MatByCols(2,3, v) sees v as " << M2x3ByCols << endl;
    cout << "MatByRows(1,6, v) sees v as " << M1x6ByRows << endl;
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
