// Copyright (c) 2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to solve linear systems (matrix equations). \n";

const string LongDescription =
  "This program shows how to solve linear systems (matrix equations). \n"
  "At the moment linear system solving is only partially implemented. \n"
  "It will work if the matrix is over a field; otherwise it may fail. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring QQ = RingQQ();

    // Here is a normal 3x3 matrix (with entries 1,2,3, 4,5,6, 7,8,9)
    // we shall use it later on.
    matrix M = NewDenseMat(QQ,3,3);
    SetEntry(M,0,0, 1);  SetEntry(M,0,1, 2);  SetEntry(M,0,2, 3);
    SetEntry(M,1,0, 4);  SetEntry(M,1,1, 5);  SetEntry(M,1,2, 6);
    SetEntry(M,2,0, 7);  SetEntry(M,2,1, 8);  SetEntry(M,2,2, 9);

    matrix V = NewDenseMat(QQ,3,1);
    SetEntry(V,0,0, 6);  SetEntry(V,1,0, 15);  SetEntry(V,2,0, 24);
    cout << "A solution to the matrix equation is " << LinSolve(M,V) << endl;

    // If no solution exists, LinSolve returns as answer a 0-by-0 matrix...
    SetEntry(V,0,0, 5);
    const matrix soln = LinSolve(M,V); // this one has no solution...
    if (IsMat0x0(soln))
      cout << "The second matrix equation was correctly identified as unsolvable." << endl;
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
