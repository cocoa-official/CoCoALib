// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example show how to call IdealOfPoints.  \n";

const string LongDescription =
  "This short example show how to call IdealOfPoints.                \n"
  "The points are specified as rows of a matrix (over some field k). \n"
  "You need a polynomial ring for the answer: its coeff ring must be \n"
  "the same field k as used in the matrix.                           \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    const int DimOfSpace = 3;
    const int NumPoints = 4;
    const ring& k = RingQQ();
    PolyRing P = NewPolyRing(k, SymbolRange("x", 1,DimOfSpace));
    matrix M = NewDenseMat(k, NumPoints, DimOfSpace);
    // Fill rows of M with the points:
    for (int i=0; i < NumPoints; ++i)
    {
      for (int j=0; j < DimOfSpace; ++j)
        SetEntry(M,i,j, RandomLong(-9,9)); //  random point
    }
    ideal I = IdealOfPoints(P, M);
    cout << "Points are the rows of " << M << endl;
    cout << "Ideal of all polys which vanish at these points:" << endl << I << endl;

    matrix M1 = NewDenseMat(RowMat(M, 0));
    I = IdealOfProjectivePoints(P, M1);
    cout << "Points are the rows of " << M1 << endl;
    cout << "Ideal of all homog polys which vanish at these points:" << endl << I << endl;

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
