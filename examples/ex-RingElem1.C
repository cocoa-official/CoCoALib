// Copyright (c) 2012,2021  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "First introduction to polynomial rings and polynomials (RingElem) \n";

const string LongDescription =
  "Create your first polynomial ring and polynomial.           \n"
  "Basic operations an printing. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Here is an intro to creating/using RingElems in C++.
  // You should consult a good introductory book.
  // Later, the books by Scott Meyers contain many useful
  // design hints (but they assume you know basic C++ already).

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // First step is to create the ring you want to use:
    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    
    // An easy way to create specific elements is with RingElem ctor
    RingElem f = RingElem(P, "x^2+y^2-3");

    cout << "f  is  " << f << endl;
    cout << "2*f +3  is  " << 2*f +3 << endl;
    cout << "Square of  f  is  " << power(f,2) << endl;
    // Recall: you CANNOT use ^ for exponentiation!!

    cout << "-------------------------------" << endl;
    // Suppose we do not know the ring of f:
    // "owner" gives the ring to which a RingElem belongs
    cout << "f is an element of the ring " << owner(f) << endl;

    // We can check whether a ring is a polynomial ring:
    if (IsPolyRing(owner(f)))
    {
      cout << "f is an element of a polynomial ring" << endl;
      cout << "the indeterminates of this ring are  " << indets(owner(f)) << endl;
      cout << "and its coefficient ring is  " << CoeffRing(owner(f)) << endl;
    }

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
