// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows some operations on ideals. \n";

const string LongDescription =
  "This short example shows how to perform various operations  \n"
  "involving ideals.  Certain operations are valid only for    \n"
  "ideals in a polynomial ring.                                \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void IdealOperations(const ideal& I)
  {
    cout << endl << "--------------------------------------------" << endl;
    cout << "Information about " << I << endl;

    // ring and generators
    const ring& R = RingOf(I);
    cout << "The ideal is in the ring  " << R << endl;
    cout << "How many generators does I have?  " << NumGens(I) << endl;
    cout << "The generators of I are  " << gens(I) << endl;

    // Various properties
    cout << "Is the ideal zero?  " << IsZero(I) << endl;
    cout << "Is the ideal one (the whole ring)?  " << IsOne(I) << endl;
//NYI    cout << "Is the ideal prime: " << IsPrime(I) << endl;
//NYI    cout << "Is the ideal maximal: " << IsMaximal(I) << endl;

    // Testing membership/containment
    cout << "Is 3 in I?  " << IsElem(3, I) << endl;
    cout << "Is I contained in ideal(0)?  " << IsContained(I, ideal(zero(R))) << endl;

    if (!IsPolyRing(R)) return;
    // Specific operations for ideals in a polynomial ring:
    cout << endl;
    cout << "Groebner basis for I is  " << GBasis(I) << endl;
    cout << "`leading term ideal' of I is  " << LT(I) << endl;
    cout << "Is I zero-dimensional?  " << IsZeroDim(I) << endl;
    cout << "Is I homogeneous?  " << IsHomog(I) << endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // First we create a PolyRing in which to compute -- see also ex-PolyRing1
    PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    // Two polynomials...
    RingElem f = RingElem(P, "x*y-1");
    RingElem g = RingElem(P, "x^2+y^2-1");

    // Now create some ideals
    ideal I1 = ideal(f);
    IdealOperations(I1);
    
    ideal I2 = ideal(f,g);
    IdealOperations(I2);

    IdealOperations(ideal(one(P)));
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
