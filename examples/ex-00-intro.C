// Copyright (c) 2022  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows some very basic features of CoCoALib.  See also \n"
  "the examples ex-c++-XXX for some basic C++ guidelines.             \n";

const string LongDescription =
  "This example shows some very basic features of CoCoALib.           \n"
  "We show types for representing \"big\" integers and rationals.     \n"
  "We also show how to do some computations with polynomials.         \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // --------------------------------------------
  // INTRODUCTION TO COCOALIB EXAMPLE PROGRAMS
  // --------------------------------------------

  // The active part of an example is the procedure called "program";
  // some examples define other functions/procedures as well (but
  // this example does not).
  
  // YOU may safely ignore the function "main" at the bottom;
  // but the COMPILER needs it!   "main" is the same in all examples.

  // Now look inside "program" below: read the comments & the code.

  void program()
  {
    GlobalManager CoCoAFoundations;  // Necessary!

    cout << ShortDescription << endl;

    // Machine integers have only limited ranges, but computations are fast.
    // CoCoALib offers the type "BigInt" which can represent large values.

    BigInt N = 1 + factorial(1000);

    // To see examples of some operations on BigInts look at the examples:
    //   ex-BigInt*        [basic operations]
    //   ex-NumTheory*     [more advanced operations from number theory]
    //   ex-ToString*      [conversion to a string]


    // --------------------------------------------
    // CoCoALib can also compute directly with rational numbers.
    // These values have type "BigRat".
    // ***WARNING***  do not write rational constants (C++ misinterprets them).
    
    BigRat q(1,3);  // It is ***WRONG*** to write: BigRat q = 1/3;

    // To see examples of operations on BigRats look at the examples:
    //   ex-BigRat*        [basic operations]
    //   ex-NumTheory4.C   [continued fractions from number theory]
    //   ex-ToString*      [conversion to a string]
    //
    // Hint: prefer to use BigInt over BigRat: BigRat arithmetic is slow.


    // --------------------------------------------
    // CoCoALib can compute with polynomials.
    // There is no type "poly" or "polynomial"; instead
    // CoCoALib uses "RingElem" (meaning "element of a ring").

    // To compute with polynomials you must first do some "set up":
    // (1) create the ring of coefficients, k  (e.g. just RinQQ())
    // (2) create the polynomial ring, p = k[x,y,z] specifying both
    //     the coefficient ring and the indeterminates;
    // (3) now create your polynomials as "RingElem" belonging to P.

    ring P = NewPolyRing(RingQQ(), symbols("x,y"));  // poly ring  QQ[x,y]
    RingElem f = RingElem(P, "x^105 - y^105");
    cout << "Factorization of f is " << factor(f) << endl;


    // Most of the examples involve polynomials in one way or another,
    // but often only as minor players while exhibiting some other feature.
    // Start by looking at these examples:
    //   ex-ring1.C
    //   ex-RingQQ1.C
    //   ex-RingElem1.C
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

  catch (const CoCoA::InterruptReceived& intr)
  {
    cerr << endl
         << "------------------------------" << endl
         << ">>>  CoCoALib interrupted  <<<" << endl
         << "------------------------------" << endl
         << "-->>  " << intr << "  <<--" << endl;
    return 2;
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
