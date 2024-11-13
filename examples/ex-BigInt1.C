// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing basic use of BigInt values: creation, printing, and  \n"
  "some simple arithmetic.                                              \n";

const string LongDescription =
  "This program illustrates basic operations on BigInt values, showing that   \n"
  "they can be used much like normal C++ ints except that there is a almost   \n"
  "no limit on the magnitude of the values.                                   \n"
  "NB If you need extreme efficiency then use the GMP library directly.       \n"
  "Contrast this example with ex-RingZZ1.                                     \n";
//-----------------------------------------------------------------------------


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    BigInt N1;       // default ctor, initial value is 0.
    BigInt N2(99);   // ctor from a machine integer.
//BigInt N2 = 99;  *** WARNING: this syntax DOES NOT WORK *** (it won't even compile!)
    BigInt N3 = BigIntFromString("12345678901234567890"); // ctor from string or C string

    // Basic arithmetic: the usual syntax works.
    // You can do arithmetic between BigInts and machine integers too.
    N1 = N2 + 1;
    N1 = N2 - 2;
    N1 = 3*N2*N3;
    N1 = N2/4;  // !!integer division!!
    N1 = N2%5;

    // *** WARNING: you cannot use ^ for powers ***
    // Instead use the function power:
    N1 = power(2, 99);
    N1 = power(N2, 99);
    N1 = power(2, N2);
    N1 = power(N2, N2);

    // The usual comparisons work.
    // There is also a function cmp(a,b): result is <,=,> 0 according as a <,=,> b
    cout << "N2 = " << N2 << endl
         << "N3 = " << N3 << endl
         << "cmp(N2,N3) = " << cmp(N2,N3) << endl;

    // Function for generating (pseudo-)random numbers in a given range:
    cout << "RandomBigInt(10,20) = " << RandomBigInt(10,20) << endl;

    // Functions for factorial, binomial coefficients, & fibonacci numbers:
    cout << "factorial(8) = "   << factorial(8) << endl
         << "binomial(10,5) = " << binomial(10,5) << endl
         << "fibonacci(0) = "   << fibonacci(0) << endl
         << "fibonacci(1) = "   << fibonacci(1) << endl;
  }

} // end of namespace CoCoA


// We write main() like this so we can handle uncaught CoCoA errors in
// a sensible way (i.e. by announcing them).
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error" << endl;
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
