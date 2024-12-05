// Copyright (c) 2009  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program illustrating basic use of BigRat values (i.e. rational numbers) \n"
  "showing that they can be used with a natural syntax.                \n";

const string LongDescription =
  "Program giving example of basic arithmetic with exact rational numbers \n"
  "represented as values of type BigRat.  The syntax recalls that used for    \n"
  "the built-in C++ numerical types.  Emphasis is on convenience rather   \n"
  "utmost execution speed.  To understand better the difference between   \n"
  "a BigRat value and an element of the ring RingQ, contrast this example     \n"
  "with ex-RingQ1.C.                                                      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Rational numbers can be constructed from a pair of integers
    // (being numerator and denominator).
    cout << "The number 3 as a rational is  " << BigRat(3) << endl
         << "Its reciprocal is  " << BigRat(1,3) << "  which is the same as  " << 1/BigRat(3) << endl
         << "The fraction is automatically simplified: e.g. BigRat(2,6) = " << BigRat(2,6) << endl
         << "You can also make a rational from a string: " << BigRatFromString("22/7") << endl
         << endl;

    // Rational numbers are always exact; they are not approximated.
    const BigRat OneThird(1,3); // OneThird = 1/3
    cout << "1/3 + 1/3 + 1/3 - 1 = " << OneThird + OneThird + OneThird - 1 << endl;
    cout << "3*(1/3) - 1 = " << 3*OneThird - 1 << endl;

    // Here is an important caveat: be very careful about rational constants.
    // One might reasonably expect  OneThird + 2/3  to produce 1, but it does not!
    // The C++ compiler interprets the expression  2/3  as an integer division.
    cout << "This value is NOT equal to one: " << OneThird + 2/3 << endl;

    // Here is one way to obtain the desired behaviour:
    cout << "This value IS equal to one: " << OneThird + BigRat(2,3) << endl;

    // The functions  num  and  den  give the numerator and denominator (as a BigInt).
    // den(Q) is always positive, and  num and den are always coprime.
    const BigRat q(123,456);
    cout << "num(" << q << ") = " << num(q) << endl;
    cout << "den(" << q << ") = " << den(q) << endl;

    // The usual arithmetic operators work as you would expect, but note that
    // operator% is NOT DEFINED as it does not make sense.
    const BigRat q1 = (2*q+1)/(4*q-3);
    const BigRat q2 = power(q,2);

    // The usual comparison operators work as you would expect.
    // There is also the function  cmp(a,b)  which returns a machine integer
    // which is <0, =0, >0 according as a<b, a=b, a>b.
    if (q1 < q2)   cout << q1 << " is smaller than " << q2 << endl;
    if (q1 == q2)  cout << q1 << " is equal to " << q2 << endl;
    if (q1 > q2)   cout << q1 << " is larger than " << q2 << endl;

    cout << "cmp(" << q1 << ", " << q2 << ") = " << cmp(q1, q2) << endl;

    // There are a few specific tests:
    if (IsZero(q))      cout << q << " is zero" << endl;
    if (IsOne(q))       cout << q << " is one" << endl;
    if (IsMinusOne(q))  cout << q << " is minus one" << endl;

    // Conversion of a rational into an integer.
    // If you want to know whether a rational number is an integer...
    if (IsOneDen(q)) cout << q << " is the integer " << num(q) << endl;
    // Three functions to convert to a nearby integer:
    cout << "floor(" << q << ") = " << floor(q) << endl;
    cout << "ceil(" << q << ") = " << ceil(q) << endl;
    cout << "round(" << q << ") = " << round(q) << endl;
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
