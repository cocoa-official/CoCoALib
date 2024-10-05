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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-BigRat1.C,v 1.5 2022/02/13 09:56:56 abbott Exp $
// $Log: ex-BigRat1.C,v $
// Revision 1.5  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2021/01/15 16:59:33  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.3  2018/04/20 18:51:24  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.2  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.1  2011/08/26 10:19:48  bigatti
// -- renamed after ZZ->BigInt, QQ->BigRat
//
// Revision 1.6  2011/08/24 10:46:32  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.5  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.4  2011/06/23 16:01:07  abbott
// Removed single arg ctor QQ(MachineInteger), & consequential changes.
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2010/03/22 11:47:33  abbott
// Added example using QQ ctor from a string; plus a warning about QQ(1/3).
//
// Revision 1.1  2009/07/08 12:26:53  abbott
// Added floor and ceil functions for QQs.
// Added example program for QQs.
// Minor correction to convert.C; minor cleaning to ex-BigInt1.C
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
