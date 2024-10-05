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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-BigInt1.C,v 1.8 2022/02/13 09:56:55 abbott Exp $
// $Log: ex-BigInt1.C,v $
// Revision 1.8  2022/02/13 09:56:55  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.7  2018/04/20 18:51:24  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.6  2017/10/17 14:18:23  abbott
// Summary: Removed GlobalRandomSource (unhelpful here)
//
// Revision 1.5  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.4  2012/12/04 09:56:47  abbott
// Improved a comment (warning about integer division).
//
// Revision 1.3  2012/03/30 10:35:15  abbott
// Tidied layout in the description strings.
//
// Revision 1.2  2012/02/10 17:19:54  abbott
// Changed RingZ into RingZZ.
//
// Revision 1.1  2011/08/26 10:19:48  bigatti
// -- renamed after ZZ->BigInt, QQ->BigRat
//
// Revision 1.8  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.7  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.6  2010/03/22 11:48:35  abbott
// Made the example simpler (& more useful?).  Previous example is now ex-ZZ2.
//
// Revision 1.5  2009/12/09 13:23:52  abbott
// Added better checks on the inputs.
//
// Revision 1.4  2009/07/08 12:27:15  abbott
// Minor cleaning.
//
// Revision 1.3  2009/06/05 12:14:56  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInteger.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
// Revision 1.2  2007/03/22 22:45:31  abbott
// Removed spaces at ends of lines.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 16:20:57  bigatti
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
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/12/09 15:08:42  cocoa
// -- added log info
//
