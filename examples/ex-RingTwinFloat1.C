// Copyright (c) 2005-2009  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <cstdlib>
// using exit

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example showing some features of RingTwinFloat.\n"
  "Program to explore the precision offered by RingTwinFloat\n";
//----------------------------------------------------------------------

#include <algorithm> // using std::swap;

namespace CoCoA
{

  RingElem SquareRoot(RingElem x)
  {
    if (!IsRingTwinFloat(owner(x)))
      CoCoA_THROW_ERROR("Argument must be element of RingTwinFloat", "SquareRoot(RingElem)");
    if (x < 0) CoCoA_THROW_ERROR("Squareroot of negative number", "SquareRoot(RingElem)");
    if (IsZero(x)) return x;
    ring RR = owner(x);
    RingElem approx(RR);
    RingElem approx2(RR, 1);

    //  while (approx != approx2)
    while (!IsPracticallyEqual(approx, approx2))
    {
      approx = (approx2 + x/approx2)/2;
      swap(approx, approx2);
    }
    return (approx2 + x/approx2)/2;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Enter the bit precision parameter for RingTwinFloat (positive integer): ";

    long BitPrec;
    cin >> BitPrec;
    if ( !cin || BitPrec < 1)
    {
      cerr << "*** ERROR *** Input must be a positive integer" << endl;
      exit(1);
    }
    const RingTwinFloat RR = NewRingTwinFloat(BitPrec);
    cout << "The ring is " << RR << endl;
    if (PrecisionBits(RR) != BitPrec)
      cout << "NOTE: Automatically increased precision to " << PrecisionBits(RR) << endl << endl;
    const RingElem three(RR, 3);
    RingElem eps(RR, 1);

    cout << endl;

    // First trial...
    cout << "FIRST TRIAL: loop comparing x+3 with (x^2-9)/(x-3) for\n";
    cout << "x = 3 + eps where  eps = 1/(2^n) for  n=1,2,...200" << endl;
    try
    {
      for (size_t i = 1; i < 200; ++i)
      {
        eps /= 2;
        const RingElem x = three + eps;
        const RingElem v1 = x + 3;
        const RingElem v2 = (x*x-9)/(x-3);
        if (v1 != v2) cout << "*****BUG: v1 and v2 DIFFER*****" << endl;
      }
      cout << "SUCCEEDED -- the program correctly said they are equal." << endl;
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      cout << "INSUFFICIENT PRECISION to complete the loop" << endl;
    }
    cout << endl << endl;

    // Second trial...
    cout << "SECOND TRIAL: see whether (sqrt(2)-1)^6 = 99-70*sqrt(2)\n";
    try
    {
      const RingElem sqrt2 = SquareRoot(RingElem(RR,2));
      if (sqrt2*sqrt2 != 2)
        cout << "*****Bad square root of 2*****" << endl;
      if (power(sqrt2-1,6) == 99-70*sqrt2)
        cout << "SUCCEEDED -- the program correctly said they are equal." << endl;
      else
        cout << "*****FAILED***** -- the program thought they were different." << endl;
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      cout << "INSUFFICIENT PRECISION to complete this trial" << endl;
    }
    cout << endl << endl;

    // Third trial...
    cout << "THIRD TRIAL: an almost equality between square roots:\n";
    cout << "sqrt(176)+sqrt(195)+sqrt(2025) =?= sqrt(190)+sqrt(398)+sqrt(1482)" << endl;
    const RingElem sum1 = SquareRoot(RingElem(RR,176)) + SquareRoot(RingElem(RR,195)) + SquareRoot(RingElem(RR,2025));
    const RingElem sum2 = SquareRoot(RingElem(RR,190)) + SquareRoot(RingElem(RR,398)) + SquareRoot(RingElem(RR,1482));
    cout << "sum1 = " << sum1 << endl;
    cout << "sum2 = " << sum2 << endl;
    bool SumsAreEqual;
    try
    {
      // The equality comparison can produce true, false, or trigger InsufficientPrecision.
      SumsAreEqual = (sum1 == sum2);
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      // The comparison triggered an exception; inform the user, and then return.
      cout << "unable to decide with adequate certainty whether" << endl
           << "the two values are equal or different." << endl;
      return;
    }
    // The comparison produced either true or false.
    if (SumsAreEqual)
      cout << "The program says they are equal (but the programmer knows they are unequal)" << endl;
    else
      cout << "SUCCEEDED -- the program correctly says that they are unequal." << endl;
    try
    {
      // The subtraction may trigger InsufficientPrecision even though the comparison above did not.
      const RingElem diff = sum1 - sum2;
      cout << endl << "and sum1 - sum2 = " << diff << endl;
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      if (SumsAreEqual) // SumsAreEqual should be false -- give error if not.
        CoCoA_THROW_ERROR("Values are equal but failed to compute difference", "ex-RingTwinFloat1");
      cout << "The two values are different, but more precision is needed to compute" << endl
           << "the difference with the required accuracy (i.e. " << PrecisionBits(RR) << " bits)." << endl;
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingTwinFloat1.C,v 1.12 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingTwinFloat1.C,v $
// Revision 1.12  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.11  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.10  2016/03/25 19:56:28  abbott
// Summary: Added some consts; improved a couple of comments
//
// Revision 1.9  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.8  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.7  2012/04/04 08:39:21  bigatti
// -- beautified printouts
//
// Revision 1.6  2011/03/14 10:23:20  abbott
// Changed size_t into long.
//
// Revision 1.5  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/01/19 17:36:32  abbott
// Added extra check on the input value.
//
// Revision 1.2  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.1  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/12 15:59:00  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.4  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.2  2006/07/17 19:10:04  cocoa
// Fixed "hanging print" when an exception was thrown.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/02/13 12:08:04  cocoa
// Fixed a problem with some missing assignment ops for certain PPMonoidElems.
// Fixed a bug in RingDistrMPoly::myIndetPower.
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
