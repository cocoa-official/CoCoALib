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
