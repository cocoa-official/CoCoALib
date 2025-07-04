// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program exhibiting a way of using ever higher precisions...            \n"
  "This example shows how failure can be a success: this pathological     \n"
  "computation produces the **same wrong result** when using normal       \n"
  "floating point arithmetic at any given finite precision!  However,     \n"
  "since twin floats are self-checking, we detect that there is a problem.\n";

const string LongDescription =
  "Example showing iterative increase of precision using RingTwinFloat    \n"
  "until the answer is found (or some maximum precision is reached).      \n"
  "This program will always fail to find the limit: J-M Muller's sequence \n"
  "actually converges to 6 (rather slowly), however it is pathological    \n"
  "because it converges to 100 using any finite precision arithmetic.     \n"
  "RingTwinFloat detects the onset of pathological convergence to 100, and\n"
  "throws an InsufficientPrecision exception.                             \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // This function computes J-M Muller's sequence in the ring R until two
  // successive values are "practically equal" (see documentation).
  // It returns the value the sequence converged to.
  // Note: we use the initial values 2 and -4 which can be represented exactly;
  // the original sequence had rational initial values.
  RingElem MullerSeq(const RingTwinFloat& RR)
  {
    RingElem Vprev2(RR, 2);  // These starting values can be represented exactly.
    RingElem Vprev1(RR, -4);

    RingElem Vcurr(RR);
    while (!IsPracticallyEqual(Vprev1, Vprev2))
    {
      Vcurr = 111 - 1130/Vprev1 + 3000/(Vprev1*Vprev2); // Muller's recurrence.
      Vprev2 = Vprev1;
      Vprev1 = Vcurr;
    }
    return Vprev1;
  }


  // This loop tries ever higher precisions until the answer is found.
  // Arbitrarily limit precision to 4096 so it cannot run forever.
  RingElem MullerLoop()
  {
    for (int BitPrec = 32; BitPrec <= 4096; BitPrec *= 2)
    {
      try { return MullerSeq(NewRingTwinFloat(BitPrec)); }
      catch (const RingTwinFloat::InsufficientPrecision&)
      {
        // Inform user about the failure...
        cout << "A bit precision of " << BitPrec << " was not sufficient." << endl;
      }
    }
    CoCoA_THROW_ERROR1("Required too much precision... giving up!");
    return zero(RingZZ()); // Never executed; just to keep compiler quiet.
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    try
    {
      const RingElem limit = MullerLoop();
      cout << "Muller's sequence converges to " << limit << endl;
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      const string& FullMesg = message(err);
      const size_t end = FullMesg.size();
      if (FullMesg.substr(end-41, end) != "Required too much precision... giving up!")  throw;
      cout << endl
           << "As expected we caught this error: " << message(err) << endl;
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
