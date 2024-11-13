// Copyright (c) 2005,2007  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a basic example about the creation and use of the ring of integers.\n"
  "It illustrates the CoCoALib \"philosophy\" of first creating a ring, and   \n"
  "then computing with values in that ring.                                   \n"
  "The C++ commands for performing arithmetic on RingElems have a natural     \n"
  "syntax (except we cannot use ^ for powers).                                \n"
  "It warns about \"mixed ring arithmetic\", which is forbidden in CoCoALib.  \n";

const string LongDescription =
  "To calculate with elements of a ring we must first create the   \n"
  "ring, then we can create C++ objects of type RingElem which     \n"
  "belong to the ring -- a RingElem can change its value but not   \n"
  "the ring to which it belongs.                                   \n";
//---------------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring ZZ = RingZZ();
    RingElem n(ZZ);    // n is a RingElem belonging to ZZ; initial value is 0.

    // To find the ring to which a RingElem belongs use the function owner:
    cout << "The initial value of n is " << n
         << "; this is an element of " << owner(n) << "." << endl;

    //----------------------------------------------------------------------
    // We can do arithmetic on RingElems.
    n = 6;
    n = n+4;
    RingElem m = 12*n + n/2; // rhs specifies the ring to which m belongs
    cout << m << " = one hundred and twenty five" << endl;

    //----------------------------------------------------------------------
    // In a ring only exact division is allowed (contrast with BigInt values!)
    n = m/5; // OK because m is exactly divisible by 5.
    //  However n = m/7;  would throw a CoCoA::ErrorInfo exception with code ERR::BadQuot

    //----------------------------------------------------------------------
    // The next block contains a caveat.
    {
      // It is no longer possible to create a "different" copy of RingZZ:
      ring ZZ2 = RingZZ();
      if (ZZ != ZZ2)
        cerr << "THIS SHOULD NEVER BE PRINTED!" << endl;
      // in general even if R and S are canonically isomorphic rings, but not identical
      // C++ objects, they are regarded as being different.
      // As an exception, all copies of RingZZ() *are* identical C++ object.

      RingElem n2(ZZ2); // This is a RingElem belonging to ZZ2.
      // Since n and n2 do belong to the same identical ring we can
      // combine them arithmetically...
      try
      {
        n+n2; // This used to cause an exception...
      }
      catch (const CoCoA::ErrorInfo& err)
      {
        cerr << "THIS SHOULD NEVER BE PRINTED!" << endl;
        if (err != ERR::MixedRings) throw;
        // It was the exception we used to have
      }
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
