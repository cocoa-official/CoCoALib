// Copyright (c)  2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

namespace CoCoA
{

//----------------------------------------------------------------------
const string ShortDescription =
  "(Advanced) example illustrating the `soft transitions' from false to uncertain \n"
  "and from uncertain to true in the equality test between twin floats (as eps    \n"
  "to zero when comparing 1+eps to 1).                                            \n";

const string LongDescription =
  "Example illustrating the `soft transitions' from false to uncertain      \n"
  "and from uncertain to true in the equality test between twin floats.     \n"
  "For many small values of eps this example computes 1+eps and compares    \n"
  "the to 1 (either by direct comparison, or by subtacting and seeing if    \n"
  "the difference is positive or zero.  Direct comparison is better than    \n"
  "computing the difference (which throws InsufficientPrecision more often).\n";
  
//----------------------------------------------------------------------

// This test

  void trial(long BitPrec)
  {
    constexpr int IterMax = 1000;
    const ring RR = NewRingTwinFloat(BitPrec);

    cout << "RING is " << RR << endl;
    BigRat eps = BigRat(1);
    for (long i = 1; i < 99; ++i)
    {
      eps /= 2;

      // Subprogram A
      vector<int> SubprogA(3);
      for (int j=0; j < IterMax; ++j)
      {
        try
        {
          const RingElem X(RR, 1);
          if (X+eps == X)
            ++SubprogA[0];
          else
            ++SubprogA[2];
        }
        catch (const RingTwinFloat::InsufficientPrecision&)
        {
          ++SubprogA[1];
        }
      }

      // Subprogram B
      vector<int> SubprogB(3);
      for (int j=0; j < IterMax; ++j)
      {
        try
        {
          const RingElem X(RR, 1);
          const RingElem Y(RR, 1+eps);
          if (X == Y)
            ++SubprogB[0];
          else
            ++SubprogB[2];
        }
        catch (const RingTwinFloat::InsufficientPrecision&)
        {
          ++SubprogB[1];
        }
      }

      // Subprogram C
      vector<int> SubprogC(3);
      for (int j=0; j < IterMax; ++j)
      {
        try
        {
          const RingElem X(RR, 1);
          const RingElem Y(RR, 1+eps);
          if (IsZero(X-Y))
            ++SubprogC[0];
          else
            ++SubprogC[2];
        }
        catch (const RingTwinFloat::InsufficientPrecision&)
        {
          ++SubprogC[1];
        }
      }


      cout << "precision=" << BitPrec
           << "  iter=" << i
           << "  SubprogA=" << SubprogA
           << "  SubprogB=" << SubprogB
           << "  SubprogC=" << SubprogC
           << endl;
    }
    cout << "---------------------------------" << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    trial(32);
  }

} // end of namespace CoCoA


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
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
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
