// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing how some RingTwinFloat values can have more precision   \n"
  "than that requested.                                                    \n";

const string LongDescription =
  "Example showing that certain RingTwinFloat values may have a precision  \n"
  "higher than that requested.  In this case we request 64 bits (i.e.      \n"
  "about 19 decimal digits), but in fact we can remove the nineteen most   \n"
  "significant digits and still have a result with the requested precision.\n"
  "So the original value of the variable third did in fact have at         \n"
  "least 127 bits correct (i.e. about 38 decimal digits).                   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring RR = NewRingTwinFloat(64); // Request 64 bits, about 19 decimal digits of precision.
    RingElem third = one(RR)/3;
    for (size_t i=0; i < 19; ++i)
    {
      cout << "Variable third=" << third << endl;
      cout << "Removing a leading significant digit..." << endl;
      third = 10*third-3;
    }
    cout << endl
         << "We conclude that the variable `third' initially had at least 93 bits of precision." << endl;
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
